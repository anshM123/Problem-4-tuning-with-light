"""Light-driven, mode- and momentum-resolved phonon dynamics of any material from its structure.

    python run_material.py <structure file | common name> --pump-eV 1.55 --fluence-mJcm2 5
                           --fwhm-fs 50 --T 300 [--thickness-nm 40] [--tmax-ps 50] ...

Pipeline (every material is treated identically; no per-material fitting):
  1. structure: ASE-readable file or built-in name; validation; MACE (medium-omat-0) relaxation of
     the internal coordinates at fixed cell (symmetry kept).
  2. phonons: harmonic and cubic force constants from MACE finite displacements (phonopy/phono3py).
  3. three-phonon collision operator on the full q-grid (no symmetry reduction), absolute prefactor
     calibrated against phono3py's imaginary self-energy (as materials/graphite/phph.py).
  4. electron-phonon coupling from a universal sp3 Slater-Koster tight-binding model:
     Harrison's universal two-centre hoppings V_ll'm = eta_ll'm hbar^2/(m d^2) (Harrison, Electronic
     Structure and the Properties of Solids, 1980; eta_ss_sigma -1.40, eta_sp_sigma 1.84,
     eta_pp_sigma 3.24, eta_pp_pi -0.81) between nearest-neighbour shells, on-site energies from
     the NIST LDA atomic eigenvalues (Kotochigova et al., NIST SRD 141). The electron-phonon vertex
     is the exact derivative of this Hamiltonian with respect to the MACE phonon displacements
     (d^-2 scaling plus angular terms). The q-resolved, energy-resolved coupling function
     W_qnu(eps) (a mode-resolved Eliashberg function) is tabulated once per material and gives the
     phonon generation coefficients A_qnu(T_e) of materials/graphite/dynamics.py for any electron
     temperature. Semiconductors: the measured or tight-binding gap is imposed by a scissor shift;
     photo-carriers keep their density (no recombination unless --recomb-ps), electrons and holes
     share one temperature with separate quasi-Fermi levels; the energy hv - E_g per pair heats
     the carriers, E_g per pair stays electronic.
  5. optics: absorbed energy from thin-film (transfer-matrix) optics with a complex refractive index
     (user value, built-in table, or the tight-binding dielectric function as fallback).
  6. Level-B kinetics (REPORT.md 2.2): thermalized electrons at T_e + every phonon mode n_qnu(t),
     with the electron-phonon and full nonlinear three-phonon collision integrals; energy conserved.
  7. observables: branch populations and effective temperatures, characteristic times, one-phonon
     diffuse-scattering intensities at high-symmetry points, Debye-Waller (mean-square displacement)
     dynamics of chosen Bragg peaks, and experiment-style fits. Saved as JSON + figure.

Intermediates are cached in materials/blind/<name>/ (re-used when the relevant settings match).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
import warnings

NTHREADS = int(os.environ.get("RM_THREADS", "4"))
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, str(NTHREADS))

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------------------------
# constants
HBAR_eVs = 6.582119569e-16
HBAR_Js = 1.054571817e-34
AMU = 1.66053906660e-27
KB_eV = 8.617333262e-5
THz2eV = 4.135667696e-3            # h * 1 THz (phonopy frequencies are nu in THz)
THz_K = 47.99243                   # h * 1 THz / kB
HB2M = 7.619964                    # hbar^2 / m_e  (eV A^2)
E2_EPS0 = 4 * np.pi * 14.3996454   # e^2 / eps0  (eV A)
EV_J = 1.602176634e-19
HC_eVnm = 1239.841984

# Harrison (1980) universal sp coefficients
ETA = {"sss": -1.40, "sps": 1.84, "pps": 3.24, "ppp": -0.81}

# ----------------------------------------------------------------------------------------------
# built-in structures and tabulated material data (inputs, not fitted quantities)
def _graphite():
    from ase import Atoms
    a = 2.464
    cell = [[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, 15.0]]
    return Atoms("C2", cell=cell, scaled_positions=[[1 / 3, 2 / 3, .5], [2 / 3, 1 / 3, .5]], pbc=True)


def _black_phosphorus():
    from ase.spacegroup import crystal
    # Cmce, Brown & Rundqvist, Acta Cryst. 19, 684 (1965); a zigzag, b out of plane, c armchair
    return crystal("P", [(0, 0.10168, 0.08056)], spacegroup=64, cellpar=[3.3136, 10.478, 4.3763, 90, 90, 90])


def _silicon():
    from ase.build import bulk
    return bulk("Si", "diamond", a=5.431)


def _hbn():
    from ase import Atoms
    a, c = 2.504, 6.66
    cell = [[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]]
    return Atoms("BNBN", cell=cell, pbc=True,
                 scaled_positions=[[1 / 3, 2 / 3, 0], [2 / 3, 1 / 3, 0], [2 / 3, 1 / 3, .5], [1 / 3, 2 / 3, .5]])


BUILTIN = {
    # graphite is modelled, as in materials/graphite, by one graphene layer per 3.35 A of thickness
    "graphite": dict(build=_graphite, layer_thickness=3.35, gap=0.0,
                     nk={1.55: 2.72 + 1.47j}, nk_note="ordinary index of graphite at 800 nm (as materials/graphite/dynamics.py)"),
    "black_phosphorus": dict(build=_black_phosphorus, gap=0.33,
                             gap_note="bulk black phosphorus band gap ~0.3 eV (Keyes, Phys. Rev. 92, 580 (1953); Morita, Appl. Phys. A 39, 227 (1986))"),
    "silicon": dict(build=_silicon, gap=1.12, nk={1.55: 3.68 + 0.0055j},
                    nk_note="Si at 800 nm (Green, Sol. Energy Mater. Sol. Cells 92, 1305 (2008))"),
    "hbn": dict(build=_hbn, gap=6.0),
}
ALIASES = {"bp": "black_phosphorus", "black-phosphorus": "black_phosphorus", "phosphorus": "black_phosphorus",
           "graphene_layer": "graphite", "si": "silicon", "h-bn": "hbn", "boron_nitride": "hbn"}


def log(*a):
    print(*a, flush=True)


def _key(d):
    return hashlib.sha1(json.dumps(d, sort_keys=True, default=str).encode()).hexdigest()[:16]


def _cached(path, key):
    if os.path.exists(path):
        try:
            with np.load(path, allow_pickle=False) as z:
                if "key" in z and str(z["key"]) == key:
                    return True
        except Exception:
            return False
    return False


# ----------------------------------------------------------------------------------------------
# 1. structure
def load_structure(spec):
    from ase.io import read
    meta = {}
    name = spec.lower().replace(" ", "_")
    name = ALIASES.get(name, name)
    if name in BUILTIN:
        atoms = BUILTIN[name]["build"]()
        meta = {k: v for k, v in BUILTIN[name].items() if k != "build"}
        meta["source"] = f"built-in '{name}'"
    elif os.path.exists(spec):
        atoms = read(spec)
        name = os.path.splitext(os.path.basename(spec))[0].lower().replace(" ", "_")
        meta["source"] = os.path.abspath(spec)
    else:
        raise SystemExit(f"error: '{spec}' is neither a readable structure file nor a known name "
                         f"({', '.join(sorted(BUILTIN))})")
    return atoms, name, meta


def periodic_dims(atoms, vac=7.0):
    """True for lattice directions without a vacuum gap > vac (A) between periodic images."""
    cell = np.array(atoms.cell)
    V = abs(np.linalg.det(cell))
    out = []
    for i in range(3):
        j, k = [x for x in range(3) if x != i]
        h = V / np.linalg.norm(np.cross(cell[j], cell[k]))
        f = np.sort(atoms.get_scaled_positions()[:, i] % 1.0)
        gaps = np.diff(np.concatenate([f, [f[0] + 1]]))
        out.append(not (gaps.max() * h > vac))
    return np.array(out)


def validate_structure(atoms):
    from pymatgen.core import Element
    problems, warn = [], []
    if len(atoms) == 0:
        problems.append("no atoms")
    if abs(np.linalg.det(np.array(atoms.cell))) < 1.0:
        problems.append("cell volume is zero or undefined (a periodic cell is required)")
    if len(atoms) > 1:
        from ase.neighborlist import neighbor_list
        d = neighbor_list("d", atoms.copy(), 0.6)
        if len(d):
            problems.append("atoms closer than 0.6 A (overlapping or duplicated sites)")
    for s in set(atoms.get_chemical_symbols()):
        el = Element(s)
        if el.Z > 89:
            problems.append(f"element {s} not covered by the MACE foundation model")
        if el.is_transition_metal or el.is_lanthanoid or el.is_actinoid:
            warn.append(f"{s} has open d/f valence shells: the universal sp tight-binding has no d/f "
                        "orbitals, so electron-phonon coupling and optics are qualitative only")
    if problems:
        raise SystemExit("invalid structure: " + "; ".join(problems))
    return warn


def mace_calc():
    import torch
    torch.set_num_threads(NTHREADS)
    from mace.calculators import mace_mp
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return mace_mp(model="medium-omat-0", default_dtype="float64", device="cpu")


def relax_positions(atoms, calc, fmax=2e-3):
    from ase.constraints import FixSymmetry
    from ase.optimize import BFGS
    at = atoms.copy()
    at.calc = calc
    f0 = np.abs(at.get_forces()).max()
    if f0 < fmax:
        return at, f0, f0
    at.set_constraint(FixSymmetry(at))
    BFGS(at, logfile=None).run(fmax=fmax, steps=300)
    f1 = np.abs(at.get_forces()).max()
    at.set_constraint()
    return at, f0, f1


# ----------------------------------------------------------------------------------------------
# 2. force constants
def supercell_dims(atoms, length, pdims):
    cell = np.array(atoms.cell)
    V = abs(np.linalg.det(cell))
    out = []
    for i in range(3):
        j, k = [x for x in range(3) if x != i]
        h = V / np.linalg.norm(np.cross(cell[j], cell[k]))
        out.append(int(math.ceil(length / h)) if pdims[i] else 1)
    return out


def stage_force_constants(atoms, outdir, calc_getter, args, pdims):
    from phonopy.structure.atoms import PhonopyAtoms
    from phono3py import Phono3py
    from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5
    from ase import Atoms
    from ase.neighborlist import neighbor_list
    sc3 = supercell_dims(atoms, args.fc3_length, pdims)
    sc2 = supercell_dims(atoms, args.fc2_length, pdims)
    dnn = neighbor_list("d", atoms, 6.0).min()
    cutoff = max(4.0, 1.8 * dnn)
    key = _key(dict(cell=np.round(np.array(atoms.cell), 5).tolist(), pos=np.round(atoms.positions, 4).tolist(),
                    sym=atoms.get_chemical_symbols(), sc3=sc3, sc2=sc2, cut=cutoff, model="medium-omat-0"))
    keyfile = os.path.join(outdir, "fc_key.npz")
    yaml = os.path.join(outdir, "phono3py_params.yaml")
    if _cached(keyfile, key) and os.path.exists(yaml):
        log(f"[fc] cached force constants ({outdir})")
        return dict(sc2=sc2, sc3=sc3, cutoff=cutoff, cached=True)
    calc = calc_getter()
    unit = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=np.array(atoms.cell),
                        scaled_positions=atoms.get_scaled_positions())
    ph3 = Phono3py(unit, supercell_matrix=sc3, phonon_supercell_matrix=sc2, primitive_matrix="auto")
    ph3.generate_displacements(distance=0.01, cutoff_pair_distance=cutoff)
    ph3.generate_fc2_displacements(distance=0.01)

    def forces(s):
        at = Atoms(s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True)
        at.calc = calc
        return at.get_forces()

    t0 = time.time()
    ph3.phonon_forces = np.array([forces(s) for s in ph3.phonon_supercells_with_displacements])
    log(f"[fc] fc2: supercell {sc2}, {len(ph3.phonon_supercells_with_displacements)} displacements "
        f"({time.time()-t0:.0f}s)")
    allsc = ph3.supercells_with_displacements
    nreal = sum(s is not None for s in allsc)
    log(f"[fc] fc3: supercell {sc3}, pair cutoff {cutoff:.2f} A, {nreal} displacements")
    F = []
    for i, s in enumerate(allsc):
        F.append(np.zeros((len(ph3.supercell), 3)) if s is None else forces(s))
        if i % 200 == 0:
            log(f"[fc]   {i}/{len(allsc)} ({time.time()-t0:.0f}s)")
    ph3.forces = np.array(F)
    ph3.produce_fc2(symmetrize_fc2=True)
    ph3.produce_fc3(symmetrize_fc3r=True)
    ph3.save(yaml)
    write_fc2_to_hdf5(ph3.fc2, filename=os.path.join(outdir, "fc2.hdf5"))
    write_fc3_to_hdf5(ph3.fc3, filename=os.path.join(outdir, "fc3.hdf5"))
    np.savez(keyfile, key=key)
    log(f"[fc] done ({time.time()-t0:.0f}s)")
    return dict(sc2=sc2, sc3=sc3, cutoff=cutoff, cached=False, n_fc3_disp=int(nreal), seconds=time.time() - t0)


def load_ph3(outdir):
    import phono3py
    return phono3py.load(os.path.join(outdir, "phono3py_params.yaml"),
                         fc2_filename=os.path.join(outdir, "fc2.hdf5"),
                         fc3_filename=os.path.join(outdir, "fc3.hdf5"),
                         produce_fc=False, is_mesh_symmetry=False, log_level=0)


def make_phonopy(ph3):
    from phonopy import Phonopy
    ph = Phonopy(ph3.unitcell, supercell_matrix=ph3.phonon_supercell_matrix, primitive_matrix=ph3.primitive_matrix)
    ph.force_constants = ph3.fc2
    return ph


def phonons_at(ph, qfrac):
    ph.run_qpoints(np.asarray(qfrac, float), with_eigenvectors=True)
    return ph.qpoints.frequencies.copy(), ph.qpoints.eigenvectors.copy()


def default_mesh(prim_cell, pdims, nb, budget=2.5e8, cap=600):
    rec = np.linalg.norm(np.linalg.inv(prim_cell).T, axis=1)  # |b_i| / 2pi
    target = min(cap, math.sqrt(budget / nb ** 3))
    d = int(pdims.sum())
    r = np.where(pdims, rec, 0.0)
    scale = (target / np.prod(r[pdims])) ** (1.0 / d)
    mesh = [max(2, int(round(scale * r[i]))) if pdims[i] else 1 for i in range(3)]
    # keep equal divisions for equal reciprocal lengths (needed for symmetric planes)
    for i in range(3):
        for j in range(i):
            if pdims[i] and pdims[j] and abs(rec[i] - rec[j]) < 1e-3 * rec[i]:
                mesh[i] = mesh[j]
    mesh = [m + (m % 2) if p and m > 2 else m for m, p in zip(mesh, pdims)]  # even: contains zone-boundary points
    return mesh


# ----------------------------------------------------------------------------------------------
# 3. three-phonon collision operator (general version of materials/graphite/phph.py)
def stage_phph(outdir, mesh, sig_dyn, args):
    from phono3py.phonon3.imag_self_energy import ImagSelfEnergy
    path = os.path.join(outdir, "phph_procs.npz")
    key = _key(dict(fc=str(np.load(os.path.join(outdir, "fc_key.npz"))["key"]), mesh=mesh, sig=sig_dyn, fmin=args.fmin_thz))
    if _cached(path, key):
        log("[phph] cached three-phonon processes")
        return dict(np.load(path))
    ph3 = load_ph3(outdir)
    ph3.mesh_numbers = mesh
    ph3.init_phph_interaction()
    itr = ph3.phph_interaction
    grid = ph3.grid
    itr.run_phonon_solver()
    freqs_bz, _, _ = itr.get_phonons()
    N = int(np.prod(mesh))
    g2b = grid.grg2bzg; b2g = grid.bzg2grg
    freq = freqs_bz[g2b[:N]].copy()
    qfrac = grid.addresses[g2b[:N]] / np.array(grid.D_diag, float)
    nb = freq.shape[1]
    valid = freq > args.fmin_thz
    SIGC = 0.25

    def gauss(x, s):
        return np.exp(-0.5 * (x / s) ** 2) / (np.sqrt(2 * np.pi) * s)

    t0 = time.time()
    procs = [[[], [], [], []] for _ in range(3)]
    cal = []
    cal_pts = sorted(set([N // 7 + 3, N // 2 + 1, (2 * N) // 3 + 5]))
    for g in range(N):
        itr.set_grid_point(g2b[g]); itr.run()
        trips, wts, _, _ = itr.get_triplets_at_q()
        q12 = b2g[trips[:, 1:]]
        P = itr.interaction_strength
        q1, q2 = q12[:, 0], q12[:, 1]
        w0 = freq[g][None, :, None, None]; w1 = freq[q1][:, None, :, None]; w2 = freq[q2][:, None, None, :]
        okv = valid[g][None, :, None, None] & valid[q1][:, None, :, None] & valid[q2][:, None, None, :]
        for kind in range(3):
            dlt = (w0 - w1 - w2) if kind == 0 else ((w0 + w1 - w2) if kind == 1 else (w0 - w1 + w2))
            ok = (np.abs(dlt) < 4 * sig_dyn) & okv
            t_, j_, j1_, j2_ = np.nonzero(ok)
            val = wts[t_] * P[t_, j_, j1_, j2_] * gauss(dlt[t_, j_, j1_, j2_], sig_dyn)
            k = val > 0
            pr = procs[kind]
            pr[0].append((g * nb + j_[k]).astype(np.int32)); pr[1].append((q1[t_[k]] * nb + j1_[k]).astype(np.int32))
            pr[2].append((q2[t_[k]] * nb + j2_[k]).astype(np.int32)); pr[3].append(val[k].astype(np.float32))
        if g in cal_pts:
            ise = ImagSelfEnergy(itr)
            ise.set_grid_point(g2b[g]); ise.set_sigma(SIGC); ise.set_temperature(300.0)
            ise.run_interaction(); ise.run_integration_weights(); ise.run()
            gam = np.array(ise.imag_self_energy)
            n1 = nB_thz(freq[q1], 300.0); n2 = nB_thz(freq[q2], 300.0)
            mine = np.zeros(nb)
            for j in range(nb):
                d1 = gauss(freq[g, j] - freq[q1][:, :, None] - freq[q2][:, None, :], SIGC)
                d2 = gauss(freq[g, j] + freq[q1][:, :, None] - freq[q2][:, None, :], SIGC) - \
                    gauss(freq[g, j] - freq[q1][:, :, None] + freq[q2][:, None, :], SIGC)
                term = (n1[:, :, None] + n2[:, None, :] + 1) * d1 + (n1[:, :, None] - n2[:, None, :]) * d2
                mine[j] = np.sum(wts[:, None, None] * P[:, j] * term)
            sel = (gam > 1e-6) & (mine > 0) & valid[g]
            cal += list(gam[sel] / mine[sel])
        if g % max(1, N // 8) == 0:
            log(f"[phph] grid point {g}/{N} ntrip={len(wts)} ({time.time()-t0:.0f}s)")
    cal = np.array(cal)
    C = 4 * np.pi * cal.mean()
    log(f"[phph] prefactor calibration vs phono3py: spread std/mean = {cal.std()/cal.mean():.2e} over {len(cal)} modes")
    out = dict(key=key, freq=freq, qfrac=qfrac, mesh=np.array(mesh), C=C, cal_spread=cal.std() / cal.mean(), sig=sig_dyn)
    for kind in range(3):
        a, b, c, v = [np.concatenate(x) for x in procs[kind]]
        out[f"a{kind}"], out[f"b{kind}"], out[f"c{kind}"], out[f"v{kind}"] = a, b, c, (C * v).astype(np.float32)
        log(f"[phph] kind {kind}: {len(v)} processes")
    np.savez(path, **out)
    log(f"[phph] done ({time.time()-t0:.0f}s)")
    return out


def nB_thz(f, T):
    f = np.asarray(f, float)
    x = f * THz_K / T
    return np.where(f > 1e-3, 1 / np.expm1(np.maximum(x, 1e-12)), 0.0)


# ----------------------------------------------------------------------------------------------
# 4. universal tight-binding model and electron-phonon coupling
def _sk_block(d, oi, oj, Vfun):
    """Slater-Koster block for orbital lists oi, oj ('s','x','y','z') and bond vector d (i -> j)."""
    r = np.linalg.norm(d)
    n = d / r
    V = {k: Vfun(k, r) for k in ("sss", "sps", "pps", "ppp")}
    ax = {"x": 0, "y": 1, "z": 2}
    B = np.zeros((len(oi), len(oj)))
    for a, o1 in enumerate(oi):
        for b, o2 in enumerate(oj):
            if o1 == "s" and o2 == "s":
                B[a, b] = V["sss"]
            elif o1 == "s":
                B[a, b] = n[ax[o2]] * V["sps"]
            elif o2 == "s":
                B[a, b] = -n[ax[o1]] * V["sps"]
            else:
                l1, l2 = n[ax[o1]], n[ax[o2]]
                B[a, b] = l1 * l2 * (V["pps"] - V["ppp"]) + (V["ppp"] if o1 == o2 else 0.0)
    return B


def harrison_V(kind, r):
    return ETA[kind] * HB2M / r ** 2


def atomic_levels(sym):
    """Valence s/p on-site energies (eV, NIST LDA atomic eigenvalues) and valence electron count."""
    from pymatgen.core import Element
    el = Element(sym)
    orb = el.atomic_orbitals  # Hartree
    Ha = 27.211386
    if sym == "H":
        return dict(orbs=["s"], es=orb["1s"] * Ha, ep=None, nval=1)
    ns = max(int(k[:-1]) for k in orb if k.endswith("s"))
    es = orb[f"{ns}s"] * Ha
    ep = orb.get(f"{ns}p")
    note = None
    if ep is None:
        ep = es * 0.5
        note = f"{sym}: unoccupied p level estimated as eps_s/2"
    else:
        ep *= Ha
    nval = 0
    for (n, l, occ) in el.full_electronic_structure:
        if n == ns and l in ("s", "p"):
            nval += occ
    return dict(orbs=["s", "x", "y", "z"], es=es, ep=ep, nval=nval, note=note)


class UniversalTB:
    def __init__(self, cell, frac, symbols, masses, pdims, shell=1.2):
        from ase import Atoms
        from ase.neighborlist import neighbor_list
        self.cell = np.array(cell, float)
        self.pos = np.asarray(frac) @ self.cell
        self.nat = len(symbols)
        self.masses = np.asarray(masses, float)
        self.pdims = pdims
        self.rec = 2 * np.pi * np.linalg.inv(self.cell).T
        lev = [atomic_levels(s) for s in symbols]
        self.notes = [l["note"] for l in lev if l.get("note")]
        self.orbs, self.off, onsite = [], [], []
        for l in lev:
            self.off.append(len(onsite))
            self.orbs.append(l["orbs"])
            for o in l["orbs"]:
                onsite.append(l["es"] if o == "s" else l["ep"])
        self.norb = len(onsite)
        self.nval = int(sum(l["nval"] for l in lev))
        self.onsite = np.array(onsite)
        at = Atoms(symbols, cell=self.cell, positions=self.pos, pbc=pdims)
        i, j, D, S = neighbor_list("ijDS", at, 6.0)
        dist = np.linalg.norm(D, axis=1)
        dmin = np.array([dist[i == a].min() for a in range(self.nat)])
        keep = dist < shell * np.minimum(dmin[i], dmin[j])
        i, j, D, S = i[keep], j[keep], D[keep], S[keep]
        self.bonds = list(zip(i, j, D, S))
        Rs = sorted(set(map(tuple, S)) | {(0, 0, 0)})
        self.Rint = np.array(Rs)
        self.R = self.Rint @ self.cell
        ridx = {r: k for k, r in enumerate(Rs)}
        nR = len(Rs)
        self.HR = np.zeros((nR, self.norb, self.norb))
        self.HR[ridx[(0, 0, 0)]] += np.diag(onsite)
        self.DP = np.zeros((nR, 3 * self.nat, self.norb, self.norb))
        self.DM = np.zeros((nR, 3 * self.nat, self.norb, self.norb))
        # bond-vector representation for optics (Peierls velocity)
        self.bd, self.bblk = [], []
        h = 1e-5
        for (a, b, d, s) in self.bonds:
            blk = _sk_block(d, self.orbs[a], self.orbs[b], harrison_V)
            r = ridx[tuple(s)]
            sa = slice(self.off[a], self.off[a] + len(self.orbs[a])); sb = slice(self.off[b], self.off[b] + len(self.orbs[b]))
            self.HR[r][sa, sb] += blk
            for al in range(3):
                e = np.zeros(3); e[al] = h
                dblk = (_sk_block(d + e, self.orbs[a], self.orbs[b], harrison_V) -
                        _sk_block(d - e, self.orbs[a], self.orbs[b], harrison_V)) / (2 * h)
                self.DP[r, 3 * b + al][sa, sb] += dblk
                self.DM[r, 3 * a + al][sa, sb] += dblk
            self.bd.append(d); self.bblk.append((sa, sb, blk))
        self.nbonds_per_atom = np.bincount([b[0] for b in self.bonds], minlength=self.nat)
        self.scissor = 0.0

    def hk(self, kc):
        ph = np.exp(1j * kc @ self.R.T)
        return np.einsum("kr,rab->kab", ph, self.HR)

    def eig(self, kc):
        e, v = np.linalg.eigh(self.hk(kc))
        if self.scissor:
            e[:, self.nval // 2:] += self.scissor
        return e, v

    def kgrid(self, nk, shift=0.5):
        axes = [((np.arange(n) + (shift if n > 1 else 0.0)) / n) for n in nk]
        f = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
        return f, f @ self.rec


def tb_optics(tb, nk, pol, hv, sigma=0.1, chunk=4000):
    """Interband epsilon(omega) along unit vector pol from the Peierls velocity of the TB model
    (scissor-corrected energies, velocities rescaled by the scissor ratio); KK for eps1.
    Returns complex refractive index at photon energy hv (eV) and the eps2 spectrum."""
    _, kc = tb.kgrid(nk)
    V = abs(np.linalg.det(tb.cell))
    edges = np.linspace(0.0, 40.0, 4001)
    w = 0.5 * (edges[1:] + edges[:-1])
    e2 = np.zeros_like(w)
    nv = tb.nval // 2
    pol = np.asarray(pol, float); pol /= np.linalg.norm(pol)
    for s in range(0, len(kc), chunk):
        k = kc[s:s + chunk]
        H = np.zeros((len(k), tb.norb, tb.norb), complex)
        dH = np.zeros_like(H)
        for d, (sa, sb, blk) in zip(tb.bd, tb.bblk):
            ph = np.exp(1j * k @ d)
            H[:, sa, sb] += ph[:, None, None] * blk
            dH[:, sa, sb] += (1j * (d @ pol) * ph)[:, None, None] * blk
        H += np.diag(tb.onsite)[None]
        E, U = np.linalg.eigh(H)
        vm = np.einsum("kam,kab,kbn->kmn", U.conj(), dH, U)
        dE0 = E[:, nv:, None] - E[:, None, :nv]
        dE = dE0 + tb.scissor
        M2 = np.abs(vm[:, nv:, :nv]) ** 2 * (dE / np.maximum(dE0, 1e-6)) ** 2
        x = dE.reshape(-1); y = (M2 / dE ** 2).reshape(-1)
        e2 += np.histogram(x, bins=edges, weights=y)[0]
    dw = w[1] - w[0]
    kern = np.exp(-0.5 * (np.arange(-int(5 * sigma / dw), int(5 * sigma / dw) + 1) * dw / sigma) ** 2)
    e2 = np.convolve(e2 / dw, kern / kern.sum(), mode="same")
    e2 *= np.pi * E2_EPS0 / V * 2 / len(kc)
    # Kramers-Kronig
    e1 = np.ones_like(w)
    for i, wi in enumerate(w):
        den = w ** 2 - wi ** 2
        den[i] = np.inf
        e1[i] += 2 / np.pi * np.sum(w * e2 / den) * dw
    eps = np.interp(hv, w, e1) + 1j * np.interp(hv, w, e2)
    return np.sqrt(eps), w, e1, e2


def film_absorption(n, d_m, lam_m):
    """Normal-incidence absorptance 1-R-T of a free-standing film (repo graphite formula)."""
    n0 = 1.0
    r01 = (n0 - n) / (n0 + n); r12 = (n - n0) / (n + n0)
    t01 = 2 * n0 / (n0 + n); t12 = 2 * n / (n + n0)
    b = 2 * np.pi * n * d_m / lam_m
    r = (r01 + r12 * np.exp(2j * b)) / (1 + r01 * r12 * np.exp(2j * b))
    t = t01 * t12 * np.exp(1j * b) / (1 + r01 * r12 * np.exp(2j * b))
    return float(1 - abs(r) ** 2 - abs(t) ** 2), float(abs(r) ** 2), float(abs(t) ** 2)


def electronic_setup(tb, args, pdims, eg_input):
    """Band edges, scissor, metal/semiconductor classification, DOS histogram."""
    nk = args.nk_dos
    _, kc = tb.kgrid(nk)
    E = np.concatenate([np.linalg.eigvalsh(tb.hk(kc[s:s + 20000])) for s in range(0, len(kc), 20000)])
    nv = tb.nval // 2
    # band edges also at the high-symmetry points (band extrema often sit there, off the shifted grid)
    sp = np.array(list(special_points(tb.cell, pdims).values())) @ tb.rec
    Esp = np.linalg.eigvalsh(tb.hk(sp))
    if tb.nval % 2:
        log("[eph] odd number of valence electrons per cell: treated as a metal")
    evbm = max(E[:, nv - 1].max(), Esp[:, nv - 1].max()) if nv > 0 else E.min()
    ecbm = min(E[:, nv].min(), Esp[:, nv].min()) if nv < tb.norb else E.max()
    eg_tb = float(ecbm - evbm)
    if eg_input is None:
        metal = (eg_tb < 0.05) or (tb.nval % 2 == 1)
        eg = 0.0 if metal else eg_tb
    else:
        metal = eg_input <= 0.0 or tb.nval % 2 == 1
        eg = max(eg_input, 0.0)
    tb.scissor = 0.0 if metal else eg - eg_tb
    if not metal:
        E[:, nv:] += tb.scissor
        ecbm += tb.scissor
    return dict(metal=bool(metal), eg_tb=eg_tb, eg=float(eg), evbm=float(evbm), ecbm=float(ecbm), E=E, nk=nk)


def dos_hist(E, lo, hi, de=0.005):
    edges = np.arange(lo, hi + de, de)
    h, _ = np.histogram(E.reshape(-1), bins=edges)
    return 0.5 * (edges[1:] + edges[:-1]), 2.0 * h / E.shape[0]   # states per cell (both spins) per bin


def fermi(x):
    return 0.5 * (1 - np.tanh(0.5 * x))


def solve_mu(ec, g, N, T, lo, hi, holes=False):
    for _ in range(100):
        mu = 0.5 * (lo + hi)
        f = fermi((ec - mu) / (KB_eV * T))
        n = np.sum(g * (1 - f)) if holes else np.sum(g * f)
        if (n > N) ^ holes:
            hi = mu
        else:
            lo = mu
    return 0.5 * (lo + hi)


def electron_tables(el, n_exc, T0):
    """mu(Te) (metal) or quasi-Fermi levels (semiconductor), and electronic energy E(Te) per cell (eV)."""
    E = el["E"]
    Tt = np.geomspace(T0, 60000, 400)
    if el["metal"]:
        ec, g = dos_hist(E, E.min() - 1, E.max() + 1)
        Nel = el["nval"]
        mus, Es = [], []
        for T in Tt:
            mu = solve_mu(ec, g, Nel, T, ec[0], ec[-1])
            mus.append(mu)
            Es.append(np.sum(g * ec * fermi((ec - mu) / (KB_eV * T))))
        return dict(T=Tt, mu=np.array(mus), E=np.array(Es) - Es[0])
    nv = el["nval"] // 2
    ecc, gc = dos_hist(E[:, nv:], el["ecbm"] - 0.01, el["ecbm"] + 6)
    ecv, gv = dos_hist(E[:, :nv], el["evbm"] - 6, el["evbm"] + 0.01)
    muc, muv, Es = [], [], []
    for T in Tt:
        mc = solve_mu(ecc, gc, n_exc, T, el["ecbm"] - 3, el["ecbm"] + 6)
        mv = solve_mu(ecv, gv, n_exc, T, el["evbm"] - 6, el["evbm"] + 3, holes=True)
        muc.append(mc); muv.append(mv)
        Es.append(np.sum(gc * (ecc - el["ecbm"]) * fermi((ecc - mc) / (KB_eV * T))) +
                  np.sum(gv * (el["evbm"] - ecv) * (1 - fermi((ecv - mv) / (KB_eV * T)))))
    return dict(T=Tt, muc=np.array(muc), muv=np.array(muv), E=np.array(Es) - Es[0])


def eph_W_point(tb, kc, Ek, Ck, bsel, q_cart, freqs, evecs, classes, ebins, sig_e, nk_total, fmin, chunk=2000):
    """Energy-resolved coupling W[class, nu, bin] (eV) for one phonon wavevector."""
    nb = len(freqs)
    ncls = len(classes)
    W = np.zeros((ncls, nb, len(ebins) - 1))
    tau = tb.pos
    phase_tau = np.exp(1j * tau @ q_cart)          # (nat,)
    ok = freqs > fmin
    if not ok.any():
        return W
    w_rad = 2 * np.pi * np.maximum(freqs, 1e-6) * 1e12
    l = np.sqrt(HBAR_Js / (2 * tb.masses[:, None] * AMU * w_rad[None, :])) * 1e10   # (nat, nb) A
    U = evecs.reshape(tb.nat, 3, nb) * l[:, None, :] * phase_tau[:, None, None]
    U = U.reshape(3 * tb.nat, nb)
    hw = freqs * THz2eV
    nR = len(tb.R)
    DPf = tb.DP.reshape(nR, -1); DMf = tb.DM.reshape(nR, -1)
    nsel = len(bsel)
    cls_mask = [np.outer(np.isin(bsel, c), np.isin(bsel, c)) for c in classes]   # (m, n)
    for s in range(0, len(kc), chunk):
        k = kc[s:s + chunk]; nk = len(k)
        kq = k + q_cart
        Hkq = np.einsum("kr,rab->kab", np.exp(1j * kq @ tb.R.T), tb.HR)
        Ekq, Ckq = np.linalg.eigh(Hkq)
        if tb.scissor:
            Ekq[:, tb.nval // 2:] += tb.scissor
        Ekq = Ekq[:, bsel]; Ckq = Ckq[:, :, bsel]
        M = (np.exp(1j * kq @ tb.R.T) @ DPf - np.exp(1j * k @ tb.R.T) @ DMf).reshape(nk, 3 * tb.nat, tb.norb, tb.norb)
        T_ = M @ Ck[s:s + chunk][:, None, :, :]                  # (k, x, norb, nsel)
        gx = np.conj(Ckq).transpose(0, 2, 1)[:, None] @ T_       # (k, x, nsel, nsel)
        g = np.einsum("kxmn,xv->kvmn", gx, U)                     # (k, nu, m, n)
        g2 = np.abs(g) ** 2
        de = Ekq[:, None, :, None] - Ek[s:s + chunk][:, None, None, :] - hw[None, :, None, None]
        wgt = g2 * np.exp(-0.5 * (de / sig_e) ** 2) / (np.sqrt(2 * np.pi) * sig_e)
        ib = np.digitize(Ek[s:s + chunk], ebins) - 1                # (k, n)
        inb = (ib >= 0) & (ib < len(ebins) - 1)
        for c in range(ncls):
            wc = (wgt * cls_mask[c][None, None]).sum(2)            # (k, nu, n)
            for nu in range(nb):
                if not ok[nu]:
                    continue
                W[c, nu] += np.bincount(ib[inb], weights=wc[:, nu, :][inb], minlength=len(ebins) - 1)[:len(ebins) - 1]
    return W * 2.0 / nk_total


def A_from_W(W, hw, ebins, el, tabs, Te_grid):
    """Generation coefficients A[q, nu, T] (1/ps) from W[q, class, nu, bin]."""
    ec = 0.5 * (ebins[1:] + ebins[:-1])
    Nq, ncls, nb, _ = W.shape
    A = np.zeros((Nq, nb, len(Te_grid)))
    for it, T in enumerate(Te_grid):
        kT = KB_eV * T
        if el["metal"]:
            mus = [np.interp(T, tabs["T"], tabs["mu"])]
        else:
            mus = [np.interp(T, tabs["T"], tabs["muv"]), np.interp(T, tabs["T"], tabs["muc"])]
        for c in range(ncls):
            f1 = fermi((ec[None, None, :] - mus[c]) / kT)
            f2 = fermi((ec[None, None, :] + hw[:, :, None] - mus[c]) / kT)
            A[:, :, it] += np.sum(W[:, c] * (f1 * (1 - f2)), axis=-1)
    return A * (2 * np.pi / HBAR_eVs) * 1e-12


def stage_eph(outdir, ph3, tb, el, mesh, freq, qfrac, args, pdims, n_exc, T0):
    path = os.path.join(outdir, "eph_W.npz")
    key = _key(dict(fc=str(np.load(os.path.join(outdir, "fc_key.npz"))["key"]), mesh=mesh, nk=args.nk_eph,
                    sig=args.sig_e, nsub=args.nsub, win=args.e_window, eg=el["eg"], metal=el["metal"],
                    refine=args.refine_cells, v=3))
    nv = tb.nval // 2
    if el["metal"]:
        mu0 = float(np.interp(T0, el["tabs"]["T"], el["tabs"]["mu"]))
        lo, hi = mu0 - args.e_window, mu0 + args.e_window
    else:
        lo, hi = el["evbm"] - args.e_window, el["ecbm"] + args.e_window
    ebins = np.arange(lo, hi + 0.02, 0.02)
    if _cached(path, key):
        log("[eph] cached coupling function W")
        z = np.load(path)
        return z["W"], z["ebins"], json.loads(str(z["info"]))
    ph = make_phonopy(ph3)
    kf, kc = tb.kgrid(args.nk_eph)
    Ek, Ck = [], []
    for s in range(0, len(kc), 20000):
        e, v = tb.eig(kc[s:s + 20000]); Ek.append(e); Ck.append(v)
    Ek = np.concatenate(Ek); Ck = np.concatenate(Ck)
    inwin = (Ek >= lo) & (Ek <= hi)
    bsel = np.where(inwin.any(0))[0]
    ksel = np.where(inwin.any(1))[0]
    if el["metal"]:
        classes = [bsel]
    else:
        classes = [bsel[bsel < nv], bsel[bsel >= nv]]
    log(f"[eph] TB: {tb.norb} orbitals, {len(tb.R)} lattice vectors, bands in window {bsel.tolist()}, "
        f"k-points {len(ksel)}/{len(kc)} in window [{lo:.2f}, {hi:.2f}] eV")
    kcs, Eks, Cks = kc[ksel], Ek[ksel][:, bsel], Ck[ksel]
    Cks = Cks[:, :, bsel]
    Nq = len(qfrac); nb = freq.shape[1]
    dims = [i for i in range(3) if pdims[i]]

    def offsets(ns):
        o = (np.arange(ns) + 0.5) / ns - 0.5
        grids = np.meshgrid(*[o] * len(dims), indexing="ij")
        out = np.zeros((ns ** len(dims), 3))
        for a, i in enumerate(dims):
            out[:, i] = grids[a].reshape(-1)
        return out / np.array(mesh)[None, :]

    def cell_W(iq, ns):
        subs = qfrac[iq][None, :] + offsets(ns)
        fs, evs = phonons_at(ph, subs)
        acc = np.zeros((len(classes), nb, len(ebins) - 1))
        for isub, q in enumerate(subs):
            acc += eph_W_point(tb, kcs, Eks, Cks, bsel, q @ tb.rec, fs[isub], evs[isub], classes, ebins,
                               args.sig_e, len(kc), args.fmin_thz)
        return acc / len(subs)

    t0 = time.time()
    W = np.zeros((Nq, len(classes), nb, len(ebins) - 1))
    for iq in range(Nq):
        W[iq] = cell_W(iq, args.nsub)
        if iq % max(1, Nq // 10) == 0:
            log(f"[eph] cell {iq}/{Nq} ({time.time()-t0:.0f}s)")
    # refinement of the cells that carry most of the emitted power at a hot electron temperature
    hw = freq * THz2eV
    Tref = 3000.0
    Aref = A_from_W(W, hw, ebins, el, el["tabs"], [Tref])[:, :, 0]
    P = (Aref * hw * np.where(freq > args.fmin_thz, 1, 0)).sum(1)
    order = np.argsort(P)[::-1]
    cum = np.cumsum(P[order]) / max(P.sum(), 1e-30)
    nref = int(min(args.refine_cells, np.searchsorted(cum, 0.5) + 1))
    for c, iq in enumerate(order[:nref]):
        old = P[iq]
        W[iq] = cell_W(iq, max(2 * args.nsub, 3))
        log(f"[eph] refined cell {c+1}/{nref} q={np.round(qfrac[iq], 3)} power share {P[iq]/P.sum():.3f} "
            f"({time.time()-t0:.0f}s)")
    info = dict(n_k_window=int(len(ksel)), n_k_total=int(len(kc)), bands=bsel.tolist(), n_refined=nref,
                seconds=time.time() - t0)
    np.savez(path, key=key, W=W, ebins=ebins, info=json.dumps(info))
    return W, ebins, info


# ----------------------------------------------------------------------------------------------
# 6. kinetics
def _make_phph_kernel():
    """Parallel three-phonon collision sum (numba); falls back to numpy when numba is unavailable."""
    try:
        import numba
    except Exception:
        return None

    @numba.njit(parallel=True, cache=True)
    def kern(nf, a, b, c, v, kind, size, nthr):
        npro = a.shape[0]
        buf = np.zeros((nthr, size))
        step = (npro + nthr - 1) // nthr
        for t in numba.prange(nthr):
            lo = t * step
            hi = min(npro, lo + step)
            for i in range(lo, hi):
                na = nf[a[i]]; n1 = nf[b[i]]; n2 = nf[c[i]]
                if kind == 0:
                    r = (na + 1) * n1 * n2 - na * (n1 + 1) * (n2 + 1)
                elif kind == 1:
                    r = (na + 1) * (n1 + 1) * n2 - na * n1 * (n2 + 1)
                else:
                    r = (na + 1) * n1 * (n2 + 1) - na * (n1 + 1) * n2
                buf[t, a[i]] += v[i] * r
        return buf.sum(0)
    try:
        numba.set_num_threads(NTHREADS)
    except Exception:
        pass
    return kern


_PHPH_KERNEL = _make_phph_kernel() if os.environ.get("RM_NUMBA", "1") == "1" else None


class Kinetics:
    def __init__(self, freq, procs, Atab, Tgrid, etab, args, E_abs_cell, frac_e, metal, fmin):
        self.freq = freq
        self.Nq, self.nb = freq.shape
        self.hw = freq * THz2eV
        self.valid = freq > fmin
        self.Atab = Atab; self.logT = np.log(Tgrid); self.Tgrid = Tgrid
        self.etab = etab
        self.procs = procs
        self.E_abs = E_abs_cell; self.frac_e = frac_e; self.metal = metal
        self.fwhm = args.fwhm_fs * 1e-3
        self.T0 = args.T
        self.trec = args.recomb_ps
        self.cw0 = None

    def Te_of(self, Ee):
        return float(np.interp(Ee, self.etab["E"], self.etab["T"]))

    def power(self, t):
        s = self.fwhm / (2 * np.sqrt(2 * np.log(2)))
        return self.E_abs * np.exp(-0.5 * (t / s) ** 2) / (np.sqrt(2 * np.pi) * s)

    def injected(self, t):
        s = self.fwhm / (2 * np.sqrt(2 * np.log(2)))
        return 0.5 * (1 + math.erf(t / (np.sqrt(2) * s)))

    def phph_rate(self, n):
        nf = np.ascontiguousarray(n.reshape(-1), dtype=np.float64)
        if _PHPH_KERNEL is not None:
            out = np.zeros(self.Nq * self.nb)
            for kind, (a, b, c, v) in enumerate(self.procs):
                out += _PHPH_KERNEL(nf, a, b, c, v, kind, self.Nq * self.nb, NTHREADS)
            return out.reshape(self.Nq, self.nb)
        out = np.zeros(self.Nq * self.nb)
        for kind, (a, b, c, v) in enumerate(self.procs):
            na, n1, n2 = nf[a], nf[b], nf[c]
            if kind == 0:
                r = (na + 1) * n1 * n2 - na * (n1 + 1) * (n2 + 1)
            elif kind == 1:
                r = (na + 1) * (n1 + 1) * n2 - na * n1 * (n2 + 1)
            else:
                r = (na + 1) * n1 * (n2 + 1) - na * (n1 + 1) * n2
            out += np.bincount(a, weights=v * r, minlength=self.Nq * self.nb)
        return out.reshape(self.Nq, self.nb)

    def eph_rate(self, n, Te, scale):
        x = np.clip((np.log(Te) - self.logT[0]) / (self.logT[1] - self.logT[0]), 0, len(self.Tgrid) - 1.000001)
        i = int(x); f = x - i
        A = ((1 - f) * self.Atab[:, :, i] + f * self.Atab[:, :, i + 1]) * scale
        boltz = np.exp(-self.hw / (KB_eV * Te))
        return np.where(self.valid, A * ((n + 1) * boltz - n), 0.0)

    def carrier_scale(self, t, y_rec):
        if self.metal:
            return 1.0
        return self.injected(t) * y_rec

    def rhs(self, t, y):
        Ee = y[0]; rec = y[1]; n = y[2:].reshape(self.Nq, self.nb)
        Te = max(self.Te_of(Ee), self.T0 * 0.5)
        r_ep = self.eph_rate(n, Te, self.carrier_scale(t, rec))
        r_pp = self.phph_rate(n)
        resid = np.sum(self.hw * r_pp)
        cw = np.where(self.valid, self.hw * n * (n + 1), 0.0)
        r_pp = r_pp - resid * cw / np.sum(self.hw * cw)
        dEe = self.power(t) * self.frac_e - np.sum(self.hw * r_ep) / self.Nq
        drec = 0.0
        if (not self.metal) and self.trec:
            # recombination releases the gap energy of the recombining pairs into the carrier bath
            drec = -rec / self.trec
            dEe += -drec * self.injected(t) * self.E_abs * (1 - self.frac_e)
        return np.concatenate([[dEe, drec], (r_ep + r_pp).reshape(-1)])

    def jac(self, t, y):
        """Jacobian of the electron-phonon part (the stiff part); the three-phonon term is left to
        the Newton iterations of the implicit integrator."""
        from scipy.sparse import coo_matrix
        Ee = y[0]; rec = y[1]; n = y[2:].reshape(self.Nq, self.nb)
        Te = max(self.Te_of(Ee), self.T0 * 0.5)
        s = self.carrier_scale(t, rec)
        if not hasattr(self, "_dTdE"):
            self._dTdE = np.gradient(self.etab["T"]) / np.maximum(np.gradient(self.etab["E"]), 1e-30)
        dTdE = float(np.interp(Ee, self.etab["E"], self._dTdE))
        x = np.clip((np.log(Te) - self.logT[0]) / (self.logT[1] - self.logT[0]), 0, len(self.Tgrid) - 1.000001)
        i = int(x); f = x - i
        A = ((1 - f) * self.Atab[:, :, i] + f * self.Atab[:, :, i + 1]) * s
        b = np.exp(-self.hw / (KB_eV * Te))
        dn_diag = np.where(self.valid, -A * (1 - b), 0.0).reshape(-1)
        dT = 1e-3 * Te
        dr_dT = ((self.eph_rate(n, Te + dT, s) - self.eph_rate(n, Te - dT, s)) / (2 * dT)).reshape(-1)
        M = self.Nq * self.nb
        hw = self.hw.reshape(-1)
        idx = np.arange(M) + 2
        rows = np.concatenate([idx, idx, np.zeros(M, int), [0]])
        cols = np.concatenate([idx, np.zeros(M, int), idx, [0]])
        vals = np.concatenate([dn_diag, dr_dT * dTdE, -hw * dn_diag / self.Nq,
                               [-np.sum(hw * dr_dT) * dTdE / self.Nq]])
        return coo_matrix((vals, (rows, cols)), shape=(M + 2, M + 2)).tocsc()


def eff_temperature(hw_branch, E_target, Tlo=1.0, Thi=1e5):
    """Temperature at which sum_q hw nB(hw,T) equals E_target (bisection, vectorized over time)."""
    lo = np.full(E_target.shape, Tlo); hi = np.full(E_target.shape, Thi)
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        Em = np.sum(hw_branch[None] / np.expm1(np.minimum(hw_branch[None] / (KB_eV * mid[:, None]), 700)), axis=1)
        big = Em > E_target
        hi = np.where(big, mid, hi); lo = np.where(big, lo, mid)
    return np.sqrt(lo * hi)


# ----------------------------------------------------------------------------------------------
# 7. observables and fits
def gauss_conv(t_uniform, y, fwhm):
    if fwhm <= 0:
        return y
    dt = t_uniform[1] - t_uniform[0]
    s = fwhm / (2 * np.sqrt(2 * np.log(2)))
    m = int(np.ceil(4 * s / dt))
    k = np.exp(-0.5 * (np.arange(-m, m + 1) * dt / s) ** 2); k /= k.sum()
    yp = np.concatenate([np.full(m, y[0]), y, np.full(m, y[-1])])
    return np.convolve(yp, k, mode="valid")


MODELS = {
    "biexp": (lambda x, A1, t1, A2, t2: np.where(x > 0, A1 * (1 - np.exp(-x / t1)) + A2 * (1 - np.exp(-x / t2)), 0.0),
              "A1*(1-exp(-t/tau1)) + A2*(1-exp(-t/tau2)) for t>0, convolved with a Gaussian of FWHM = irf"),
    "single": (lambda x, A, t1: np.where(x > 0, A * (1 - np.exp(-x / t1)), 0.0),
               "A*(1-exp(-t/tau)) for t>0, convolved with a Gaussian of FWHM = irf"),
    "rise_decay": (lambda x, A, tr, td, C: np.where(x > 0, A * (1 - np.exp(-x / tr)) * np.exp(-x / td) + C * (1 - np.exp(-x / td)), 0.0),
                   "A*(1-exp(-t/tau_r))*exp(-t/tau_d) + C*(1-exp(-t/tau_d)) for t>0, convolved with a Gaussian of FWHM = irf"),
}


def fit_trace(t, y, model, irf_ps, probe_ps, window=None):
    """Fit y(t) (simulated, pump-limited) with model (x) Gaussian IRF, after convolving the simulated
    curve with the probe duration. Returns parameters dict."""
    from scipy.optimize import curve_fit
    tu = np.arange(t[0], t[-1], 0.005 if t[-1] < 30 else 0.01)
    yu = np.interp(tu, t, y)
    yu = gauss_conv(tu, yu, probe_ps)
    if window:
        m = (tu >= window[0]) & (tu <= window[1])
        tu, yu = tu[m], yu[m]
    f, _ = MODELS[model]
    fun = lambda x, *p: gauss_conv(x, f(x, *p), irf_ps) if irf_ps > 0 else f(x, *p)
    span = tu[-1]
    ymax = yu[np.argmax(np.abs(yu))]
    best = None
    if model == "biexp":
        starts = [[0.5 * yu[-1], 0.3, 0.5 * yu[-1], 5.0], [0.7 * ymax, 1.0, yu[-1] - 0.7 * ymax, 20.0],
                  [ymax, 0.2, yu[-1] - ymax, 3.0]]
        bounds = ([-np.inf, 0.005, -np.inf, 0.005], [np.inf, span * 5, np.inf, span * 20])
    elif model == "single":
        starts = [[yu[-1], 1.0], [ymax, 0.3], [ymax, 10.0]]
        bounds = ([-np.inf, 0.005], [np.inf, span * 20])
    else:
        starts = [[ymax, 0.2, 2.0, yu[-1]], [ymax, 0.5, 10.0, yu[-1]], [ymax, 0.1, 1.0, yu[-1]]]
        bounds = ([-np.inf, 0.005, 0.01, -np.inf], [np.inf, span, span * 20, np.inf])
    for p0 in starts:
        try:
            p, _ = curve_fit(fun, tu, yu, p0=p0, bounds=bounds, maxfev=20000)
            r = np.sqrt(np.mean((fun(tu, *p) - yu) ** 2)) / max(np.abs(yu).max(), 1e-30)
            if best is None or r < best[1]:
                best = (p, r)
        except Exception:
            continue
    if best is None:
        return None
    p, r = best
    names = {"biexp": ["A1", "tau1_ps", "A2", "tau2_ps"], "single": ["A", "tau_ps"],
             "rise_decay": ["A", "tau_rise_ps", "tau_decay_ps", "C"]}[model]
    out = {k: float(v) for k, v in zip(names, p)}
    if model == "biexp" and out["tau1_ps"] > out["tau2_ps"]:
        out = {"A1": out["A2"], "tau1_ps": out["tau2_ps"], "A2": out["A1"], "tau2_ps": out["tau1_ps"]}
    out["rel_rms_residual"] = float(r)
    return out


def rise_time(t, y, frac=1 - np.exp(-1)):
    y = y - y[0]
    pk = np.argmax(y)
    if y[pk] <= 0:
        return None, None
    k = np.where(y >= frac * y[pk])[0][0]
    return float(t[k]), float(t[pk])


def decay_time(t, y):
    """1/e time of the excess above the final value, measured from the peak."""
    y = y - y[0]
    pk = np.argmax(y)
    yf = y[-1]
    rel = (y[pk:] - yf) / max(y[pk] - yf, 1e-30)
    k = np.where(rel <= np.exp(-1))[0]
    return float(t[pk + k[0]] - t[pk]) if len(k) else None


def scattering_factor(sym, Q):
    """Screened-Coulomb (Thomas-Fermi/Wentzel) first Born electron scattering amplitude, arbitrary units."""
    from pymatgen.core import Element
    Z = Element(sym).Z
    aTF = 0.8853 * 0.529177 * Z ** (-1 / 3)
    return Z / (np.asarray(Q) ** 2 + aTF ** -2)


def special_points(prim_cell, pdims):
    from ase.cell import Cell
    c = Cell(prim_cell)
    try:
        bp = c.bandpath(npoints=0, pbc=pdims)
        sp = bp.special_points
    except Exception:
        sp = {"G": np.zeros(3)}
    return {k: np.asarray(v, float) for k, v in sp.items()}


def star_cart(q_cart, rots_cart):
    pts = []
    for R in rots_cart:
        for s in (1, -1):
            v = s * (R @ q_cart)
            if not any(np.allclose(v, p, atol=1e-6) for p in pts):
                pts.append(v)
    return pts


def nearest_grid(qc, rec, mesh, pdims):
    f = np.linalg.solve(rec.T, qc)
    m = np.array(mesh)
    idx = np.round(f * m).astype(int)
    f_grid = idx / m
    shift = f_grid - f
    i = idx % m
    gi = int(i[0] + i[1] * m[0] + i[2] * m[0] * m[1])
    return gi, f_grid @ rec, float(np.linalg.norm(shift @ rec))


# ----------------------------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("structure")
    p.add_argument("--pump-eV", type=float, required=True)
    p.add_argument("--fluence-mJcm2", type=float, required=True, help="incident fluence")
    p.add_argument("--fwhm-fs", type=float, required=True)
    p.add_argument("--T", type=float, default=300.0, help="base temperature (K)")
    p.add_argument("--thickness-nm", type=float, default=None, help="film thickness (free-standing); default semi-infinite")
    p.add_argument("--tmax-ps", type=float, default=50.0)
    p.add_argument("--absorbed-fluence-mJcm2", type=float, default=None, help="override optics")
    p.add_argument("--nk", type=str, default=None, help="complex refractive index at the pump, e.g. 3.19+0.29j")
    p.add_argument("--pol", type=str, default="1,0,0", help="pump polarization (Cartesian) for TB optics")
    p.add_argument("--gap-eV", type=float, default=None, help="band gap for the scissor (default: built-in value or TB)")
    p.add_argument("--layer-thickness", type=float, default=None, help="A per layer, for 2D cells with vacuum")
    p.add_argument("--recomb-ps", type=float, default=None, help="carrier recombination time (default: none)")
    p.add_argument("--name", type=str, default=None)
    p.add_argument("--outdir", type=str, default=None)
    p.add_argument("--qmesh", type=str, default=None, help="e.g. 24,24,1")
    p.add_argument("--sigma-thz", type=float, default=None, help="energy-conservation width of the 3-phonon kernel")
    p.add_argument("--fc2-length", type=float, default=18.0)
    p.add_argument("--fc3-length", type=float, default=10.0)
    p.add_argument("--fmin-thz", type=float, default=0.05)
    p.add_argument("--nk-eph", type=str, default=None)
    p.add_argument("--nk-dos", type=str, default=None)
    p.add_argument("--sig-e", type=float, default=None, help="energy-conservation width of e-ph deltas (eV)")
    p.add_argument("--nsub", type=int, default=None, help="q sub-points per periodic direction per cell")
    p.add_argument("--refine-cells", type=int, default=None, help="default 16 (2D), 8 (3D)")
    p.add_argument("--e-window", type=float, default=None)
    p.add_argument("--no-relax", action="store_true")
    p.add_argument("--beam", type=str, default=None, help="electron beam direction [uvw] of the input cell")
    p.add_argument("--qmax", type=float, default=7.0, help="max |Q| (1/A) for diffuse scattering")
    p.add_argument("--bragg", type=str, default=None, help="Bragg peaks hkl of the input cell, e.g. '4 0 0;0 0 4'")
    p.add_argument("--diffuse-zone", type=str, default=None, help="Bragg peak hkl whose zone is used for diffuse traces")
    p.add_argument("--irf-fs", type=float, default=None, help="instrument response FWHM for fits (default pump FWHM)")
    p.add_argument("--eph-scale", type=float, default=1.0, help="sensitivity analysis only: scales all A_qnu")
    p.add_argument("--tag", type=str, default="", help="suffix of output files (sensitivity runs)")
    p.add_argument("--reuse-dynamics", action="store_true", help="re-use dynamics<tag>.npz if its settings match")
    p.add_argument("--stop-after", type=str, default=None, choices=[None, "fc", "phph", "eph"])
    return p.parse_args(argv)


def validate_args(a):
    errs = []
    if not (0.05 < a.pump_eV < 20): errs.append("--pump-eV must be in (0.05, 20) eV")
    if not (0 < a.fluence_mJcm2 < 1000): errs.append("--fluence-mJcm2 must be in (0, 1000)")
    if not (0 < a.fwhm_fs < 5000): errs.append("--fwhm-fs must be in (0, 5000)")
    if not (1 <= a.T <= 3000): errs.append("--T must be in [1, 3000] K")
    if a.thickness_nm is not None and not (0 < a.thickness_nm < 1e6): errs.append("--thickness-nm must be positive")
    if not (0.5 <= a.tmax_ps <= 2000): errs.append("--tmax-ps must be in [0.5, 2000]")
    if a.gap_eV is not None and a.gap_eV < 0: errs.append("--gap-eV must be >= 0")
    if errs:
        raise SystemExit("invalid arguments: " + "; ".join(errs))


def parse_triplet(s, typ=float):
    return [typ(x) for x in s.replace(",", " ").split()]


def main(argv=None):
    args = parse_args(argv)
    validate_args(args)
    T_start = time.time()
    atoms, name, meta = load_structure(args.structure)
    name = args.name or name
    outdir = args.outdir or os.path.join(ROOT, "materials", "blind", name)
    os.makedirs(outdir, exist_ok=True)
    warns = validate_structure(atoms)
    for w in warns:
        log("WARNING:", w)
    pdims = periodic_dims(atoms)
    dim = int(pdims.sum())
    log(f"[structure] {name}: {atoms.get_chemical_formula()}, {len(atoms)} atoms, periodic dims {pdims.tolist()}")
    timings = {}
    calc_holder = {}

    def calc_getter():
        if "c" not in calc_holder:
            calc_holder["c"] = mace_calc()
        return calc_holder["c"]

    relax_info = None
    relaxed_path = os.path.join(outdir, "structure_relaxed.json")
    if not args.no_relax:
        from ase.io import read, write
        rkey = _key(dict(cell=np.round(np.array(atoms.cell), 5).tolist(), pos=np.round(atoms.positions, 5).tolist(),
                         sym=atoms.get_chemical_symbols()))
        if os.path.exists(relaxed_path) and _cached(os.path.join(outdir, "relax_key.npz"), rkey):
            atoms = read(relaxed_path)
            relax_info = "cached"
        else:
            t0 = time.time()
            atoms, f0, f1 = relax_positions(atoms, calc_getter())
            write(relaxed_path, atoms)
            np.savez(os.path.join(outdir, "relax_key.npz"), key=rkey)
            relax_info = dict(max_force_before=float(f0), max_force_after=float(f1), seconds=time.time() - t0)
            log(f"[relax] max force {f0:.2e} -> {f1:.2e} eV/A")
    t0 = time.time()
    fcinfo = stage_force_constants(atoms, outdir, calc_getter, args, pdims)
    timings["force_constants"] = time.time() - t0
    if args.stop_after == "fc":
        return
    ph3 = load_ph3(outdir)
    prim = ph3.primitive
    prim_cell = np.array(prim.cell)
    nat = len(prim.masses); nb = 3 * nat
    mesh = parse_triplet(args.qmesh, int) if args.qmesh else default_mesh(prim_cell, pdims, nb)
    ph = make_phonopy(ph3)
    fG, _ = phonons_at(ph, [[0, 0, 0]])
    cg = [np.arange(4) / 4 if p else np.zeros(1) for p in pdims]
    qc_ = np.stack(np.meshgrid(*cg, indexing="ij"), -1).reshape(-1, 3)
    fmax = float(phonons_at(ph, qc_)[0].max())
    Neff = np.prod([m for m in mesh if m > 1]) ** (1.0 / dim)
    sig_dyn = args.sigma_thz or 0.405 * fmax / Neff
    log(f"[phonons] primitive cell {nat} atoms, {nb} branches, mesh {mesh}, max freq ~{fmax:.1f} THz, "
        f"3-phonon width {sig_dyn:.3f} THz; Gamma: {np.round(fG[0], 2).tolist()} THz")
    t0 = time.time()
    P = stage_phph(outdir, mesh, sig_dyn, args)
    timings["three_phonon"] = time.time() - t0
    if args.stop_after == "phph":
        return
    freq, qfrac = P["freq"], P["qfrac"]
    Nq = len(qfrac)
    n_invalid = int((freq <= args.fmin_thz).sum())
    # ---------------- electrons
    tb = UniversalTB(prim_cell, prim.scaled_positions, prim.symbols, prim.masses, pdims)
    for nt in tb.notes:
        warns.append(nt); log("WARNING:", nt)
    if args.nk_dos is None:
        args.nk_dos = default_kmesh(tb, pdims, 2.0e5 if dim == 3 else 1.2e5)
    else:
        args.nk_dos = parse_triplet(args.nk_dos, int)
    if args.nk_eph is None:
        args.nk_eph = default_kmesh(tb, pdims, 1.5e4 if dim == 3 else 9.0e4)
    else:
        args.nk_eph = parse_triplet(args.nk_eph, int)
    if args.nsub is None:
        args.nsub = 1 if dim == 3 else 3
    if args.refine_cells is None:
        args.refine_cells = 8 if dim == 3 else 16
    if args.sig_e is None:
        args.sig_e = 0.10 if dim == 3 else 0.05
    eg_in = args.gap_eV if args.gap_eV is not None else meta.get("gap")
    el = electronic_setup(tb, args, pdims, eg_in)
    el["nval"] = tb.nval
    if args.e_window is None:
        args.e_window = 2.0 if el["metal"] else float(np.clip(args.pump_eV - el["eg"] + 0.5, 1.0, 2.5))
    log(f"[eph] universal TB: {tb.nval} valence electrons, {tb.norb} orbitals, bonds/atom {tb.nbonds_per_atom.tolist()}; "
        f"TB gap {el['eg_tb']:.3f} eV -> {'metal/semimetal' if el['metal'] else 'semiconductor, gap %.3f eV (scissor %+.3f)' % (el['eg'], tb.scissor)}")
    # ---------------- optics and absorbed energy
    V_prim = abs(np.linalg.det(prim_cell))
    layer_t = None
    if dim == 2:
        layer_t = args.layer_thickness or meta.get("layer_thickness")
        if layer_t is None:
            ax = int(np.where(~pdims)[0][0])
            z = prim.scaled_positions[:, ax] * np.linalg.norm(prim_cell[ax])
            layer_t = float(z.max() - z.min() + 3.35)
        V_cell = V_prim / np.linalg.norm(prim_cell[int(np.where(~pdims)[0][0])]) * layer_t
    elif dim == 3:
        V_cell = V_prim
    else:
        raise SystemExit("only 2D (with vacuum) and 3D periodic structures are supported")
    lam = HC_eVnm / args.pump_eV * 1e-9
    nk_src = None
    try:
        pol = np.array(parse_triplet(args.pol))
        nk_tb, w_opt, e1, e2 = tb_optics(tb, default_kmesh(tb, pdims, 4e4 if dim == 3 else 6e4), pol, args.pump_eV)
        if dim == 2:
            # 2D sheet: rescale eps-1 from the vacuum cell to the layer volume
            ratio = V_prim / V_cell
            eps = (nk_tb ** 2 - 1) * ratio + 1
            nk_tb = np.sqrt(eps)
    except Exception as exc:
        nk_tb = None
        log("WARNING: TB optics failed:", exc)
    if args.nk:
        nk = complex(args.nk.replace("i", "j")); nk_src = "user (--nk)"
    elif "nk" in meta and any(abs(e - args.pump_eV) < 0.08 for e in meta["nk"]):
        e = [e for e in meta["nk"] if abs(e - args.pump_eV) < 0.08][0]
        nk = meta["nk"][e]; nk_src = "built-in table: " + meta.get("nk_note", "")
    else:
        nk = nk_tb; nk_src = "universal tight-binding dielectric function (interband, scissor corrected)"
    if args.absorbed_fluence_mJcm2 is not None:
        F_abs = args.absorbed_fluence_mJcm2; Aabs = F_abs / args.fluence_mJcm2; R_ = T_ = None
        depth_m = (args.thickness_nm or 0) * 1e-9
    elif args.thickness_nm is not None:
        Aabs, R_, T_ = film_absorption(nk, args.thickness_nm * 1e-9, lam)
        F_abs = args.fluence_mJcm2 * Aabs
        depth_m = args.thickness_nm * 1e-9
    else:
        r = (1 - nk) / (1 + nk); R_ = float(abs(r) ** 2); Aabs = 1 - R_; T_ = 0.0
        F_abs = args.fluence_mJcm2 * Aabs
        depth_m = lam / (4 * np.pi * max(nk.imag, 1e-6))   # surface energy density F_abs * alpha
    if depth_m <= 0:
        raise SystemExit("--absorbed-fluence-mJcm2 requires --thickness-nm")
    E_abs_cell = (F_abs * 1e-3 / 1e-4) / depth_m * (V_cell * 1e-30) / EV_J
    n_exc = E_abs_cell / args.pump_eV
    frac_e = 1.0 if el["metal"] else max(args.pump_eV - el["eg"], 0.0) / args.pump_eV
    log(f"[optics] n = {nk:.3f} ({nk_src}); TB estimate {nk_tb if nk_tb is None else np.round(nk_tb, 3)}; "
        f"absorptance {Aabs:.3f}; absorbed {F_abs:.2f} mJ/cm^2 -> {E_abs_cell*1e3:.1f} meV/cell, "
        f"{n_exc:.4f} e-h pairs/cell")
    el["tabs"] = electron_tables(el, n_exc, args.T)
    t0 = time.time()
    W, ebins, einfo = stage_eph(outdir, ph3, tb, el, mesh, freq, qfrac, args, pdims, n_exc, args.T)
    timings["electron_phonon"] = time.time() - t0
    if args.stop_after == "eph":
        return
    Tgrid = np.geomspace(max(args.T * 0.5, 5.0), 50000, 48)
    hw = freq * THz2eV
    Atab = A_from_W(W, hw, ebins, el, el["tabs"], Tgrid)
    Atab[freq <= args.fmin_thz] = 0.0
    Atab *= args.eph_scale
    # mode-resolved coupling summary: lambda-like strength from the high-Te limit
    procs = [(P[f"a{k}"], P[f"b{k}"], P[f"c{k}"], P[f"v{k}"] if _PHPH_KERNEL is not None else P[f"v{k}"].astype(float)) for k in range(3)]
    kin = Kinetics(freq, procs, Atab, Tgrid, el["tabs"], args, E_abs_cell, frac_e, el["metal"], args.fmin_thz)
    from scipy.integrate import solve_ivp
    n0 = np.where(kin.valid, nB_thz(freq, args.T), 0.0)
    y0 = np.concatenate([[0.0, 1.0], n0.reshape(-1)])
    tmax = args.tmax_ps
    t_start = -3 * args.fwhm_fs * 1e-3 - 0.05
    tout = np.unique(np.concatenate([np.linspace(t_start, min(2.0, tmax), 120),
                                     np.linspace(min(2.0, tmax), min(10.0, tmax), 120),
                                     np.geomspace(min(10.0, tmax), tmax, 120)]))
    tout = tout[(tout >= t_start) & (tout <= tmax)]
    A300 = Atab[:, :, np.searchsorted(Tgrid, 3000)]
    log(f"[dyn] {Nq} q x {nb} branches, {sum(len(p[3]) for p in procs)} three-phonon processes; "
        f"max A(3000 K) = {A300.max():.2f}/ps; integrating to {tmax} ps")
    t0 = time.time()
    # explicit RK45 when the electron-phonon rates are moderate (graphite), implicit BDF with the
    # analytic electron-phonon Jacobian when they make the system stiff
    stiff = A300.max() > 20.0 or os.environ.get("RM_STIFF") == "1"
    dyn_path = os.path.join(outdir, "dynamics" + args.tag + ".npz")
    dkey = _key(dict(P=str(P["key"]), A=float(Atab.sum()), E=E_abs_cell, fe=frac_e, T=args.T, tmax=tmax,
                     fwhm=args.fwhm_fs, rec=args.recomb_ps))
    if args.reuse_dynamics and _cached(dyn_path, dkey):
        z = np.load(dyn_path)
        t, Ee, n = z["t"], z["Ee"], z["n"].astype(float)
        nfev = int(z["nfev"])
        log("[dyn] re-using cached dynamics")
    else:
        if stiff:
            sol = solve_ivp(kin.rhs, (t_start, tmax), y0, method="BDF", t_eval=tout, rtol=1e-4, atol=1e-7,
                            jac=kin.jac, first_step=1e-3)
        else:
            sol = solve_ivp(kin.rhs, (t_start, tmax), y0, method="RK45", t_eval=tout, rtol=1e-4, atol=1e-7)
        timings["dynamics"] = time.time() - t0
        log(f"[dyn] {sol.message} nfev={sol.nfev} ({timings['dynamics']:.0f}s)")
        t = sol.t
        Ee = sol.y[0]; n = sol.y[2:].T.reshape(len(t), Nq, nb)
        nfev = int(sol.nfev)
    Te = np.array([kin.Te_of(e) for e in Ee])
    Eph = (n * hw[None]).sum((1, 2)) / Nq
    E_in = E_abs_cell * frac_e
    energy_check = dict(absorbed_to_carriers_eV_cell=float(E_in),
                        final_carrier_plus_phonon_excess_eV_cell=float(Ee[-1] + Eph[-1] - Eph[0]))
    # ---------------- observables
    res = observables(t, n, Te, freq, qfrac, mesh, ph, prim, prim_cell, pdims, atoms, args, kin)
    res.update(dict(
        material=name, structure_source=meta.get("source"), formula=atoms.get_chemical_formula(),
        n_atoms_primitive=nat, periodic_dims=pdims.tolist(),
        conditions=dict(pump_eV=args.pump_eV, incident_fluence_mJcm2=args.fluence_mJcm2, fwhm_fs=args.fwhm_fs,
                        T_K=args.T, thickness_nm=args.thickness_nm, tmax_ps=tmax),
        optics=dict(n_complex=[nk.real, nk.imag], source=nk_src,
                    n_complex_tb=None if nk_tb is None else [float(nk_tb.real), float(nk_tb.imag)],
                    absorptance=Aabs, reflectance=R_, transmittance=T_, absorbed_fluence_mJcm2=F_abs,
                    absorbed_energy_meV_per_cell=E_abs_cell * 1e3, energy_density_J_cm3=E_abs_cell * EV_J / (V_cell * 1e-24),
                    eh_pairs_per_cell=n_exc, eh_density_cm3=n_exc / (V_cell * 1e-24),
                    fraction_to_carrier_heat=frac_e),
        electrons=dict(model="universal sp3 Slater-Koster (Harrison 1980) + NIST LDA term values",
                       metal=el["metal"], gap_tb_eV=el["eg_tb"], gap_used_eV=el["eg"], scissor_eV=tb.scissor,
                       valence_electrons=tb.nval, Te_peak_K=float(Te.max()), t_Te_peak_ps=float(t[np.argmax(Te)])),
        phonons=dict(model="MACE medium-omat-0", mesh=list(mesh), fc2_supercell=fcinfo["sc2"], fc3_supercell=fcinfo["sc3"],
                     fc3_pair_cutoff_A=fcinfo["cutoff"], three_phonon_width_THz=sig_dyn,
                     calibration_spread=float(P["cal_spread"]), n_modes_excluded=n_invalid,
                     gamma_THz=np.round(fG[0], 3).tolist()),
        numerics=dict(nk_eph=list(args.nk_eph), nk_dos=list(args.nk_dos), sig_e_eV=args.sig_e, nsub=args.nsub,
                      e_window_eV=args.e_window, eph_info=einfo, ode_nfev=nfev,
                      integrator="BDF (analytic e-ph Jacobian)" if stiff else "RK45"),
        tag=args.tag, eph_scale=args.eph_scale,
        energy_check=energy_check, warnings=warns, timings_s=timings, relax=relax_info,
    ))
    res["timings_s"]["total"] = time.time() - T_start
    np.savez_compressed(dyn_path, key=dkey, nfev=nfev, t=t, n=n.astype(np.float32), Te=Te, Ee=Ee,
                        freq=freq, qfrac=qfrac, E_abs_cell=E_abs_cell)
    np.savez_compressed(os.path.join(outdir, "eph_rates" + args.tag + ".npz"), qfrac=qfrac, freq=freq, Te=Tgrid, A=Atab.astype(np.float32))
    with open(os.path.join(outdir, "results" + args.tag + ".json"), "w") as f:
        json.dump(_jsonable(res), f, indent=1)
    make_figure(res, outdir, name)
    log(f"[done] results in {outdir} (total {res['timings_s']['total']:.0f}s)")
    return res


def default_kmesh(tb, pdims, target):
    rec = np.linalg.norm(tb.rec, axis=1)
    d = int(pdims.sum())
    scale = (target / np.prod(rec[pdims])) ** (1.0 / d)
    return [max(4, int(round(scale * rec[i]))) if pdims[i] else 1 for i in range(3)]


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, complex):
        return [o.real, o.imag]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def observables(t, n, Te, freq, qfrac, mesh, ph, prim, prim_cell, pdims, atoms, args, kin):
    import spglib
    Nq, nb = freq.shape
    hw = freq * THz2eV
    valid = freq > args.fmin_thz
    dn = n - n[0][None]
    irf = (args.irf_fs if args.irf_fs is not None else args.fwhm_fs) * 1e-3
    probe = math.sqrt(max(irf ** 2 - (args.fwhm_fs * 1e-3) ** 2, 0.0))
    out = {"fit_settings": dict(irf_fwhm_ps=irf, probe_convolution_fwhm_ps=probe,
                                note="simulated traces (pump-limited) are convolved with the probe duration and fitted "
                                     "with model (x) Gaussian IRF")}
    # branch-resolved populations and effective temperatures
    ac = np.arange(nb) < 3
    E_br = (n * hw[None]).sum(1) / Nq                        # (t, nb) eV per cell
    Tb = np.zeros((len(t), nb))
    for j in range(nb):
        hwj = hw[valid[:, j], j]
        Tb[:, j] = eff_temperature(hwj, (n[:, valid[:, j], j] * hwj[None]).sum(1)) if len(hwj) else np.nan
    hw_op = hw[:, ~ac][valid[:, ~ac]]; hw_ac = hw[:, ac][valid[:, ac]]
    T_op = eff_temperature(hw_op, (n[:, :, ~ac][:, valid[:, ~ac]] * hw_op).sum(1))
    T_ac = eff_temperature(hw_ac, (n[:, :, ac][:, valid[:, ac]] * hw_ac).sum(1))
    nop = np.array([dn[i][:, ~ac][valid[:, ~ac]].mean() for i in range(len(t))])
    nac = np.array([dn[i][:, ac][valid[:, ac]].mean() for i in range(len(t))])
    # hot-phonon set: top 5 % of modes by excess occupation at the electron-temperature peak + 0.3 ps
    it_h = int(np.argmin(abs(t - (t[np.argmax(Te)] + 0.3))))
    score = np.where(valid, dn[it_h], -np.inf).reshape(-1)
    nh = max(1, int(0.05 * valid.sum()))
    hot = np.argsort(score)[::-1][:nh]
    hot_trace = dn.reshape(len(t), -1)[:, hot].mean(1)
    hq, hj = np.unravel_index(hot, (Nq, nb))
    r63, tpk = rise_time(t, hot_trace)
    times = dict(
        Te_peak_K=float(Te.max()),
        Te_cooling_1e_ps=decay_time(t, Te),
        hot_phonon_rise_63_ps=r63, hot_phonon_peak_ps=tpk, hot_phonon_decay_1e_ps=decay_time(t, hot_trace),
        hot_phonon_fit_rise_decay=fit_trace(t, hot_trace, "rise_decay", 0.0, 0.0, (t[0], min(t[-1], 30.0))),
        optical_T_rise_63_ps=rise_time(t, T_op)[0], optical_T_peak_ps=rise_time(t, T_op)[1],
        acoustic_T_rise_63_ps=rise_time(t, T_ac)[0],
        acoustic_T_90_ps=_t_frac(t, T_ac, 0.9),
    )
    spread = np.nanmax(Tb, 1) - np.nanmin(Tb, 1)
    k = np.argmax(spread)
    below = np.where(spread[k:] < 0.1 * spread[k])[0]
    times["lattice_equilibration_ps"] = float(t[k + below[0]]) if len(below) else None
    times["lattice_equilibration_definition"] = ("time after which the spread of branch effective temperatures "
                                                 "falls below 10% of its maximum")
    times["hot_set"] = dict(n_modes=int(nh), branches=np.bincount(hj, minlength=nb).tolist(),
                            mean_freq_THz=float(freq[hq, hj].mean()))
    out["characteristic_times"] = times
    sel_t = [x for x in (0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 50, 100, 200) if x <= t[-1]]
    idx = [int(np.argmin(abs(t - x))) for x in sel_t]
    out["branches"] = dict(
        times_ps=sel_t, Te_K=Te[idx].tolist(), T_optical_K=T_op[idx].tolist(), T_acoustic_K=T_ac[idx].tolist(),
        T_branch_K=Tb[idx].tolist(), excess_energy_meV_per_cell=((E_br - E_br[0]) * 1e3)[idx].tolist(),
        mean_excess_occupation_optical=nop[idx].tolist(), mean_excess_occupation_acoustic=nac[idx].tolist(),
        branch_order="branches sorted by frequency at each q; 0-2 acoustic")
    out["traces"] = dict(t_ps=t.tolist(), Te_K=Te.tolist(), T_optical_K=T_op.tolist(), T_acoustic_K=T_ac.tolist(),
                         T_branch_K=Tb.tolist(), dn_optical=nop.tolist(), dn_acoustic=nac.tolist(),
                         dn_hot=hot_trace.tolist())
    # ---------------- eigenvectors on the grid
    fq, ev = phonons_at(ph, qfrac)
    rec = 2 * np.pi * np.linalg.inv(prim_cell).T
    rec_in = 2 * np.pi * np.linalg.inv(np.array(atoms.cell)).T
    masses = np.array(prim.masses)
    syms = list(prim.symbols)
    tau = np.array(prim.scaled_positions) @ prim_cell
    nat = len(masses)
    # Debye-Waller: per-atom mean-square displacement tensors U_s(t) (A^2)
    w_rad = 2 * np.pi * np.maximum(freq, 1e-6) * 1e12
    e3 = ev.reshape(Nq, nat, 3, nb)
    pref = HBAR_Js / (2 * masses[None, :, None] * AMU) * 1e20 / Nq     # (1,nat,1)
    U = np.zeros((len(t), nat, 3, 3))
    outer = np.einsum("qsan,qsbn->qsabn", e3, e3.conj()).real            # (q,s,a,b,nu)
    for it in range(len(t)):
        wgt = np.where(valid, (2 * n[it] + 1) / w_rad, 0.0)               # (q, nu)
        U[it] = np.einsum("qsabn,qn->sab", outer, wgt) * pref[0, :, 0][:, None, None]
    dU = U - U[0][None]
    out["debye_waller"] = dict(U0_A2=U[0].tolist(), dU_avg_A2_at_times=[dU[i].mean(0).tolist() for i in idx],
                               times_ps=sel_t)
    if args.bragg:
        bl = []
        for s in args.bragg.split(";"):
            hkl = np.array(parse_triplet(s))
            G = hkl @ rec_in
            Gh = G / np.linalg.norm(G)
            f = scattering_factor_list(syms, np.linalg.norm(G))
            def inten(Ut):
                return abs(np.sum(f * np.exp(-1j * tau @ G) * np.exp(-0.5 * np.einsum("a,sab,b->s", G, Ut, G)))) ** 2
            I = np.array([inten(U[i]) for i in range(len(t))])
            lnI = -np.log(I / I[0])
            dUG = np.einsum("a,tsab,b->ts", Gh, dU, Gh).mean(1)
            bl.append(dict(hkl=hkl.tolist(), G_invA=float(np.linalg.norm(G)),
                           neg_ln_I_over_I0=lnI[idx].tolist(), dU_along_G_A2=dUG[idx].tolist(),
                           fit_dU_biexp=fit_trace(t, dUG, "biexp", irf, probe),
                           fit_neg_ln_I_biexp=fit_trace(t, lnI, "biexp", irf, probe),
                           trace_dU_along_G_A2=dUG.tolist()))
            log(f"[obs] Bragg {hkl}: dU(G-direction) at {sel_t[-1]} ps = {dUG[idx[-1]]:.2e} A^2; "
                f"biexp {bl[-1]['fit_dU_biexp']}")
        out["bragg"] = bl
    # ---------------- one-phonon diffuse scattering at high-symmetry points
    ds = spglib.get_symmetry_dataset((prim_cell, np.array(prim.scaled_positions), np.array(prim.numbers)))
    A_ = prim_cell.T
    rots = [A_ @ R @ np.linalg.inv(A_) for R in ds.rotations]
    sp = special_points(prim_cell, pdims)
    beam = None
    if args.beam:
        beam = np.array(parse_triplet(args.beam)) @ np.array(atoms.cell)
        beam /= np.linalg.norm(beam)
    elif not pdims.all():
        ax = int(np.where(~pdims)[0][0])
        beam = np.cross(prim_cell[(ax + 1) % 3], prim_cell[(ax + 2) % 3]); beam /= np.linalg.norm(beam)
    zoneG = None
    if args.diffuse_zone:
        zoneG = np.array(parse_triplet(args.diffuse_zone)) @ rec_in
    rng = np.arange(-12, 13)
    Gs = []
    pr = [rng if pdims[i] else np.array([0]) for i in range(3)]
    for h in pr[0]:
        for k_ in pr[1]:
            for l_ in pr[2]:
                Gs.append(np.array([h, k_, l_]) @ rec)
    Gs = np.array(Gs)
    Gs = Gs[np.linalg.norm(Gs, axis=1) < args.qmax + 2]
    diff = {}
    for lab, qf_ in sp.items():
        qc0 = qf_ @ rec
        pts = []
        for qs in star_cart(qc0, rots):
            gi, qg, err = nearest_grid(qs, rec, mesh, pdims)
            # use the exact wavevector at which the eigenvectors were computed: with phonopy's phase
            # convention e(q + G') differs from e(q) by exp(-i G'.r_s), so G must be Q - q_used
            qg = qfrac[gi] @ rec
            for G in Gs:
                Q = qg + G
                Qn = np.linalg.norm(Q)
                if Qn < 0.3 or Qn > args.qmax:
                    continue
                if beam is not None and abs(Q @ beam) > 0.05:
                    continue
                # G must be a reciprocal-lattice vector of the primitive cell for the phase convention
                Gf = np.linalg.solve(rec.T, Q - qg)
                if not np.allclose(Gf, np.round(Gf), atol=1e-6):
                    continue
                pts.append((gi, Q, Q - qg, err))
        if pts and zoneG is not None:
            # only the equivalent points closest to the chosen Bragg peak (the zone around it)
            dz = np.array([np.linalg.norm(p[1] - zoneG) for p in pts])
            pts = [p for p, dd in zip(pts, dz) if dd <= dz.min() + 1e-3]
        if not pts:
            continue
        Wm = np.zeros((len(pts), nb))
        for kp, (gi, Q, G, err) in enumerate(pts):
            f = scattering_factor_list(syms, np.linalg.norm(Q))
            for j in range(nb):
                if not valid[gi, j]:
                    continue
                e = ev[gi][:, j].reshape(nat, 3)
                F = np.sum(f / np.sqrt(masses) * (e @ Q) * np.exp(-1j * tau @ G))
                Wm[kp, j] = abs(F) ** 2 / freq[gi, j]
        gis = np.array([p[0] for p in pts])
        I = np.array([np.sum((n[it][gis] + 0.5) * Wm) for it in range(len(t))])
        rel = I / I[0] - 1
        lab2 = "Gamma" if lab == "G" else lab
        if lab == "G":
            continue
        diff[lab2] = dict(q_frac_primitive=np.round(qf_, 4).tolist(), n_Q_points=len(pts),
                          example_Q_invA=np.round(pts[0][1], 3).tolist(), grid_offset_invA=float(max(p[3] for p in pts)),
                          dI_over_I0_at_times=rel[idx].tolist(),
                          fit_rise_decay=fit_trace(t, rel, "rise_decay", irf, probe),
                          fit_biexp=fit_trace(t, rel, "biexp", irf, probe),
                          fit_single=fit_trace(t, rel, "single", irf, probe),
                          rise_63_ps=rise_time(t, rel)[0], peak_ps=rise_time(t, rel)[1],
                          branch_weights=(Wm.sum(0) / max(Wm.sum(), 1e-30)).round(3).tolist(),
                          trace=rel.tolist())
        log(f"[obs] diffuse {lab2}: {len(pts)} Q points, rise63 {diff[lab2]['rise_63_ps']}, "
            f"peak {diff[lab2]['peak_ps']}, rise_decay fit {diff[lab2]['fit_rise_decay']}")
    out["diffuse"] = diff
    out["diffuse_definition"] = ("one-phonon intensity I1(Q) = sum_nu (n+1/2)/omega |sum_s f_s(|Q|)/sqrt(M_s) (Q.e_s) "
                                 "exp(-i G.r_s)|^2 at Q = q + G for all symmetry-equivalent q of the special point "
                                 "(nearest grid point) and all G with |Q| < qmax in the plane normal to the beam; "
                                 "dI/I0 = I(t)/I(t<0) - 1")
    return out


def scattering_factor_list(syms, Q):
    return np.array([scattering_factor(s, Q) for s in syms])


def _t_frac(t, y, frac):
    y = y - y[0]
    if y[-1] <= 0:
        return None
    k = np.where(y >= frac * y[-1])[0]
    return float(t[k[0]]) if len(k) else None


def make_figure(res, outdir, name):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tr = res["traces"]
    t = np.array(tr["t_ps"])
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(t, tr["Te_K"], "k", label="electrons $T_e$")
    ax[0].plot(t, tr["T_optical_K"], "C3", label="optical branches")
    ax[0].plot(t, tr["T_acoustic_K"], "C0", label="acoustic branches")
    Tb = np.array(tr["T_branch_K"])
    for j in range(Tb.shape[1]):
        ax[0].plot(t, Tb[:, j], lw=0.5, color="0.6")
    ax[0].set_xscale("symlog", linthresh=1); ax[0].set_yscale("log"); ax[0].set_xlabel("time (ps)")
    ax[0].set_ylabel("effective temperature (K)"); ax[0].legend(fontsize=8); ax[0].set_title(f"{name}: temperatures")
    ax[1].plot(t, tr["dn_hot"], "C3", label="hot-phonon set (top 5%)")
    ax[1].plot(t, tr["dn_optical"], "C1", label="optical, mean")
    ax[1].plot(t, tr["dn_acoustic"], "C0", label="acoustic, mean")
    ax[1].set_xscale("symlog", linthresh=1); ax[1].set_xlabel("time (ps)"); ax[1].set_ylabel("excess occupation")
    ax[1].legend(fontsize=8); ax[1].set_title("phonon populations")
    for k, (lab, d) in enumerate(res.get("diffuse", {}).items()):
        ax[2].plot(t, d["trace"], label=lab)
    for b in res.get("bragg", []):
        pass
    ax[2].set_xscale("symlog", linthresh=1); ax[2].set_xlabel("time (ps)"); ax[2].set_ylabel(r"$\Delta I/I_0$ (one-phonon)")
    ax[2].legend(fontsize=7, ncol=2); ax[2].set_title("diffuse intensity at high-symmetry points")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "results" + res.get("tag", "") + ".png"), dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
