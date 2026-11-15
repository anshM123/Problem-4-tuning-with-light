"""Anharmonic energy surface of LiNbO3 along the driven A1(TO4) coordinate (MACE).

Normal coordinates Q (Å amu^1/2) of the zone-centre modes: atomic displacement
u_k = Q e_k / sqrt(m_k), so V = 1/2 w^2 Q^2 + ... with V in eV per primitive cell.

1. Frozen-phonon scan along Q_IR = Q(TO4), fitted with the polynomial of von Hoegen et al.
   (V = 1/2 w^2 Q^2 + a3 Q^3/3 + a4 Q^4/4 + a5 Q^5/5, |Q| <= 1.5) for comparison with their DFT
   coefficients, plus an 8th-order fit used for the dynamics.
2. Two-dimensional scans (Q_IR, Q_j) for the other three A1 modes j = TO1, TO2, TO3 (the only
   modes that couple linearly to Q_IR^2 by symmetry), fitted with all coupling terms
   Q_IR^m Q_j^n, m, n >= 1, m + n <= 6, and the one-dimensional potentials V_j(Q_j).
3. Check: the cubic couplings g_j (coefficient of Q_IR^2 Q_j) and a3 against the phono3py
   cubic force constants summed at Gamma (the same fc3 used for the three-phonon rates).
4. Adiabatic potential: V along Q_IR with the other three A1 coordinates relaxed.
5. Mode effective charges: DFPT polarities (params.py) with the sign fixed for the MACE
   eigenvectors; nominal-charge estimate printed for comparison.

Writes anharmonic.npz and anharmonic.log.  Run time about 6 min.
"""
import os, json
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
import torch
torch.set_num_threads(4)
from ase import Atoms
from scipy.optimize import minimize
from mace.calculators import mace_mp
import h5py
from params import ZMODE, ZZZ_SIGN

HERE = os.path.dirname(os.path.abspath(__file__))
calc = mace_mp(model="medium-omat-0", default_dtype="float64", device="cpu")
g = np.load(os.path.join(HERE, "gamma_modes.npz"))
freq, ev, labels, masses = g["freq"], g["eigvec"], g["labels"], g["masses"]
nat = len(masses)
atoms0 = Atoms(numbers=g["numbers"], cell=g["cell"], scaled_positions=g["scaled_positions"], pbc=True)
atoms0.calc = calc
E0 = atoms0.get_potential_energy()
pos0 = atoms0.get_positions().copy()
sym = atoms0.get_chemical_symbols()
sqm = np.sqrt(np.repeat(masses, 3))
W2 = (2 * np.pi * freq * 1e12) ** 2 * (1e-20 * 1.66053906660e-27) / 1.602176634e-19  # eV/(Å^2 amu)

# ---- identify the A1 modes: polar along z with non-zero z displacement of Li and Nb
evr = ev.T.reshape(30, nat, 3)            # mode, atom, xyz (mass-weighted eigenvectors)
ia1 = [b for b in range(3, 30) if str(labels[b]).startswith("A1")]
if len(ia1) != 4:  # fall back: invariant under the 3-fold axis => Li/Nb move only along z
    lin = [k for k, s in enumerate(sym) if s in ("Li", "Nb")]
    score = [np.abs(evr[b, lin, 2]).sum() / (np.abs(evr[b, lin, :2]).sum() + 1e-9) for b in range(30)]
    ia1 = sorted([b for b in range(3, 30) if score[b] > 50 and np.abs(evr[b, lin, 2]).sum() > 1e-3])
assert len(ia1) == 4, (ia1, labels)
ia1 = sorted(ia1, key=lambda b: freq[b])
iIR = ia1[-1]
out = open(os.path.join(HERE, "anharmonic.log"), "w")


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); out.write(s + "\n")


log("A1 modes (MACE):", [f"{freq[b]:.3f} THz = {freq[b]*33.356:.0f} cm-1" for b in ia1])
log("measured A1 TO (Raman/IR): 252, 275, 332, 632 cm-1; DFPT-LDA 243, 287, 355, 617")

# ---- signs: mode polarity with nominal charges -> align with DFPT polarity (positive)
qnom = np.array([ZZZ_SIGN[s] for s in sym])
zsign = []
evec = {}
for k, b in enumerate(ia1):
    e = evr[b].copy()
    pnom = np.sum(qnom * e[:, 2] / np.sqrt(masses))   # e / sqrt(amu)
    if pnom < 0:
        e = -e
    evec[k] = e
    log(f"  A1 TO{k+1}: nominal-charge polarity {abs(pnom):.3f} e/amu^1/2, DFPT {ZMODE[k]:.3f}")


def energy(disp):
    a = atoms0.copy(); a.calc = calc
    a.set_positions(pos0 + disp)
    return a.get_potential_energy() - E0


def disp_of(coeffs):
    d = np.zeros((nat, 3))
    for k, q in coeffs.items():
        d += q * evec[k] / np.sqrt(masses)[:, None]
    return d


# ---- 1. scan along Q_IR
qs = np.round(np.arange(-3.0, 3.0001, 0.05), 4)
V1 = np.array([energy(disp_of({3: q})) for q in qs])
# sign convention for Q_IR: soft side at negative Q (as in von Hoegen et al., Fig. 5)
c3 = np.polyfit(qs[np.abs(qs) <= 0.5], V1[np.abs(qs) <= 0.5], 3)[0]
sgn_ir = 1.0
if c3 < 0:
    sgn_ir = -1.0; evec[3] = -evec[3]; qs_ = -qs[::-1]; V1 = V1[::-1]
log(f"Q_IR orientation relative to the polar-charge convention: {sgn_ir:+.0f} "
    f"(polarity sign of the driven mode = {sgn_ir:+.0f})")
m = np.abs(qs) <= 1.5
A = np.vstack([qs[m] ** 2 / 2, qs[m] ** 3 / 3, qs[m] ** 4 / 4, qs[m] ** 5 / 5]).T
w2fit, a3, a4, a5 = np.linalg.lstsq(A, V1[m], rcond=None)[0]
nu_fit = np.sqrt(w2fit * 1.602176634e-19 / (1e-20 * 1.66053906660e-27)) / (2 * np.pi) / 1e12
log(f"Q_IR potential fit |Q|<=1.5 (paper form): w = {nu_fit:.2f} THz (harmonic {freq[iIR]:.2f}), "
    f"a3 = {a3*1e3:.0f}, a4 = {a4*1e3:.0f}, a5 = {a5*1e3:.0f} meV  [DFT of von Hoegen et al.: "
    f"a3 = 1568, a4 = 901, a5 = 7.1 meV, units amu^-n/2 Å^-n]")
p8 = np.polyfit(qs, V1, 10)   # name kept: high-order fit used by the dynamics
log(f"10th-order fit rms {np.sqrt(np.mean((np.polyval(p8, qs)-V1)**2))*1e3:.2f} meV over |Q|<=3")

# ---- 2. 2D scans with the other A1 modes
q4g = np.round(np.arange(-2.0, 2.0001, 0.25), 3)
qjg = np.round(np.arange(-2.5, 2.5001, 0.5), 3)
terms = [(mm, nn) for mm in range(1, 6) for nn in range(1, 6) if mm + nn <= 6]
coup = {}
selfpot = {}
for k in range(3):
    # 1D potential of mode j
    qq = np.round(np.arange(-3.0, 3.0001, 0.25), 3)
    Vj = np.array([energy(disp_of({k: q})) for q in qq])
    selfpot[k] = np.polyfit(qq, Vj, 8)
    # coupling energy surface
    rows, rhs = [], []
    for a in q4g:
        for b in qjg:
            if a == 0 or b == 0:
                continue
            Eab = energy(disp_of({3: a, k: b}))
            Ec = Eab - np.polyval(p8, a) - np.polyval(selfpot[k], b)
            rows.append([a ** mm * b ** nn for mm, nn in terms]); rhs.append(Ec)
    cf, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
    res = np.array(rows) @ cf - np.array(rhs)
    coup[k] = cf
    g21 = cf[terms.index((2, 1))]
    w2j = np.polyfit(qq[np.abs(qq) <= 0.5], Vj[np.abs(qq) <= 0.5], 2)[0] * 2
    log(f"TO{k+1}: frozen-phonon freq {np.sqrt(max(w2j,0)*1.602176634e-19/(1e-20*1.66053906660e-27))/(2*np.pi)/1e12:.2f} THz; "
        f"g (Q_IR^2 Q_j) = {g21*1e3:.1f} meV; Q_IR Q_j^2: {cf[terms.index((1,2))]*1e3:.1f}; "
        f"Q_IR^2 Q_j^2: {cf[terms.index((2,2))]*1e3:.1f} meV; fit rms {np.sqrt(np.mean(res**2))*1e3:.2f} meV")

# ---- 3. cubic couplings from phono3py fc3 at Gamma
HAVE_FC3 = os.path.exists(os.path.join(HERE, "fc3.hdf5"))
if HAVE_FC3:  # noqa
    with h5py.File(os.path.join(HERE, "fc3.hdf5"), "r") as f:
        fc3 = f["fc3"][:]
    import phono3py
    ph3 = phono3py.load(os.path.join(HERE, "phono3py_params.yaml"), produce_fc=False, log_level=0)
    s2p = ph3.supercell.s2p_map if hasattr(ph3.supercell, "s2p_map") else ph3.primitive.s2p_map
    p2s = ph3.primitive.p2s_map
    s2pp = np.array([list(p2s).index(x) for x in s2p])      # supercell atom -> primitive index
    if fc3.shape[0] != nat:  # full fc3
        fc3 = fc3[p2s]
    Phi = np.zeros((nat, nat, nat, 3, 3, 3))
    for jb in range(fc3.shape[1]):
        for jc in range(fc3.shape[2]):
            Phi[:, s2pp[jb], s2pp[jc]] += fc3[:, jb, jc]
    # map primitive ordering of phono3py onto ours (same cell, same order expected)
    def cubic(e1, e2, e3):
        x1 = e1 / np.sqrt(masses)[:, None]; x2 = e2 / np.sqrt(masses)[:, None]; x3 = e3 / np.sqrt(masses)[:, None]
        return np.einsum("abcijk,ai,bj,ck->", Phi, x1, x2, x3)
    a3_fc3 = 0.5 * cubic(evec[3], evec[3], evec[3])
    log(f"fc3 check: a3 (=Phi_444/2) = {a3_fc3*1e3:.0f} meV vs frozen phonon {a3*1e3:.0f}")
    for k in range(3):
        log(f"  g_TO{k+1} from fc3 = {0.5*cubic(evec[3], evec[3], evec[k])*1e3:.1f} meV vs frozen phonon "
            f"{coup[k][terms.index((2,1))]*1e3:.1f}")

# ---- 4. adiabatic (relaxed) potential along Q_IR: the three other A1 coordinates (the only
# ones a fully symmetric distortion couples to) relaxed at every Q_IR
others = [0, 1, 2]
basis = np.array([evec[k] for k in others])          # mass-weighted unit vectors


# continuation from Q_IR = 0 outward keeps the relaxed coordinates on the connected branch
qad = np.round(np.arange(-1.6, 1.61, 0.2), 3)
Vad = np.zeros(len(qad)); cad = np.zeros((len(qad), len(others)))
i0 = int(np.argmin(np.abs(qad)))
for rng_ in (range(i0, len(qad)), range(i0, -1, -1)):
    c_prev = np.zeros(len(others))
    for i in rng_:
        d0 = disp_of({3: qad[i]})
        def fun(c):
            d = d0 + np.tensordot(c, basis, 1) / np.sqrt(masses)[:, None]
            a = atoms0.copy(); a.calc = calc; a.set_positions(pos0 + d)
            f = a.get_forces()
            return a.get_potential_energy() - E0, -np.einsum("kai,ai->k", basis, f / np.sqrt(masses)[:, None])
        r = minimize(fun, c_prev, jac=True, method="L-BFGS-B",
                     bounds=[(cp - 0.5, cp + 0.5) for cp in c_prev], options=dict(gtol=1e-5))
        Vad[i], cad[i] = r.fun, r.x
        c_prev = r.x
log("adiabatic potential (Q, V_frozen, V_relaxed eV; relaxed TO1, TO2, TO3 coordinates):")
for q, v, c in zip(qad, Vad, cad):
    log(f"  {q:+.2f} {np.polyval(p8, q):.4f} {v:.4f}   {c[0]:+.3f} {c[1]:+.3f} {c[2]:+.3f}")

# ---- direction of the spontaneous polarization: Nb is off-centre toward the O3 face with
# the three short Nb-O bonds, i.e. along +P_S
D = atoms0.get_all_distances(mic=True, vector=True)
dzs = []
for k, s_ in enumerate(sym):
    if s_ != "Nb":
        continue
    io = [i for i, t in enumerate(sym) if t == "O"]
    d = np.linalg.norm(D[k, io], axis=1)
    near = np.array(io)[np.argsort(d)[:3]]
    dzs.append(D[k, near, 2].mean())
sgnP = float(np.sign(np.mean(dzs)))
log(f"P_S along {'+' if sgnP > 0 else '-'}z (mean dz of the three short Nb-O bonds {np.mean(dzs):+.3f} Å)")

np.savez(os.path.join(HERE, "anharmonic.npz"), sgnP=sgnP, qs=qs, V1=V1, p8=p8, a3=a3, a4=a4, a5=a5, w2fit=w2fit,
         terms=np.array(terms), coup=np.array([coup[k] for k in range(3)]),
         selfpot=np.array([selfpot[k] for k in range(3)]), ia1=np.array(ia1), freqA1=freq[ia1],
         evecA1=np.array([evec[k] for k in range(4)]), sgn_ir=sgn_ir, qad=qad, Vad=Vad,
         zmode=ZMODE * np.array([1, 1, 1, sgn_ir]))
out.close()
