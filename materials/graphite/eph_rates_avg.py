"""Brillouin-zone cell averaging of the electron-phonon generation coefficients near
Gamma, K and K' (where the emission is concentrated in a region smaller than one grid
cell).  For grid points within 2 grid spacings of these points, A_{nu q} is replaced by its
average over NSUB x NSUB points uniformly covering the cell.  Branches are matched by
frequency rank.  Writes eph_rates.npz (cell-averaged); eph_rates_point.npz keeps the
point-sampled table.
"""
import os, sys, time
import numpy as np, h5py
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from eph_tb import bands, coupling, frac_to_cart, B, HBAR_eVs, KB_eV, THz2eV

HERE = os.path.dirname(os.path.abspath(__file__))
NSUB = int(sys.argv[1]) if len(sys.argv) > 1 else 4
Nk = int(sys.argv[2]) if len(sys.argv) > 2 else 450
P = np.load(os.path.join(HERE, "eph_rates_point.npz"))
qf, freq, Te, A = P["qfrac"], P["freq"], P["Te"], P["A"].copy()
Nq1 = int(round(np.sqrt(len(qf))))
SIG = 0.03; BETA2 = 0.0994 / 0.00454
a_ = 2.464
cell = [[a_, 0, 0], [-a_ / 2, a_ * np.sqrt(3) / 2, 0], [0, 0, 15.0]]
ph = Phonopy(PhonopyAtoms(symbols=["C", "C"], cell=cell, scaled_positions=[[1/3, 2/3, .5], [2/3, 1/3, .5]]),
             supercell_matrix=[8, 8, 1])
ph.force_constants = h5py.File(os.path.join(HERE, "fc2.hdf5"))["force_constants"][:]
f = (np.arange(Nk) + 0.5) / Nk
kc = np.stack(np.meshgrid(f, f, indexing="ij"), -1).reshape(-1, 2) @ B
ek, _ = bands(kc)
kT = KB_eV * Te
M2 = np.array([[1, 0], [-0.5, np.sqrt(3) / 2]])


def near(q, q0, r):
    best = 1e9
    for s1 in (-1, 0, 1):
        for s2 in (-1, 0, 1):
            d = (q - np.array(q0) - [s1, s2]) @ M2
            best = min(best, np.linalg.norm(d))
    return best <= r + 1e-9


sel = [i for i, q in enumerate(qf) if any(near(q, c, 2.0 / Nq1) for c in ([0, 0], [1/3, 1/3], [2/3, 2/3]))]
offs = (np.arange(NSUB) + 0.5) / NSUB - 0.5
t0 = time.time()
for cnt, iq in enumerate(sel):
    subs = np.array([qf[iq] + np.array([o1, o2]) / Nq1 for o1 in offs for o2 in offs])
    ph.run_qpoints(np.c_[subs, np.zeros(len(subs))], with_eigenvectors=True)
    fs, evs = ph.qpoints.frequencies, ph.qpoints.eigenvectors
    acc = np.zeros((6, len(Te)))
    for isub, q in enumerate(subs):
        qc = frac_to_cart(q)
        ekq, _ = bands(kc + qc)
        for nu in range(6):
            hw = fs[isub, nu] * THz2eV
            if hw < 1e-3:
                continue
            g = coupling(kc, qc, evs[isub][:, nu].reshape(2, 3), fs[isub, nu], 1.0)
            g2 = BETA2 * np.abs(g) ** 2
            for m in range(2):
                for n in range(2):
                    de = ekq[:, m] - ek[:, n] - hw
                    s_ = np.abs(de) < 4 * SIG
                    if not s_.any():
                        continue
                    w = g2[s_, m, n] * np.exp(-0.5 * (de[s_] / SIG) ** 2) / (np.sqrt(2 * np.pi) * SIG)
                    en = ek[s_, n][:, None] / kT[None]; em = ekq[s_, m][:, None] / kT[None]
                    occ = 1 / (np.exp(np.clip(en, -700, 700)) + 1) * (1 - 1 / (np.exp(np.clip(em, -700, 700)) + 1))
                    acc[nu] += (w[:, None] * occ).sum(0)
    acc *= (2 * np.pi / HBAR_eVs) * 2 / len(kc) * 1e-12 / len(subs)
    ratio = acc[5, 15] / max(A[iq, 5, 15], 1e-30)
    A[iq] = acc
    print(f"{cnt+1}/{len(sel)} q={qf[iq]} A_top(avg/point)={ratio:.2f} {time.time()-t0:.0f}s", flush=True)
np.savez(os.path.join(HERE, "eph_rates.npz"), qfrac=qf, freq=freq, Te=Te, A=A)
print("done")
