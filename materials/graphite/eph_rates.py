"""Tabulate phonon generation coefficients A_{nu q}(Te) on the Nq x Nq phonon grid.
Output: eph_rates.npz with qfrac (Nq^2,2), freq (Nq^2,6) [THz], Te (nT,), A (Nq^2,6,nT) [1/ps].
"""
import os, sys, time
import numpy as np, h5py
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from eph_tb import (bands, coupling, frac_to_cart, B, HBAR_eVs, KB_eV, THz2eV)

HERE = os.path.dirname(os.path.abspath(__file__))
Nq = int(sys.argv[1]) if len(sys.argv) > 1 else 24
Nk = int(sys.argv[2]) if len(sys.argv) > 2 else 450
SIG = 0.03                    # eV Gaussian broadening of the energy delta
BETA2 = 0.0994 / 0.00454      # beta^2 (eV/A)^2 fixed by <g^2_K> = 0.0994 eV^2
Te = np.geomspace(300, 40000, 32)

a_ = 2.464
cell = [[a_, 0, 0], [-a_ / 2, a_ * np.sqrt(3) / 2, 0], [0, 0, 15.0]]
unit = PhonopyAtoms(symbols=["C", "C"], cell=cell,
                    scaled_positions=[[1 / 3, 2 / 3, .5], [2 / 3, 1 / 3, .5]])
ph = Phonopy(unit, supercell_matrix=[8, 8, 1])
ph.force_constants = h5py.File(os.path.join(HERE, "fc2.hdf5"))["force_constants"][:]
qf = np.array([[i / Nq, j / Nq] for i in range(Nq) for j in range(Nq)])
ph.run_qpoints(np.c_[qf, np.zeros(len(qf))], with_eigenvectors=True)
freq = ph.qpoints.frequencies; evec = ph.qpoints.eigenvectors

f = (np.arange(Nk) + 0.5) / Nk
kc = np.stack(np.meshgrid(f, f, indexing="ij"), -1).reshape(-1, 2) @ B
ek, _ = bands(kc)
A = np.zeros((len(qf), 6, len(Te)))
t0 = time.time()
kT = KB_eV * Te
for iq, q in enumerate(qf):
    qc = frac_to_cart(q)
    ekq, _ = bands(kc + qc)
    for nu in range(6):
        hw = freq[iq, nu] * THz2eV
        if hw < 1e-3:
            continue
        g = coupling(kc, qc, evec[iq][:, nu].reshape(2, 3), freq[iq, nu], 1.0)
        g2 = BETA2 * np.abs(g) ** 2                        # (Nk, m, n)
        for m in range(2):
            for n in range(2):
                de = ekq[:, m] - ek[:, n] - hw
                sel = np.abs(de) < 4 * SIG
                if not sel.any():
                    continue
                w = g2[sel, m, n] * np.exp(-0.5 * (de[sel] / SIG) ** 2) / (np.sqrt(2 * np.pi) * SIG)
                en = ek[sel, n][:, None] / kT[None]; em = ekq[sel, m][:, None] / kT[None]
                occ = 1 / (np.exp(np.clip(en, -700, 700)) + 1) * (1 - 1 / (np.exp(np.clip(em, -700, 700)) + 1))
                A[iq, nu] += (w[:, None] * occ).sum(0)
    if iq % 48 == 0:
        print(f"q {iq}/{len(qf)} {time.time()-t0:.0f}s", flush=True)
A *= (2 * np.pi / HBAR_eVs) * 2 / len(kc) * 1e-12     # 1/ps
np.savez(os.path.join(HERE, "eph_rates.npz"), qfrac=qf, freq=freq, Te=Te, A=A)
print("done", time.time() - t0)
