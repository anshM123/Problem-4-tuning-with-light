"""Energy-resolved electron-phonon transition kernel for non-thermal electrons in graphite.

W_{q nu}(i) = (2 pi / hbar) (2 / N_k) sum_{k, m, n} |g^nu_mn(k, q)|^2 delta(e_m(k+q) - e_n(k) - hw)
              [e_n(k) in energy bin i]
counts transitions that lift an electron from bin i by one phonon energy (absorption; the
reverse is emission).  With any electron distribution f,
  dn_q/dt = sum_i W(i) [ f(e_i + hw)(1 - f(e_i))(n+1) - f(e_i)(1 - f(e_i + hw)) n ],
and for a Fermi-Dirac f this reduces exactly to A_q(Te)[(n+1)e^{-hw/kTe} - n] (eph_rates.py),
which is checked at the end.  Same couplings, GW scaling (applied later) and cell averaging
near Gamma, K, K' as the thermal tables.  Output: eph_kernel.npz (W[Nq, 6, Nbins] in 1/ps).
"""
import os, sys, time
import numpy as np, h5py
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from eph_tb import bands, coupling, frac_to_cart, B, HBAR_eVs, THz2eV

HERE = os.path.dirname(os.path.abspath(__file__))
Nk = 450; SIG = 0.03; BETA2 = 0.0994 / 0.00454; NSUB = 4
DE, EMAX = 0.002, 3.0
edges = np.arange(-EMAX, EMAX + DE / 2, DE); centers = 0.5 * (edges[1:] + edges[:-1]); NB = len(centers)
P = np.load(os.path.join(HERE, "eph_rates_point.npz"))
qf, freq = P["qfrac"], P["freq"]; Nq1 = int(round(np.sqrt(len(qf))))
a_ = 2.464
cell = [[a_, 0, 0], [-a_ / 2, a_ * np.sqrt(3) / 2, 0], [0, 0, 15.0]]
f_ = (np.arange(Nk) + 0.5) / Nk
kc = np.stack(np.meshgrid(f_, f_, indexing="ij"), -1).reshape(-1, 2) @ B
ek, _ = bands(kc)
M2 = np.array([[1, 0], [-0.5, np.sqrt(3) / 2]])
offs = (np.arange(NSUB) + 0.5) / NSUB - 0.5
_ph = None


def phonons():
    global _ph
    if _ph is None:
        _ph = Phonopy(PhonopyAtoms(symbols=["C", "C"], cell=cell, scaled_positions=[[1/3, 2/3, .5], [2/3, 1/3, .5]]),
                      supercell_matrix=[8, 8, 1])
        _ph.force_constants = h5py.File(os.path.join(HERE, "fc2.hdf5"))["force_constants"][:]
    return _ph


def near(q, q0, r):
    return min(np.linalg.norm((q - np.array(q0) - [s1, s2]) @ M2) for s1 in (-1, 0, 1) for s2 in (-1, 0, 1)) <= r + 1e-9


def task(iq):
    q = qf[iq]
    if any(near(q, c, 2.0 / Nq1) for c in ([0, 0], [1/3, 1/3], [2/3, 2/3])):
        subs = np.array([q + np.array([o1, o2]) / Nq1 for o1 in offs for o2 in offs])
    else:
        subs = q[None]
    ph = phonons()
    ph.run_qpoints(np.c_[subs, np.zeros(len(subs))], with_eigenvectors=True)
    fs, evs = ph.qpoints.frequencies, ph.qpoints.eigenvectors
    W = np.zeros((6, NB))
    for isub, qq in enumerate(subs):
        qc = frac_to_cart(qq); ekq, _ = bands(kc + qc)
        for nu in range(6):
            hw = fs[isub, nu] * THz2eV
            if hw < 1e-3:
                continue
            g2 = BETA2 * np.abs(coupling(kc, qc, evs[isub][:, nu].reshape(2, 3), fs[isub, nu], 1.0)) ** 2
            for m in range(2):
                for n in range(2):
                    de = ekq[:, m] - ek[:, n] - hw
                    s_ = np.abs(de) < 4 * SIG
                    if not s_.any():
                        continue
                    w = g2[s_, m, n] * np.exp(-0.5 * (de[s_] / SIG) ** 2) / (np.sqrt(2 * np.pi) * SIG)
                    ib = np.clip(((ek[s_, n] + EMAX) / DE).astype(int), 0, NB - 1)
                    W[nu] += np.bincount(ib, weights=w, minlength=NB)
    return W * (2 * np.pi / HBAR_eVs) * 2 / len(kc) * 1e-12 / len(subs)


if __name__ == "__main__":
    from multiprocessing import Pool
    t0 = time.time()
    with Pool(int(os.environ.get("NPROC", "20"))) as pool:
        W = np.array(pool.map(task, range(len(qf)), chunksize=1))
    np.savez(os.path.join(HERE, "eph_kernel.npz"), W=W, centers=centers, DE=DE, qfrac=qf, freq=freq)
    # consistency check with the thermal table
    T = np.load(os.path.join(HERE, "eph_rates.npz"))
    kB = 8.617333262e-5
    for Te in (3000.0, 10000.0):
        it = np.argmin(abs(T["Te"] - Te)); Tg = T["Te"][it]
        fe = 1 / (np.exp(centers / (kB * Tg)) + 1)
        A = np.zeros((len(qf), 6))
        for j in range(6):
            s = np.rint(freq[:, j] * THz2eV / DE).astype(int)
            for iq in range(len(qf)):
                ii = np.arange(NB - s[iq]) if s[iq] > 0 else np.arange(0)
                A[iq, j] = np.sum(W[iq, j, ii] * fe[ii] * (1 - fe[ii + s[iq]])) if s[iq] > 0 else 0
        rel = np.abs(A - T["A"][:, :, it]).max() / T["A"][:, :, it].max()
        print(f"check vs thermal table at Te = {Tg:.0f} K: max relative deviation {rel:.3f}")
    print("done", time.time() - t0)
