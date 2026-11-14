"""Optical conductivity of the half-filled open Hubbard chain from phase-space
trajectories: Gaussian fluctuations around the Hartree-Fock antiferromagnet, each
trajectory evolved with its own TDHF Hamiltonian; kicked minus unkicked runs with the same
samples give the linear response including fluctuation (beyond-RPA) effects.
Usage: python chain_sigma_twa.py U [L] [nsamp]
"""
import sys, os, time
import numpy as np
from scipy.special import j1
from scipy.integrate import quad
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lightmat.ftwa import HubbardTWA

U = float(sys.argv[1]); L = int(sys.argv[2]) if len(sys.argv) > 2 else 32
ns = int(sys.argv[3]) if len(sys.argv) > 3 else 400
gap = U - 4 + 8 * quad(lambda w: j1(w) / (w * (1 + np.exp(w * U / 2))), 0, 200, limit=400)[0]
h = np.zeros((L, L))
for i in range(L - 1):
    h[i, i + 1] = h[i + 1, i] = -1.0
stag = np.array([(-1) ** i for i in range(L)], float)
nu = 0.5 + 0.4 * stag; nd = 1 - nu
for it in range(3000):
    eu, vu = np.linalg.eigh(h + U * np.diag(nd)); ed, vd = np.linalg.eigh(h + U * np.diag(nu))
    ru = vu[:, :L // 2] @ vu[:, :L // 2].T; rd = vd[:, :L // 2] @ vd[:, :L // 2].T
    if max(abs(np.diag(ru) - nu).max(), abs(np.diag(rd) - nd).max()) < 1e-11:
        break
    nu, nd = 0.5 * nu + 0.5 * np.diag(ru), 0.5 * nd + 0.5 * np.diag(rd)
A0 = 2e-3
occ = np.r_[np.ones(L // 2), np.zeros(L - L // 2)]
W = np.zeros((2 * L, 2 * L)); W[:L, :L] = vu; W[L:, L:] = vd
res = []
t0 = time.time()
for kick in (0.0, A0):
    s = HubbardTWA(np.eye(L), U, occ, occ, basis="spin", exchange=1.0, nsamp=ns, seed=11)
    s.rho = np.einsum("ab,sbc,dc->sad", W, s.rho, W)
    hk = h.astype(complex).copy()
    for i in range(L - 1):
        hk[i, i + 1] = -np.exp(-1j * kick); hk[i + 1, i] = -np.exp(1j * kick)
    s.hfun = lambda t, hk=hk: hk
    dt, t, jj = 0.02, 0.0, []
    while t < 60 - 1e-9:
        r = s.rho.mean(0)
        j = 0.0
        for i in range(L - 1):
            for sh in (0, L):
                j += -2 * np.imag(np.exp(-1j * kick) * r[i + 1 + sh, i + sh]) * (-1)
        jj.append(j); s.step(t, dt); t += dt
    res.append(np.array(jj))
tt = np.arange(len(res[0])) * dt
j = res[1] - res[0]
w = np.linspace(0.02, 1.5 * U + 6, 700)
sig = np.array([-np.sum(j * np.exp((1j * x - 0.1) * tt)) * dt for x in w]) / (A0 * L)
re = sig.real
thr = w[np.argmax(re > 0.1 * re.max())]
np.savez(os.path.join(os.path.dirname(__file__), "data", f"chain_sigma_twa_U{U:g}.npz"), w=w, sigma=sig, gap=gap)
print(f"U={U} TWA: Bethe-ansatz gap {gap:.3f}; onset(10%) {thr:.2f}, peak {w[np.argmax(re)]:.2f}; {time.time()-t0:.0f}s", flush=True)
