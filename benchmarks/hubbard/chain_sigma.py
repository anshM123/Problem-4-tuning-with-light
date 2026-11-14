"""Optical conductivity of the half-filled open Hubbard chain (cf. Fig. 10 of
Schroedter et al.).  Start: self-consistent Hartree-Fock antiferromagnet at U.  Two runs
with identical fluctuation samples, one with a weak vector-potential step A0 at t = 0
(E = -A0 delta(t)) and one without; their current difference is the linear response
including the correlation dynamics.  sigma(w) = -j(w)/(A0 L), Lorentzian eta = 0.1 J.
Compared with the exact Bethe-ansatz Mott gap of the infinite chain:
  Delta(U) = U - 4J + 8J int_0^inf J1(w) / (w (1 + exp(w U / 2J))) dw.
Usage: python chain_sigma.py U [L] [nsamp] [kernel] [collide 0/1]
"""
import sys, os, time
import numpy as np
from scipy.special import j1
from scipy.integrate import quad
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lightmat.fluct_ph import PHFluct

U = float(sys.argv[1]); L = int(sys.argv[2]) if len(sys.argv) > 2 else 32
ns = (None if sys.argv[3] == "exact" else int(sys.argv[3])) if len(sys.argv) > 3 else 200
kernel = sys.argv[4] if len(sys.argv) > 4 else "gwx"
collide = bool(int(sys.argv[5])) if len(sys.argv) > 5 else True
gap = U - 4 + 8 * quad(lambda w: j1(w) / (w * (1 + np.exp(w * U / 2))), 0, 200, limit=400)[0]
A0 = 2e-3
h = np.zeros((L, L))
for i in range(L - 1):
    h[i, i + 1] = h[i + 1, i] = -1.0
stag = np.array([(-1) ** i for i in range(L)], float)
nu = 0.5 + 0.4 * stag; nd = 1 - nu
for it in range(2000):                                   # HF antiferromagnet
    eu, vu = np.linalg.eigh(h + U * np.diag(nd)); ed, vd = np.linalg.eigh(h + U * np.diag(nu))
    ru = vu[:, :L // 2] @ vu[:, :L // 2].T; rd = vd[:, :L // 2] @ vd[:, :L // 2].T
    nu2, nd2 = np.diag(ru).copy(), np.diag(rd).copy()
    if max(abs(nu2 - nu).max(), abs(nd2 - nd).max()) < 1e-10:
        break
    nu, nd = 0.5 * nu + 0.5 * nu2, 0.5 * nd + 0.5 * nd2
gap_hf = (eu[L // 2] - eu[L // 2 - 1])
rho0 = np.zeros((2 * L, 2 * L)); rho0[:L, :L] = ru; rho0[L:, L:] = rd
bonds = [(i, i + 1, 1.0, -1.0) for i in range(L - 1)]
t0 = time.time()
runs = []
for kick in (0.0, A0):
    s = PHFluct(L, bonds, U, rho0, kernel=kernel, nsamp=ns, seed=7, peierls=lambda t, k=kick: k)
    s.collide = collide
    dt, t, jj, tt = 0.02, 0.0, [], []
    while t < 60.0 - 1e-9:
        tt.append(t); jj.append(s.current(t))
        s.step(t, dt); t += dt
    runs.append(np.array(jj))
tt = np.array(tt); j = runs[1] - runs[0]
w = np.linspace(0.02, 1.5 * U + 6, 700)
eta = 0.1
sig = np.array([-np.sum(j * np.exp((1j * x - eta) * tt)) * dt for x in w]) / (A0 * L)
re = sig.real
thr = w[np.argmax(re > 0.1 * re.max())]
tag = f"U{U:g}_{kernel}_c{int(collide)}"
np.savez(os.path.join(os.path.dirname(__file__), "data", f"chain_sigma_{tag}.npz"),
         w=w, sigma=sig, gap=gap, gap_hf=gap_hf, t=tt, j=j)
print(f"U={U} [{tag}]: Bethe-ansatz gap {gap:.3f}; HF gap {gap_hf:.3f}; sigma onset(10%) {thr:.2f}, "
      f"peak {w[np.argmax(re)]:.2f}; {time.time()-t0:.0f}s", flush=True)
