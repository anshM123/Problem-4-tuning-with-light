"""Ladder ground-state spin correlations with the fermionic phase-space trajectory method
(adiabatic switching of U from the noninteracting ground state).
Usage: python ladder_twa.py basis exchange nsamp Tramp
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lightmat.ftwa import HubbardTWA
from lightmat.ed import ladder_bonds

basis = sys.argv[1] if len(sys.argv) > 1 else "spin"
x = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
ns = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
Tr = float(sys.argv[4]) if len(sys.argv) > 4 else 20.0
Lx = 6; L = 12
legs, rungs = ladder_bonds(Lx)
h = np.zeros((L, L))
for i, j in legs:
    h[i, j] = h[j, i] = -1.0
for i, j in rungs:
    h[i, j] = h[j, i] = -0.5
e, v = np.linalg.eigh(h)
occ = np.r_[np.ones(6), np.zeros(6)]
rows = []
for U in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.5]:
    t0 = time.time()
    s = HubbardTWA(np.diag(e), 0.0, occ, occ, basis=basis, exchange=x, nsamp=ns, seed=1)
    W = np.zeros((2 * L, 2 * L)); W[:L, :L] = v; W[L:, L:] = v
    s.rho = np.einsum("ab,sbc,dc->sad", W, s.rho, W)
    s.hfun = lambda t: h
    t, dt, acc = 0.0, 0.02, []
    while t < Tr + 10:
        s.U = U * np.sin(0.5 * np.pi * min(t / Tr, 1)) ** 2
        s.step(t, dt); t += dt
        if t > Tr:
            a = np.real(s.rho[:, np.arange(L), np.arange(L)])
            b = np.real(s.rho[:, np.arange(L) + L, np.arange(L) + L])
            if basis == "nambu":
                b = 1 - b
            sz = 0.5 * (a - b); sz = sz - sz.mean(0)
            SS = sz.T @ sz / len(sz)
            loc = np.mean(np.diag(SS))
            leg = -np.mean([SS[xx + Lx * r, (xx + 1) % Lx + Lx * r] for xx in range(Lx) for r in (0, 1)])
            rung = -np.mean([SS[xx, xx + Lx] for xx in range(Lx)])
            acc.append((loc, leg, rung))
    c = np.mean(acc, 0)
    rows.append((U,) + tuple(c))
    print(f"TWA[{basis},x={x}] U={U}: C=({c[0]:.4f},{c[1]:.4f},{c[2]:.4f}) {time.time()-t0:.0f}s", flush=True)
np.savetxt(os.path.join(os.path.dirname(__file__), "data", f"ladder_gs_twa_{basis}_x{x:g}.txt"), rows,
           header="U C_loc C_leg C_rung")
