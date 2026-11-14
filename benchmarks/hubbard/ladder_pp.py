"""Ladder ground-state spin correlations from the factorized pp-ladder method with
adiabatic switching.  <dSz_i dSz_j> = 1/4 [ C^uu_ij + C^dd_ij - C^ud_ij - C^du_ij ] with
C^ss_ij = n_i delta_ij - |rho_ij|^2 (same spin, exchange hole) and C^ud_ij = G_{ij,ij}
(pair correlation of the T-matrix ladder)."""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lightmat.fluct_hubbard import PPLadder
from lightmat.ed import ladder_bonds

Lx, L = 6, 12
Tr = float(sys.argv[1]) if len(sys.argv) > 1 else 40.0
legs, rungs = ladder_bonds(Lx)
h = np.zeros((L, L))
for i, j in legs: h[i, j] = h[j, i] = -1.0
for i, j in rungs: h[i, j] = h[j, i] = -0.5
e, v = np.linalg.eigh(h); r0 = v[:, :6] @ v[:, :6].T
rows = []
for U in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.5]:
    t0 = time.time()
    s = PPLadder(h, U, r0, r0, Ufun=lambda t, U=U: U * np.sin(0.5 * np.pi * min(t / Tr, 1)) ** 2)
    t, dt, acc = 0.0, 0.01, []
    while t < Tr + 10 - 1e-9:
        s.step(t, dt); t += dt
        if t > Tr and int(round(t / dt)) % 20 == 0:
            ru, rd = s.rhos(); G = s.pair_corr()
            Cuu = np.diag(np.real(np.diag(ru))) - np.abs(ru) ** 2
            Cdd = np.diag(np.real(np.diag(rd))) - np.abs(rd) ** 2
            SS = 0.25 * (Cuu + Cdd - G - G.T)
            loc = np.mean(np.diag(SS))
            leg = -np.mean([SS[x + Lx * r, (x + 1) % Lx + Lx * r] for x in range(Lx) for r in (0, 1)])
            rung = -np.mean([SS[x, x + Lx] for x in range(Lx)])
            acc.append((loc, leg, rung))
    c = np.mean(acc, 0); rows.append((U,) + tuple(c))
    print(f"pp-ladder U={U}: C=({c[0]:.4f},{c[1]:.4f},{c[2]:.4f}) {time.time()-t0:.0f}s", flush=True)
np.savetxt(os.path.join(os.path.dirname(__file__), "data", "ladder_gs_ppladder.txt"), rows, header="U C_loc C_leg C_rung")
