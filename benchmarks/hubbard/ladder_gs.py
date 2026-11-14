"""Ground-state spin correlations of the half-filled periodic 6x2 Hubbard ladder
(J_leg = 1, J_rung = 0.5): exact diagonalization vs the particle-hole fluctuation method
with adiabatic switching of U.
C_loc = <dSz_j dSz_j>, C_leg = -<dSz_j dSz_{j+1}>, C_rung = -<dSz_j dSz_{j+Lx}> (site-averaged).
Usage: python ladder_gs.py [ed|fast] [kernel] [Tramp]
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lightmat.ed import Hubbard, ladder_bonds
from lightmat.fluct_ph import PHFluct

Lx = 6; L = 2 * Lx
legs, rungs = ladder_bonds(Lx, periodic=True)
Jmap = {b: 1.0 for b in legs}; Jmap.update({b: 0.5 for b in rungs})
Us = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.5]
what = sys.argv[1] if len(sys.argv) > 1 else "ed"
kernel = sys.argv[2] if len(sys.argv) > 2 else "gwx"
Tramp = float(sys.argv[3]) if len(sys.argv) > 3 else 40.0
out = os.path.join(os.path.dirname(__file__), "data")


def corr(SS):
    loc = np.mean([SS[j, j] for j in range(L)])
    leg = -np.mean([SS[x + Lx * r, (x + 1) % Lx + Lx * r] for x in range(Lx) for r in (0, 1)])
    rung = -np.mean([SS[x, x + Lx] for x in range(Lx)])
    return loc, leg, rung


rows = []
if what == "ed":
    for U in Us:
        t0 = time.time()
        m = Hubbard(L, Lx, Lx, legs + rungs, U, J=Jmap)
        E, psi = m.ground_state()
        SS = m.szsz(psi); sz = m.sz(psi)
        SS = SS - np.outer(sz, sz)
        rows.append((U,) + corr(SS))
        print(f"ED U={U}: E={E:.5f} C=({rows[-1][1]:.4f},{rows[-1][2]:.4f},{rows[-1][3]:.4f}) {time.time()-t0:.0f}s", flush=True)
    np.savetxt(os.path.join(out, "ladder_gs_ed.txt"), rows, header="U C_loc C_leg C_rung")
else:
    # noninteracting ground state
    h = np.zeros((L, L))
    bonds = []
    for (i, j) in legs:
        dx = -1.0 if (j == (i + 1) % Lx + Lx * (i // Lx)) else 1.0
        bonds.append((i, j, 1.0, dx)); h[i, j] = h[j, i] = -1.0
    for (i, j) in rungs:
        bonds.append((i, j, 0.5, 0.0)); h[i, j] = h[j, i] = -0.5
    e, v = np.linalg.eigh(h)
    r1 = v[:, :Lx] @ v[:, :Lx].T
    rho0 = np.zeros((2 * L, 2 * L)); rho0[:L, :L] = r1; rho0[L:, L:] = r1
    for U in Us:
        t0 = time.time()
        ramp = lambda t, U=U: U * (np.sin(0.5 * np.pi * min(t / Tramp, 1.0)) ** 2)
        s = PHFluct(L, bonds, U, rho0, kernel=kernel, Ufun=ramp)
        dt = 0.02; t = 0.0
        acc = []
        while t < Tramp + 10 - 1e-9:
            s.step(t, dt); t += dt
            if t > Tramp:
                acc.append(corr(s.szsz()))
        c = np.mean(acc, axis=0)
        rows.append((U,) + tuple(c))
        print(f"fast[{kernel}] U={U}: C=({c[0]:.4f},{c[1]:.4f},{c[2]:.4f}) S={s.S} {time.time()-t0:.0f}s", flush=True)
    np.savetxt(os.path.join(out, f"ladder_gs_fast_{kernel}.txt"), rows, header="U C_loc C_leg C_rung")
