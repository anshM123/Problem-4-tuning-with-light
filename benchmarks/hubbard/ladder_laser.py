"""Laser-driven 6x2 Hubbard ladder at U/J = 1 (protocol of Schroedter et al., Fig. 9):
A(t) = A_L cos(w_L (t - t_L)) exp(-(t - t_L)^2 / 2 s^2), A_L = 0.5, w_L = 2, t_L = 10, s = 2,
Peierls phase on the leg bonds.  Observable dC_i(t) = C_i(t) - C_i(0), C as in ladder_gs.py.
Usage: python ladder_laser.py ed|fast
"""
import sys, os, time
import numpy as np
from scipy.sparse.linalg import expm_multiply
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lightmat.ed import Hubbard, ladder_bonds
from lightmat.fluct_ph import PHFluct

Lx, L, U = 6, 12, 1.0
AL, wL, tL, sg = 0.5, 2.0, 10.0, 2.0
A = lambda t: AL * np.cos(wL * (t - tL)) * np.exp(-(t - tL) ** 2 / (2 * sg ** 2))
legs, rungs = ladder_bonds(Lx, periodic=True)
Jmap = {b: 1.0 for b in legs}; Jmap.update({b: 0.5 for b in rungs})
tout = np.arange(0, 100.0001, 0.1)
D = os.path.join(os.path.dirname(__file__), "data")


def corr(SS):
    loc = np.mean([SS[j, j] for j in range(L)])
    leg = -np.mean([SS[x + Lx * r, (x + 1) % Lx + Lx * r] for x in range(Lx) for r in (0, 1)])
    rung = -np.mean([SS[x, x + Lx] for x in range(Lx)])
    return loc, leg, rung


what = sys.argv[1] if len(sys.argv) > 1 else "ed"
t0 = time.time()
if what == "ed":
    m = Hubbard(L, Lx, Lx, legs + rungs, U, J=Jmap)
    E, psi = m.ground_state()
    def obs(p):
        SS = m.szsz(p); sz = m.sz(p)
        return corr(SS - np.outer(sz, sz))
    Hs = m.H()
    out = [obs(psi)]
    for k in range(1, len(tout)):
        ta, tb = tout[k - 1], tout[k]
        if tb < tL + 6 * sg + 0.2:
            nsub = 4; dt = (tb - ta) / nsub
            for s in range(nsub):
                tm = ta + (s + 0.5) * dt
                ph = {b: (np.exp(1j * A(tm)) if b in legs else 1.0) for b in legs + rungs}
                psi = expm_multiply(-1j * m.H(ph) * dt, psi)
        else:
            psi = expm_multiply(-1j * Hs * (tb - ta), psi)
        out.append(obs(psi))
        if k % 100 == 0:
            print(f"t={tb:.1f} {time.time()-t0:.0f}s", flush=True)
    np.savez(os.path.join(D, "ladder_laser_ed.npz"), t=tout, C=np.array(out))
elif what == "pp":
    from lightmat.fluct_hubbard import PPLadder
    Tr = 40.0
    h0 = np.zeros((L, L), complex)
    for (i, j) in legs: h0[i, j] = h0[j, i] = -1.0
    for (i, j) in rungs: h0[i, j] = h0[j, i] = -0.5
    e, v = np.linalg.eigh(h0.real); r0 = v[:, :Lx] @ v[:, :Lx].T
    def hfun(t):
        a = A(t - Tr) if t > Tr else 0.0
        hh = h0.copy()
        for (i, j) in legs:      # Peierls phase on c^+_i c_j along the leg
            hh[i, j] = -np.exp(1j * a); hh[j, i] = -np.exp(-1j * a)
        return hh
    s = PPLadder(hfun, U, r0, r0, Ufun=lambda t: U * np.sin(0.5 * np.pi * min(t / Tr, 1)) ** 2)
    def cpp():
        ru, rd = s.rhos(); G = s.pair_corr()
        Cuu = np.diag(np.real(np.diag(ru))) - np.abs(ru) ** 2
        Cdd = np.diag(np.real(np.diag(rd))) - np.abs(rd) ** 2
        return corr(0.25 * (Cuu + Cdd - G - G.T))
    dt, t = 0.01, 0.0
    while t < Tr - 1e-9:
        s.step(t, dt); t += dt
    out = [cpp()]
    for k in range(1, len(tout)):
        while t < Tr + tout[k] - 1e-9:
            s.step(t, dt); t += dt
        out.append(cpp())
    np.savez(os.path.join(D, "ladder_laser_pp.npz"), t=tout, C=np.array(out))
else:
    kernel = sys.argv[2] if len(sys.argv) > 2 else "gwx"
    Tr = 30.0
    h = np.zeros((L, L)); bonds = []
    for (i, j) in legs:
        bonds.append((i, j, 1.0, -1.0)); h[i, j] = h[j, i] = -1.0
    for (i, j) in rungs:
        bonds.append((i, j, 0.5, 0.0)); h[i, j] = h[j, i] = -0.5
    e, v = np.linalg.eigh(h); r1 = v[:, :Lx] @ v[:, :Lx].T
    rho0 = np.zeros((2 * L, 2 * L)); rho0[:L, :L] = r1; rho0[L:, L:] = r1
    s = PHFluct(L, bonds, U, rho0, kernel=kernel,
                Ufun=lambda t: U * np.sin(0.5 * np.pi * min(t / Tr, 1)) ** 2,
                peierls=lambda t: A(t - Tr) if t > Tr else 0.0)
    dt, t = 0.02, 0.0
    while t < Tr - 1e-9:
        s.step(t, dt); t += dt
    out = [corr(s.szsz())]
    for k in range(1, len(tout)):
        while t < Tr + tout[k] - 1e-9:
            s.step(t, dt); t += dt
        out.append(corr(s.szsz()))
    np.savez(os.path.join(D, f"ladder_laser_fast_{kernel}.npz"), t=tout, C=np.array(out))
print("done", time.time() - t0)
