"""Checks of the pp-ladder solver against exact dynamics on tiny chains."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from lightmat.ed import Hubbard, chain_bonds, evolve
from lightmat.fluct_hubbard import PPLadder


def run(L, up, dn, U, tmax=4.0, dt=0.005, mode="factorized"):
    m = Hubbard(L, len(up), len(dn), chain_bonds(L), U)
    psi = m.product_state(up, dn)
    times = np.arange(0, tmax + 1e-9, 0.5)
    ex = evolve(m, psi, times, observe=lambda p, t: (m.densities(p), m.double_occ(p)))
    h = np.zeros((L, L))
    for i, j in chain_bonds(L):
        h[i, j] = h[j, i] = -1.0
    ru = np.diag([1.0 if i in up else 0.0 for i in range(L)])
    rd = np.diag([1.0 if i in dn else 0.0 for i in range(L)])
    s = PPLadder(h, U, ru, rd, mode=mode)
    t = 0.0; out = []
    nsteps = int(round(0.5 / dt))
    for k, T in enumerate(times):
        nu, nd = s.densities()
        out.append((nu + nd, nu * nd + s.dcorr()))
        for _ in range(nsteps):
            s.step(t, dt); t += dt
    err_n = max(np.abs(o[0] - (e[0][0] + e[0][1])).max() for o, e in zip(out, ex))
    err_d = max(np.abs(o[1] - e[1]).max() for o, e in zip(out, ex))
    return err_n, err_d


if __name__ == "__main__":
    print("1up1dn, U=4 (two-particle limit):", run(5, [1], [3], 4.0))
    print("1up1dn, U=4, G1-G2 mode       :", run(5, [1], [3], 4.0, mode="g1g2"))
    print("CDW L=6 U=0.5 (weak)          :", run(6, [0, 2, 4], [0, 2, 4], 0.5, tmax=2))
    print("CDW L=6 U=0.5 G1-G2           :", run(6, [0, 2, 4], [0, 2, 4], 0.5, tmax=2, mode="g1g2"))
