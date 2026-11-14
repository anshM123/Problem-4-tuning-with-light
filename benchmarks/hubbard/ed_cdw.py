"""Exact reference: CDW melting in an 11-site open Hubbard chain.
Initial state: doubly occupied odd sites (1,3,5,7,9), empty even sites (N=10).
By particle-hole symmetry of the bipartite chain this is equivalent to doublons on the
six even sites (N=12): I(t) and d_corr(t) are identical.
Outputs I(t), d(t), d_corr(t) (site-resolved and totals) to data/ed_cdw_U{U}.npz.
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lightmat.ed import Hubbard, chain_bonds, evolve

L = 11
occ = [1, 3, 5, 7, 9]
emp = [i for i in range(L) if i not in occ]
Us = [float(u) for u in sys.argv[1:]] or [5.0, 15.0]
tmax, dt = 50.0, 0.1
times = np.arange(0, tmax + 1e-9, dt)
os.makedirs(os.path.join(os.path.dirname(__file__), "data"), exist_ok=True)
for U in Us:
    t0 = time.time()
    m = Hubbard(L, len(occ), len(occ), chain_bonds(L), U)
    psi0 = m.product_state(occ, occ)

    def obs(psi, t):
        nu, nd = m.densities(psi)
        d = m.double_occ(psi)
        return dict(nu=nu, nd=nd, d=d)

    res = evolve(m, psi0, times, observe=obs)
    nu = np.array([r["nu"] for r in res]); nd = np.array([r["nd"] for r in res])
    d = np.array([r["d"] for r in res])
    n = nu + nd
    I = n[:, occ].mean(1) - n[:, emp].mean(1)
    dcorr = d - nu * nd
    np.savez(os.path.join(os.path.dirname(__file__), "data", f"ed_cdw_U{U:g}.npz"),
             t=times, n=n, d=d, dcorr=dcorr, I=I)
    print(f"U={U}: done in {time.time()-t0:.0f}s; I(50)={I[-1]:.3f} "
          f"sum dcorr(50)={dcorr[-1].sum():.3f} mean d(50)={d[-1].mean():.4f}", flush=True)
