"""(b) Level A limit: the unified solver with phonons switched off (Q = 0) against
lightmat/fluct_hubbard.PPLadder on the 11-site CDW of benchmarks/hubbard (U/J = 5, same
protocol as fast_cdw.py: t = 0..50, dt = 0.002), for the T-matrix ('factorized') and the
self-interaction-free ('sic') pair ladders.  Also compared with the stored Level A results
benchmarks/hubbard/data/fast_cdw_U5_<mode>.npz.

Usage: python levelA_limit.py [U] [tmax] [dt]   -> data/levelA_limit.npz
"""
import os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, ROOT)
from lightmat.unified import Model, UnifiedSolver
from lightmat.fluct_hubbard import PPLadder

L = 11
occ = [1, 3, 5, 7, 9]
emp = [i for i in range(L) if i not in occ]


def run(U, mode, tmax, dt, which, tout=0.1):
    h = np.zeros((L, L))
    for i in range(L - 1):
        h[i, i + 1] = h[i + 1, i] = -1.0
    r0 = np.diag([1.0 if i in occ else 0.0 for i in range(L)])
    if which == "levelA":
        s = PPLadder(h, U, r0, r0, mode=mode)
        dens = s.densities; dc = s.dcorr
    else:
        s = UnifiedSolver(Model(h, U=U), r0, r0, ee=mode)
        dens = s.densities; dc = lambda: s._dcorr(s.G, s.x, s.y)
    times = np.arange(0, tmax + 1e-9, tout)
    nsub = int(round(tout / dt))
    I, D = [], []
    t = 0.0
    for T in times:
        nu, nd = dens()
        n = nu + nd
        I.append(n[occ].mean() - n[emp].mean()); D.append(nu * nd + dc())
        for _ in range(nsub):
            if which == "levelA":
                s.step(t, dt)
            else:
                s.step(dt)
            t += dt
    return times, np.array(I), np.array(D)


if __name__ == "__main__":
    U = float(sys.argv[1]) if len(sys.argv) > 1 else 5.0
    tmax = float(sys.argv[2]) if len(sys.argv) > 2 else 50.0
    dt = float(sys.argv[3]) if len(sys.argv) > 3 else 0.002
    out = {}
    for mode in ("factorized", "sic"):
        t0 = time.time(); ta, Ia, Da = run(U, mode, tmax, dt, "levelA"); tA = time.time() - t0
        t0 = time.time(); tu, Iu, Du = run(U, mode, tmax, dt, "unified"); tU = time.time() - t0
        line = (f"U/J = {U:g}, {mode:10s}: max |I_unified - I_levelA| = {np.abs(Iu - Ia).max():.1e}, "
                f"max |d_unified - d_levelA| = {np.abs(Du - Da).max():.1e}  [Level A {tA:.0f}s, unified {tU:.0f}s]")
        ref = os.path.join(ROOT, "benchmarks", "hubbard", "data", f"fast_cdw_U{U:g}_{mode}.npz")
        if os.path.exists(ref) and abs(dt - 0.002) < 1e-12:
            R = np.load(ref)
            k = min(len(R["I"]), len(Iu))
            line += f"; vs stored {os.path.basename(ref)}: {np.abs(Iu[:k] - R['I'][:k]).max():.1e}"
        print(line, flush=True)
        out[f"{mode}_t"] = tu; out[f"{mode}_I_unified"] = Iu; out[f"{mode}_I_levelA"] = Ia
        out[f"{mode}_d_unified"] = Du; out[f"{mode}_d_levelA"] = Da
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    np.savez(os.path.join(HERE, "data", "levelA_limit.npz"), U=U, **out)
