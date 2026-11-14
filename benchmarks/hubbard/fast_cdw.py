"""Fast-method CDW melting on the 11-site chain (same protocol as ed_cdw.py).
Usage: python fast_cdw.py U [mode] [dt]
mode: factorized (default) | g1g2 | hf
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lightmat.fluct_hubbard import PPLadder

L = 11
occ = [1, 3, 5, 7, 9]
emp = [i for i in range(L) if i not in occ]
U = float(sys.argv[1]); mode = sys.argv[2] if len(sys.argv) > 2 else "factorized"
dt = float(sys.argv[3]) if len(sys.argv) > 3 else 0.002
tmax, tout = 50.0, 0.1
h = np.zeros((L, L))
for i in range(L - 1):
    h[i, i + 1] = h[i + 1, i] = -1.0
r0 = np.diag([1.0 if i in occ else 0.0 for i in range(L)])
t0 = time.time()
s = PPLadder(h, U if mode != "hf" else U, r0, r0, mode="factorized" if mode == "hf" else mode)
if mode == "hf":   # time-dependent Hartree(-Fock): drop correlations
    s.x = s.x[:, :0]; s.y = s.y[:, :0]
    s._gsel_orig = s._gsel
    s._gsel = lambda x, y, ru, rd: {k: 0.0 * np.eye(L) if k != "diag" else np.zeros(L)
                                    for k in ("iiji", "ijjj", "iiij", "jijj", "diag")}
    s._pauli_P = lambda v, ru, rd: 0.0 * v
times = np.arange(0, tmax + 1e-9, tout)
n_out, d_out, dc_out = [], [], []
t = 0.0
nsub = int(round(tout / dt))
for T in times:
    nu, nd = s.densities()
    dc = s.dcorr()
    n_out.append(nu + nd); dc_out.append(dc); d_out.append(nu * nd + dc)
    if not np.all(np.isfinite(nu)):
        print("diverged at", T); break
    for _ in range(nsub):
        s.step(t, dt); t += dt
n = np.array(n_out); d = np.array(d_out); dc = np.array(dc_out)
I = n[:, occ].mean(1) - n[:, emp].mean(1)
out = os.path.join(os.path.dirname(__file__), "data", f"fast_cdw_U{U:g}_{mode}.npz")
np.savez(out, t=times[:len(I)], n=n, d=d, dcorr=dc, I=I, runtime=time.time() - t0)
print(f"U={U} mode={mode}: {time.time()-t0:.1f}s  I(50)={I[-1]:.3f}  sum dcorr(50)={dc[-1].sum():.3f}")
