"""Direct cost comparison: two-time Kadanoff-Baym equations (second Born, lightmat/kbe.py)
versus the time-linear positive-factorized pp-ladder (lightmat/fluct_hubbard.py), same
system (11-site Hubbard chain, CDW initial state, U/J = 5), same time step dt = 0.02,
same machine, single process.  Wall times are measured for increasing numbers of time
steps N_t and the KBE cost is fitted to a N_t^3 law and extrapolated to the benchmark
propagation time t = 50 hbar/J (N_t = 2500)."""
import os, sys, time, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lightmat.kbe import KBE2B
from lightmat.fluct_hubbard import PPLadder

L, U, dt = 11, 5.0, 0.02
occ = [1, 3, 5, 7, 9]
h = np.zeros((L, L))
for i in range(L - 1):
    h[i, i + 1] = h[i + 1, i] = -1.0
r0 = np.diag([1.0 if i in occ else 0.0 for i in range(L)])
res = {"kbe": [], "fast": []}
for N in (50, 100, 200, 400):
    k = KBE2B(h, U, r0, dt, N); t0 = time.time()
    for _ in range(N):
        k.step()
    tk = time.time() - t0
    s = PPLadder(h, U, r0, r0); t0 = time.time(); t = 0.0
    for _ in range(N):
        s.step(t, dt); t += dt
    tf = time.time() - t0
    res["kbe"].append((N, tk)); res["fast"].append((N, tf))
    print(f"N_t={N:4d}: KBE {tk:8.2f} s   time-linear {tf:6.2f} s   ratio {tk/tf:7.1f}", flush=True)
Nk = np.array(res["kbe"]); Nf = np.array(res["fast"])
pk = np.polyfit(np.log(Nk[-2:, 0]), np.log(Nk[-2:, 1]), 1)
pf = np.polyfit(np.log(Nf[-2:, 0]), np.log(Nf[-2:, 1]), 1)
Nt = int(50 / dt)
tk50 = Nk[-1, 1] * (Nt / Nk[-1, 0]) ** 3
tf50 = Nf[-1, 1] * (Nt / Nf[-1, 0])
mem = 4 * (Nt + 1) ** 2 * L * L * 16 / 1e9
print(f"measured scaling exponents: KBE {pk[0]:.2f}, time-linear {pf[0]:.2f}")
print(f"t = 50 hbar/J (N_t = {Nt}): KBE ~ {tk50/3600:.0f} h (N_t^3 extrapolation), memory {mem:.0f} GB; "
      f"time-linear {tf50:.0f} s; speed-up ~ {tk50/tf50:.0f}x")
json.dump({"kbe": res["kbe"], "fast": res["fast"], "kbe_t50_s": tk50, "fast_t50_s": tf50, "kbe_mem_GB": mem,
           "exp_kbe": pk[0], "exp_fast": pf[0]}, open(os.path.join(os.path.dirname(__file__), "timing.json"), "w"), indent=1)
