"""Window check: pump-only dynamics to tau = 2 ps with wider level windows (W, N scaled at fixed
level spacing) to test the finite-window error at the largest vector potentials.
usage: python window_check.py E_kV_cm W_pi W_sigma"""
import sys
import time
import numpy as np
import realtime as rt

E = float(sys.argv[1]); Wp = float(sys.argv[2]); Ws = float(sys.argv[3])
N = (int(round(300 * Wp / 15.0)), int(round(250 * Ws / 18.0)))
t0 = time.time()
s = rt.Sim(T0=7.0, N=N, W=(Wp, Ws), seed=1)
Ef = rt.broadband(E * 1e5, 1.6)
Amax = [0.0]
ts, Efs, Ds, As = rt.run(s, Ef, 3.6, 0.01)
D = s.gaps(s.R)
print(f"E = {E:g} kV/cm, W = ({Wp:g}, {Ws:g}) meV, N = {N}: absorbed {s.Uabs*17.46e-6:.3f} J/mol, "
      f"max |A| = {np.max(abs(As)):.2e} V s/m, 2Dpi(2 ps) = {2*abs(D[0])/rt.THZ:.3f}, 2Dsig(2 ps) = {2*abs(D[1])/rt.THZ:.3f} THz "
      f"({time.time()-t0:.0f} s)", flush=True)
