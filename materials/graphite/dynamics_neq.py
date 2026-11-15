"""Graphite with all three interactions explicit in one calculation.

Electrons are split into the photoexcited non-thermal population and the thermal bath.
  * Non-thermal carriers: N_nt (electrons + holes per cell) injected by the 800 nm pump at
    e = +/- 0.775 eV (35 fs envelope, width DELTA = 0.1 eV; absorbed energy as in dynamics.py).
    - Electron-electron scattering moves them into the thermal bath at the measured ~15 fs
      thermalization rate of photoexcited carriers in graphite/graphene (Breusing et al.
      PRB 83, 153410 (2011); Brida et al. Nat. Commun. 4, 1987 (2013)), depositing their
      remaining energy there (exact energy conservation).
    - Electron-phonon scattering: each carrier emits phonons (q, nu) at the golden-rule rate
      given by the energy-resolved kernel of eph_kernel.py at its energy (final states empty),
          dn_q/dt = N_nt (n_q + 1) S_q / (g_inj DELTA),  S_q = sum_{j in window} W_q(j - s_q),
      g_inj = pi-band density of states at 0.775 eV; the carriers lose hw per emission.
  * Thermal bath: temperature Te(E_e); emission/absorption A_q(Te) exactly as in dynamics.py.
Phonons: full three-phonon operator (dynamics.py).
Energy: absorbed = E_nt + E_e(excess) + E_ph(excess) at all times (checked at the end).
"""
import os, sys, time
import numpy as np
from scipy.integrate import solve_ivp
sys.argv = sys.argv[:1]
import dynamics as D
from eph_tb import bands, B

HERE = os.path.dirname(os.path.abspath(__file__))
K = np.load(os.path.join(HERE, "eph_kernel.npz"))
W, cen, DE = K["W"].copy(), K["centers"], float(K["DE"])
if D.GW:
    W *= D.scale[:, None, None]
Nq, nb = D.freq.shape
G_EE = 1 / 0.015                     # 1/ps
E_INJ, DELTA = 0.775, 0.1
f_ = (np.arange(900) + 0.5) / 900
kk = np.stack(np.meshgrid(f_, f_, indexing="ij"), -1).reshape(-1, 2) @ B
ek, _ = bands(kk)
g_inj = 2 * np.mean(np.abs(ek[:, 1] - E_INJ) < DELTA / 2) / DELTA        # states / (eV cell), both spins
win = np.nonzero(np.abs(cen - E_INJ) < DELTA / 2)[0]
shift = np.rint(D.hw / DE).astype(int)
S = np.zeros((Nq, nb))
for q in range(Nq):
    for j in range(nb):
        if D.valid[q, j] and shift[q, j] > 0:
            S[q, j] = W[q, j, win - shift[q, j]].sum()
Wem = S / (g_inj * DELTA)                                               # 1/ps per carrier, per mode
rate0 = np.sum(Wem * (D.nB(D.T0) + 1)) / Nq
print(f"phonon-emission rate of a 0.775 eV carrier: {rate0:.1f} /ps (time {1e3/rate0:.0f} fs)", flush=True)


def rhs(t, y):
    Ee, Nnt, Ent = y[0], y[1], y[2]
    n = y[3:].reshape(Nq, nb)
    Te = max(D.Te_of(Ee), D.T0 * 0.5)
    P = D.power(t)
    gen = P / (2 * E_INJ)                       # carriers per cell per ps
    em = max(Nnt, 0.0) * Wem * (n + 1)          # dn/dt per mode from non-thermal carriers
    Eloss_em = np.sum(D.hw * em) / Nq
    mean_E = Ent / Nnt if Nnt > 1e-12 else E_INJ
    th = G_EE * Nnt
    r_ep = D.eph_rate(n, Te)
    r_pp = D.phph_rate(n)
    resid = np.sum(D.hw * r_pp)
    cw = np.where(D.valid, D.hw * n * (n + 1), 0.0)
    r_pp = r_pp - resid * cw / np.sum(D.hw * cw)
    dEe = th * mean_E - np.sum(D.hw * r_ep) / Nq
    return np.concatenate([[dEe, gen - th, P - th * mean_E - Eloss_em], (r_ep + r_pp + em).reshape(-1)])


if __name__ == "__main__":
    t0 = time.time()
    y0 = np.concatenate([[0.0, 0.0, 0.0], D.nB(D.T0).reshape(-1)])
    tout = np.concatenate([np.linspace(-0.1, 2.0, 211), np.linspace(2.05, 10, 160), np.linspace(10.5, 100, 180)])
    sol = solve_ivp(rhs, (-0.1, 100.0), y0, t_eval=tout, method="RK45", rtol=1e-5, atol=1e-8)
    print("integrated", sol.status, f"{time.time()-t0:.0f}s", flush=True)
    Ee, Nnt, Ent = sol.y[0], sol.y[1], sol.y[2]
    n = sol.y[3:].T.reshape(len(sol.t), Nq, nb)
    Eph = (n * D.hw[None]).sum((1, 2)) / Nq
    np.savez(os.path.join(HERE, "dynamics_neq.npz"), t=sol.t, n=n.astype(np.float32), Te=D.Te_of(Ee), Nnt=Nnt,
             Ent=Ent, freq=D.freq, qfrac=D.E["qfrac"])
    i05 = np.argmin(abs(sol.t - 0.3))
    print("energy check (eV/cell): absorbed %.4f, final e + nt + ph excess %.4f" % (D.E_abs_cell, Ee[-1] + Ent[-1] + Eph[-1] - Eph[0]))
