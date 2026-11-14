"""Graphite with all three interactions explicit in one calculation: non-thermal electrons
(electron-electron scattering), electron-phonon scattering with the energy-resolved kernel,
and three-phonon scattering.

Electrons: energy-resolved occupation f(e) on the pi-band density of states g(e) (both spins,
per cell), e in [-3, 3] eV, bins of 5 meV.  Terms in df/dt:
  pump   : 800 nm (1.55 eV) interband injection e -> e + 1.55 eV around -0.775 eV, Pauli blocked,
           35 fs Gaussian envelope, normalised to the absorbed energy of dynamics.py;
  e-e    : energy-conserving relaxation towards the Fermi-Dirac distribution with the same
           energy (chemical potential fixed at the charge-neutral point; impact ionisation /
           Auger processes change the carrier number), rate
           Gamma_ee(e) = c_ee max(|e|, pi k T*) / hbar, c_ee = 0.06, which gives the measured
           ~15 fs thermalization of 0.8 eV carriers in graphene/graphite (Breusing et al.
           PRB 83, 153410 (2011); Brida et al. Nat. Commun. 4, 1987 (2013));
  e-ph   : sum_q,nu emission/absorption between bins i and i + s_q,nu (s = hw/de) with the kernel
           W of eph_kernel.py (GW scaling as in dynamics.py), exactly the generalisation of
           the thermal rates.
Phonons: as in dynamics.py (e-ph from the same kernel, full three-phonon operator).
After 2 ps the electrons are thermal (checked: |f - f_FD| < 1e-3) and the run continues with
the thermal model of dynamics.py, which is then exact.
"""
import os, sys, time
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
sys.argv = sys.argv[:1]
import dynamics as D                                 # phonon data, three-phonon operator, GW scale
from eph_tb import bands, B, KB_eV, THz2eV

HERE = os.path.dirname(os.path.abspath(__file__))
K = np.load(os.path.join(HERE, "eph_kernel.npz"))
W, cen, DE = K["W"].copy(), K["centers"], float(K["DE"])
NB = len(cen)
if D.GW:
    W *= D.scale[:, None, None]
Nq, nb = D.freq.shape
# density of states per bin (states per cell, both spins)
f_ = (np.arange(600) + 0.5) / 600
kk = np.stack(np.meshgrid(f_, f_, indexing="ij"), -1).reshape(-1, 2) @ B
ek, _ = bands(kk)
gbin = 2 * np.histogram(ek.ravel(), bins=NB, range=(cen[0] - DE / 2, cen[-1] + DE / 2))[0] / len(kk)
gbin = np.maximum(gbin, 1e-12)
# flattened e-ph transition lists
shift = np.rint(D.hw / DE).astype(int)
QI, JI, II, WW = [], [], [], []
for q in range(Nq):
    for j in range(nb):
        s = shift[q, j]
        if s <= 0 or not D.valid[q, j]:
            continue
        ii = np.nonzero(W[q, j, :NB - s] > 0)[0]
        QI.append(np.full(len(ii), q * nb + j)); II.append(ii); WW.append(W[q, j, ii])
QI = np.concatenate(QI); II = np.concatenate(II); WW = np.concatenate(WW)
SS = shift.reshape(-1)[QI]
kB = KB_eV; T0 = D.T0
HW_PUMP = 1.55
C_EE = 0.06
HBAR = 6.582119569e-16 * 1e12            # eV ps


def fermi(T):
    return 1 / (np.exp(np.clip(cen / (kB * T), -700, 700)) + 1)


f0 = fermi(T0)
E0 = np.sum(gbin * cen * f0) * DE
# pump profile: injection from e_v = -0.775 eV (Gaussian 50 meV) to e_v + 1.55 eV
sp = int(round(HW_PUMP / DE))
prof = np.exp(-0.5 * ((cen + HW_PUMP / 2) / 0.05) ** 2); prof[NB - sp:] = 0
prof /= np.sum(prof * gbin * DE)                 # per unit transitions per cell
def pump_rate(t):                                # transitions / (cell ps)
    return D.power(t) / HW_PUMP


def T_star(f):
    E = np.sum(gbin * cen * f) * DE
    return brentq(lambda T: np.sum(gbin * cen * fermi(T)) * DE - E, 1.0, 2e5)


def rhs(t, y):
    f = y[:NB]; n = y[NB:].reshape(Nq, nb)
    nf = n.reshape(-1)
    df = np.zeros(NB)
    # pump (Pauli blocked)
    pr = pump_rate(t)
    if pr > 1e-12:
        tr = pr * prof * gbin * DE * f * (1 - np.roll(f, -sp))
        tr[NB - sp:] = 0
        df -= tr / (gbin * DE)
        df += np.roll(tr, sp) / (gbin * DE)
    # e-e relaxation (energy conserving)
    Ts = T_star(f)
    Gam = C_EE * np.maximum(np.abs(cen), np.pi * kB * Ts) / HBAR
    fFD = fermi(Ts)
    # weight rates so that energy is conserved exactly: adjust T until sum g e Gam (f - fFD) = 0
    def econs(T):
        return np.sum(gbin * cen * Gam * (f - fermi(T)))
    try:
        Tc = brentq(econs, 1.0, 2e5)
    except ValueError:
        Tc = Ts
    df -= Gam * (f - fermi(Tc))
    # e-ph
    fi = f[II]; fj = f[II + SS]; nn = nf[QI]
    R = WW * (fj * (1 - fi) * (nn + 1) - fi * (1 - fj) * nn)        # emission, electron j -> i
    dn_ep = np.bincount(QI, weights=R, minlength=Nq * nb)
    dN = np.bincount(II, weights=R, minlength=NB) - np.bincount(II + SS, weights=R, minlength=NB)
    df += dN / Nq / (gbin * DE)
    # phonon-phonon
    r_pp = D.phph_rate(n)
    resid = np.sum(D.hw * r_pp)
    cw = np.where(D.valid, D.hw * n * (n + 1), 0.0)
    r_pp = r_pp - resid * cw / np.sum(D.hw * cw)
    return np.concatenate([df, dn_ep + r_pp.reshape(-1)])


if __name__ == "__main__":
    t0 = time.time()
    y = np.concatenate([f0, D.nB(T0).reshape(-1)])
    tout = np.concatenate([np.linspace(-0.1, 2.0, 106)])
    dt = 0.002; t = -0.1
    rec_t, rec_y = [], []
    k_out = 0
    while t < 2.0 - 1e-9:
        if k_out < len(tout) and t >= tout[k_out] - 1e-9:
            rec_t.append(t); rec_y.append(y.copy()); k_out += 1
        k1 = rhs(t, y); k2 = rhs(t + dt / 2, y + dt / 2 * k1)
        k3 = rhs(t + dt / 2, y + dt / 2 * k2); k4 = rhs(t + dt, y + dt * k3)
        y = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt
    rec_t.append(t); rec_y.append(y.copy())
    f_end = y[:NB]; Tend = T_star(f_end)
    print(f"non-thermal stage done {time.time()-t0:.0f}s; T* at 2 ps = {Tend:.0f} K; "
          f"max |f - f_FD| = {np.abs(f_end - fermi(Tend)).max():.2e}", flush=True)
    # thermal continuation
    y2 = np.concatenate([[np.interp(Tend, D.Etab_T, D.Etab_E) - D.E0], y[NB:]])
    t2 = np.concatenate([np.linspace(2.05, 10, 160), np.linspace(10.5, 100, 180)])
    sol = solve_ivp(D.rhs, (2.0, 100.0), y2, t_eval=t2, rtol=1e-4, atol=1e-7)
    ts = np.concatenate([rec_t, sol.t])
    n_all = np.concatenate([np.array(rec_y)[:, NB:], sol.y[1:].T]).reshape(len(ts), Nq, nb)
    Te = np.concatenate([[T_star(yy[:NB]) for yy in rec_y], D.Te_of(sol.y[0])])
    fs = np.array(rec_y)[:, :NB]
    np.savez(os.path.join(HERE, "dynamics_neq.npz"), t=ts, n=n_all.astype(np.float32), Te=Te, f=fs.astype(np.float32),
             f_t=np.array(rec_t), centers=cen, freq=D.freq, qfrac=D.E["qfrac"])
    Eph = (n_all * D.hw[None]).sum((1, 2)) / Nq
    print(f"energy: absorbed {D.E_abs_cell:.4f} eV/cell; phonon excess at 100 ps {Eph[-1]-Eph[0]:.4f}; total time {time.time()-t0:.0f}s")
