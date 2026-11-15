"""Photoexcited graphite end to end with the unified solver (lightmat/unified.py).

light -> coherent Bloch electrons -> e-e and e-ph collisions -> phonon populations -> three-phonon decay.

Ingredients (all from materials/graphite, unchanged):
  * pi bands: nearest-neighbour tight binding of eph_tb.py (hbar v_F = 6.4 eV A), full BZ k-grid;
  * light: 35 fs, 800 nm, 12 mJ/cm^2 pulse; the in-plane field inside the 70 nm film from the same
    thin-film optics (n = 2.72 + 1.47i) as dynamics.py, depth-averaged |E|^2; Peierls coupling
    h(k + eA/hbar), no assumption about the absorbed energy;
  * e-ph: the first-principles energy-resolved kernel eph_kernel.npz (GW-scaled as in dynamics.py);
  * ph-ph: the MACE three-phonon process lists of dynamics.py (24 x 24 grid, 6 branches);
  * e-e, two options (argument EE):
      'rta15' (default): conserving relaxation toward the Fermi-Dirac distribution with the same
        energy and density, at the measured 15 fs thermalization time of photoexcited graphite
        (Breusing et al. 2011; Brida et al. 2013; the value used by dynamics_neq.py), applied per k
        to populations and coherences (dephasing 1/tau);
      'U<x>': the solver's own Markov limit of the second-order G2 equation with a local
        interaction x * U*, U* = 1.6 |t| = 4.8 eV (optimal local U of graphene, Schueler et al.,
        PRL 111, 036601 (2013)).  It gives a 5 ps out-scattering time at 0.775 eV and 3000 K:
        a local U misses the long-range Coulomb scattering that thermalizes graphite's carriers.
Stages (one code, lightmat/unified.py):
  1. t = -100..+150 fs: coherent Bloch density matrices rho_k (BlochCarriers) on a 600 x 600 grid,
     exact 2x2 propagation in the field, Markov e-e and e-ph collisions (populations and dephasing)
     from energy-binned occupations, phonon populations from the same e-ph kernel;
  2. 0.15..1 ps: random-phase (isotropic) electrons F(e) on the energy grid, e-e + e-ph + ph-ph;
  3. 1..100 ps: thermalized electrons (fixed point of the e-e operator, checked), e-ph + ph-ph.
Output: unified_out.npz (same layout as materials/graphite/dynamics_out.npz) and run.log.

Usage: python run.py [Nk=450] [EE=rta15 | U<factor>] [field_scale=1.0] [tag]
"""
import os, sys, time
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
GDIR = os.path.join(ROOT, "materials", "graphite")
sys.path.insert(0, ROOT); sys.path.insert(0, GDIR)
from lightmat.unified import markov_ee_local, EphKernel, markov_phph, BlochCarriers
from eph_tb import bands, B, KB_eV, THz2eV, T_HOP

NK = int(sys.argv[1]) if len(sys.argv) > 1 else 450
EE = sys.argv[2] if len(sys.argv) > 2 else "rta15"
UF = float(EE[1:]) if EE.startswith("U") else 0.0
TAU_EE = float(EE[3:]) if EE.startswith("rta") else None
FS = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
TAG = sys.argv[4] if len(sys.argv) > 4 else ""
argv = sys.argv; sys.argv = sys.argv[:1]
import dynamics as D                 # graphite Level B tables: procs, freq, GW scale, T0 (no run)
sys.argv = argv
HB = 0.6582119569                    # eV fs
T0 = D.T0
U = 1.6 * abs(T_HOP) * UF            # eV, optimal local interaction (Schueler et al. 2013)

# ---------------------------------------------------------------- energy grid, DOS
DEc = 0.02
eg = np.arange(-9.2, 9.2 + 1e-9, DEc)
nb = len(eg)
f_ = (np.arange(1200) + 0.5) / 1200
kd = np.stack(np.meshgrid(f_, f_, indexing="ij"), -1).reshape(-1, 2) @ B
ed, _ = bands(kd)
ibd = np.clip(np.rint((ed - eg[0]) / DEc).astype(int), 0, nb - 1)
dos_cell = np.bincount(ibd.ravel(), minlength=nb) / (len(kd) * DEc)      # per cell, per spin
rho_site = dos_cell / 2
gstates = 2 * dos_cell * DEc                                             # states per bin per cell, both spins
del kd, ed, ibd

# ---------------------------------------------------------------- e-ph kernel on the coarse grid
K = np.load(os.path.join(GDIR, "eph_kernel.npz"))
Wf = K["W"] * (D.scale[:, None, None] if D.GW else 1.0)                 # 1/ps, (Nq, 6, 3000)
Nq, nbr = D.freq.shape
ic = np.clip(np.rint((K["centers"] - eg[0]) / DEc).astype(int), 0, nb - 1)
Wc = np.zeros((Nq * nbr, nb))
Wflat = Wf.reshape(Nq * nbr, -1)
for j in range(Wflat.shape[1]):
    Wc[:, ic[j]] += Wflat[:, j]
valid = D.valid.reshape(-1)
Wc[~valid] = 0.0
hw = D.hw.reshape(-1)
ephk = EphKernel(Wc / 1000.0, np.where(valid, hw, 0.0), DEc, gstates, Nq)   # rates in 1/fs
procs = D.procs                                                          # 1/ps weights


def fermi(T, mu=0.0):
    return 1 / (np.exp(np.clip((eg - mu) / (KB_eV * T), -700, 700)) + 1)


_Tt = np.geomspace(1.0, 80000, 4000)
_Et = None


def ee_rates(F):
    """e-e in/out rates per state (1/fs)."""
    global _Et
    if TAU_EE is None:
        gin, gout = markov_ee_local(F, rho_site, U, DEc)
        return gin / HB, gout / HB
    if _Et is None:
        _Et = np.array([np.sum(gstates * eg * fermi(T)) for T in _Tt])
    T = float(np.interp(np.sum(gstates * eg * F), _Et, _Tt))
    Feq = fermi(T)                       # same energy; density conserved by particle-hole symmetry (mu = 0)
    return Feq / TAU_EE, (1 - Feq) / TAU_EE


def e_of_F(F):
    return float(np.sum(gstates * eg * F))


# ---------------------------------------------------------------- light
LAM, TAU, F_INC = 800e-9, 35.0, 12.0
w_L = 2 * np.pi * 2.99792458e8 / LAM * 1e-15                           # rad/fs
hw_L = HB * w_L                                                          # eV
n_c, d_f = 2.72 + 1.47j, 70e-9
r01 = (1 - n_c) / (1 + n_c); r12 = (n_c - 1) / (n_c + 1); t01 = 2 / (1 + n_c)
beta = 2 * np.pi * n_c * d_f / LAM
z = np.linspace(0, 1, 401)
Ez = t01 * (np.exp(1j * beta * z) + r12 * np.exp(1j * beta * (2 - z))) / (1 + r01 * r12 * np.exp(2j * beta))
E2avg = float(np.mean(np.abs(Ez) ** 2))
eps2 = float(np.imag(n_c ** 2))
EPS0, C0 = 8.8541878128e-12, 2.99792458e8
A_check = EPS0 * (w_L * 1e15) * eps2 * E2avg * d_f / (EPS0 * C0)       # absorbed fraction (linear)
E0_inc = np.sqrt(2 * F_INC * 10 / (EPS0 * C0 * TAU * 1e-15 * np.sqrt(np.pi / (4 * np.log(2)))))  # V/m
E0 = E0_inc * np.sqrt(E2avg) * 1e-10 * FS                               # V/A, depth-rms in-film field
V_cell = 5.258 * 3.35                                                    # A^3 per 2-atom cell
E0_in_SI = E0_inc * np.sqrt(E2avg) * FS                                # V/m
E_lin = (0.5 * EPS0 * eps2 * (w_L * 1e15) * E0_in_SI ** 2 * TAU * 1e-15 * np.sqrt(np.pi / (4 * np.log(2)))
         * V_cell * 1e-30 / 1.602176634e-19)                             # eV/cell absorbed, linear optics
Efield = lambda t: E0 * np.exp(-2 * np.log(2) * t ** 2 / TAU ** 2) * np.cos(w_L * t)

log = open(os.path.join(HERE, f"run{TAG}.log"), "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); log.write(s + "\n"); log.flush()
say(f"e-e: {'relaxation time %g fs' % TAU_EE if TAU_EE else 'local Markov, U = %.2f eV' % U}; Nk = {NK}^2; "
    f"field scale {FS}; hbar w_L = {hw_L:.3f} eV")
say(f"thin film: depth-averaged |E/E_inc|^2 = {E2avg:.4f}; linear absorptance from it {A_check:.3f} "
    f"(dynamics.py: {D.Aabs:.3f}); in-film peak field {E0:.4f} V/A; linear-optics absorbed energy "
    f"{E_lin*1e3:.1f} meV/cell (dynamics.py: {D.E_abs_cell*1e3:.1f})")

# ---------------------------------------------------------------- stage 1: coherent Bloch electrons
t_start, t_coh, dt = -100.0, 150.0, 0.05
nsub_c = 10                                   # collision update every 0.5 fs
f_ = (np.arange(NK) + 0.5) / NK
kc = np.stack(np.meshgrid(f_, f_, indexing="ij"), -1).reshape(-1, 2) @ B
car = BlochCarriers(kc, bands, KB_eV * T0, eg)
N = np.where(valid, 1 / np.expm1(np.maximum(hw / (KB_eV * T0), 1e-12)), 0.0)
N0 = N.copy()
a_vec = np.zeros(2)
E_e0 = car.energy(a_vec)
t = t_start
steps = int(round((t_coh - t_start) / dt))
tlog, Ee_log, Eph_log, coh_log, Fsnap = [], [], [], [], {}
t0 = time.time()
for it in range(steps):
    tm = t + 0.5 * dt
    a_mid = a_vec + np.array([-Efield(t + 0.25 * dt) * 0.5 * dt / HB, 0.0])
    car.step_coherent(a_mid, dt / HB)
    a_vec = a_vec + np.array([-Efield(tm) * dt / HB, 0.0])
    t += dt
    if (it + 1) % nsub_c == 0:
        s, cnt = car.binned(a_vec)
        F = np.where(cnt > 0, s / np.maximum(cnt, 1), fermi(T0))
        gin_ee, gout_ee = ee_rates(F)
        dN, gin_ep, gout_ep = ephk.rates(F, N)
        Gin = gin_ee + gin_ep; Gout = gout_ee + gout_ep
        car.step_collisions(a_vec, Gin, Gout, nsub_c * dt)
        if (it + 1) % (20 * nsub_c) == 0:
            N = N + 20 * nsub_c * dt * (dN + markov_phph(N, procs) / 1000.0)
            Ee = car.energy(a_vec) - E_e0
            Eph = float(np.sum(hw * (N - N0)) / Nq)
            tlog.append(t); Ee_log.append(Ee); Eph_log.append(Eph); coh_log.append(car.coherence(a_vec))
            if any(abs(t - x) < 0.5 * nsub_c * dt * 20 for x in (0.0, 25.0, 50.0, 100.0)):
                Fsnap[round(t)] = F.copy()
            say(f"t = {t:7.1f} fs  E_e = {Ee*1e3:7.2f}  E_ph = {Eph*1e3:6.2f} meV/cell  "
                f"|rho_cv| = {coh_log[-1]:.2e}  [{time.time()-t0:.0f}s]")
s, cnt = car.binned(np.zeros(2))
F = np.where(cnt > 0, s / np.maximum(cnt, 1), fermi(T0))
E_abs = (car.energy(np.zeros(2)) - E_e0) + float(np.sum(hw * (N - N0)) / Nq)
say(f"absorbed energy (electrons + phonons at {t_coh:.0f} fs): {E_abs*1e3:.1f} meV/cell; "
    f"linear optics {E_lin*1e3:.1f}; Level B {D.E_abs_cell*1e3:.1f}")
if os.environ.get("STAGE1_ONLY"):
    say(f"E_abs / field_scale^2 = {E_abs / FS ** 2 * 1e3:.1f} meV/cell; E_lin / field_scale^2 = {E_lin / FS ** 2 * 1e3:.1f}")
    sys.exit(0)
coh_stage = dict(t_coh=np.array(tlog), Ee_coh=np.array(Ee_log), Eph_coh=np.array(Eph_log), rho_cv=np.array(coh_log))
del car, kc

# ---------------------------------------------------------------- stage 2: isotropic electrons F(e)
E_tot = E_abs
t_iso_end = 1000.0 if TAU_EE else 3000.0
out_t, out_n, out_Te, out_dev = [], [], [], []
Tt = np.geomspace(T0, 60000, 2000)
Et = np.array([e_of_F(fermi(T)) for T in Tt]) - e_of_F(fermi(T0))
Te_of = lambda E: float(np.interp(E, Et, Tt))
E_F0 = e_of_F(fermi(T0))
nextout = t_coh
pp = markov_phph(N, procs) / 1000.0
k = 0
while t < t_iso_end - 1e-9:
    dt2 = 0.5 if t < 400 else 2.0
    if TAU_EE is None and t > 1000:
        dt2 = 5.0
    gin_ee, gout_ee = ee_rates(F)
    dN, gin_ep, gout_ep = ephk.rates(F, N)
    Gin = gin_ee + gin_ep; Gout = gout_ee + gout_ep
    tot = Gin + Gout
    feq = np.where(tot > 0, Gin / np.where(tot > 0, tot, 1), F)
    F = feq + (F - feq) * np.exp(-tot * dt2)
    if k % 20 == 0:
        pp = markov_phph(N, procs) / 1000.0
    N = N + dt2 * (dN + pp)
    t += dt2; k += 1
    if t >= nextout - 1e-9:
        Ee = e_of_F(F) - E_F0
        Te = Te_of(Ee)
        dev = float(np.max(np.abs(F - fermi(Te))))
        out_t.append(t / 1000); out_n.append(N.reshape(Nq, nbr).astype(np.float32)); out_Te.append(Te); out_dev.append(dev)
        Eph = float(np.sum(hw * (N - N0)) / Nq)
        if any(abs(t - x) < 1 for x in (200, 300, 500, 1000, 2000, 3000)):
            Fsnap[round(t)] = F.copy()
            say(f"t = {t:6.0f} fs  Te = {Te:6.0f} K  max|F - FD(Te)| = {dev:.2e}  E_e = {Ee*1e3:6.2f}  "
                f"E_ph = {Eph*1e3:6.2f}  total {1e3*(Ee+Eph):6.2f} meV/cell  [{time.time()-t0:.0f}s]")
        nextout += 10.0 if t < 1000 else 50.0

# ---------------------------------------------------------------- stage 3: thermal electrons
Ee_sw = e_of_F(F) - E_F0
say(f"switch to thermal electrons at {t/1000:.2f} ps: max |F - FD(Te)| = {out_dev[-1]:.2e}")


def rhs(tt, y):
    Ee, n = y[0], y[1:]
    Te = max(Te_of(Ee), 0.5 * T0)
    dN, _, _ = ephk.rates(fermi(Te), n)
    r_ep = dN * 1000.0                                  # 1/ps
    r_pp = markov_phph(n, procs)
    resid = np.sum(hw * r_pp)                           # same energy-residual projection as dynamics.py
    cw = np.where(valid, hw * n * (n + 1), 0.0)
    r_pp = r_pp - resid * cw / np.sum(hw * cw)
    return np.concatenate([[-np.sum(hw * r_ep) / Nq], r_ep + r_pp])


tev = np.concatenate([np.linspace(1.05, 2.0, 20), np.linspace(2.05, 10, 160), np.linspace(10.5, 100, 180)])
tev = tev[tev > t / 1000 + 1e-6]
sol = solve_ivp(rhs, (t / 1000, 100.0), np.concatenate([[Ee_sw], N]), t_eval=tev, method="RK45", rtol=1e-4, atol=1e-7)
say(f"thermal stage: {sol.message} nfev = {sol.nfev} [{time.time()-t0:.0f}s]")
tt = np.concatenate([out_t, sol.t])
nn = np.concatenate([np.array(out_n), sol.y[1:].T.reshape(-1, Nq, nbr).astype(np.float32)])
Te_all = np.concatenate([out_Te, [Te_of(e) for e in sol.y[0]]])
Eph_end = float(np.sum(hw * (sol.y[1:, -1] - N0)) / Nq); Ee_end = sol.y[0, -1]
say(f"energy: absorbed {E_abs*1e3:.2f} meV/cell; at 100 ps electrons {Ee_end*1e3:.2f} + phonons {Eph_end*1e3:.2f} "
    f"= {1e3*(Ee_end+Eph_end):.2f}")
# prepend the initial state at t = -0.1 ps for the analysis (n(t0) = thermal)
tt = np.concatenate([[-0.1], tt]); nn = np.concatenate([N0.reshape(1, Nq, nbr).astype(np.float32), nn])
Te_all = np.concatenate([[T0], Te_all])
np.savez(os.path.join(HERE, f"unified_out{TAG}.npz"), t=tt, n=nn, Te=Te_all, freq=D.freq, qfrac=D.E["qfrac"],
         E_abs_cell=E_abs, E_lin=E_lin, U=U, eg=eg, F_snap_t=np.array(sorted(Fsnap)),
         F_snap=np.array([Fsnap[k] for k in sorted(Fsnap)]), dos_cell=dos_cell, **coh_stage)
say("written", f"unified_out{TAG}.npz")
