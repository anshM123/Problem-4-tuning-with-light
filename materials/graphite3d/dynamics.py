"""Coupled electron / mode-resolved phonon dynamics of photoexcited AB graphite (4 atoms, 12 branches;
same equations as ../graphite/dynamics.py, with two pi-electron layers per cell).

Electrons: thermalized (Te) after absorption (e-e scattering ~ tens of fs), energy
  dE_e/dt = P(t) - (1/Nq) sum_{q,j} hw_qj (dn_qj/dt)|_e-ph
Phonons (every branch j and grid wavevector q):
  dn/dt = A_qj(Te) [ (n+1) e^{-hw/kTe} - n ]                     (electron-phonon)
        + C sum_t w_t sum_{j'j''} |Phi|^2 { [(n+1)n'n'' - n(n'+1)(n''+1)] g(w-w'-w'')
             + [(n+1)(n'+1)n'' - n n'(n''+1)] g(w+w'-w'')
             + [(n+1)n'(n''+1) - n(n'+1)n''] g(w-w'+w'') }       (three-phonon)
A_qj(Te) from eph_rates.npz (tight-binding pi bands, DFT-calibrated coupling),
|Phi|^2 and C from phph.npz (MACE cubic force constants, prefactor calibrated against
phono3py).  Energy is conserved exactly by both collision terms (up to delta broadening).

Usage: python dynamics.py [absorbed_fluence_mJcm2]
"""
import os, sys, time
import numpy as np
from scipy.integrate import solve_ivp
import sys as _s, os as _o
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "graphite"))
from eph_tb import electron_energy, KB_eV, THz2eV
NLAYER = 2   # AB graphite: two layers per cell

HERE = os.path.dirname(os.path.abspath(__file__))
SIG_DYN = 0.8          # THz Gaussian width for energy conservation on the 24x24 grid
T0 = 295.0
FWHM = 0.035           # ps pump duration
# ---- absorbed energy -------------------------------------------------------------
# 12 mJ/cm^2 incident at 800 nm on a 70 nm free-standing graphite film; thin-film optics
# with n~ = 2.72 + 1.47i (ordinary index) gives A = 1 - R - T (computed below).
def film_absorption(nr=2.72 + 1.47j, d=70e-9, lam=800e-9):
    n0 = 1.0
    r01 = (n0 - nr) / (n0 + nr); r12 = (nr - n0) / (nr + n0)
    t01 = 2 * n0 / (n0 + nr); t12 = 2 * nr / (nr + n0)
    b = 2 * np.pi * nr * d / lam
    r = (r01 + r12 * np.exp(2j * b)) / (1 + r01 * r12 * np.exp(2j * b))
    t = t01 * t12 * np.exp(1j * b) / (1 + r01 * r12 * np.exp(2j * b))
    return 1 - abs(r) ** 2 - abs(t) ** 2


F_inc = 12.0  # mJ/cm^2
Aabs = film_absorption()
F_abs = float(sys.argv[1]) if len(sys.argv) > 1 else F_inc * Aabs * float(os.environ.get("FABS_SCALE", "1"))
V_cell = 5.258 * 6.711 * 1e-30                    # m^3 per 4-atom AB graphite cell
E_abs_cell = (F_abs * 1e-3 / 1e-4) / 70e-9 * V_cell / 1.602176634e-19   # eV per cell
print(f"absorptance {Aabs:.3f}; absorbed fluence {F_abs:.2f} mJ/cm^2; {E_abs_cell*1e3:.1f} meV/cell")

# ---- data ------------------------------------------------------------------------
E = np.load(os.path.join(HERE, "eph_rates.npz"))
P = np.load(os.path.join(HERE, "phph.npz"))
freq = P["freq"]; Nq, nb = freq.shape
assert np.allclose(freq, E["freq"], atol=1e-6)
hw = freq * THz2eV                          # eV
Tgrid = E["Te"]; Atab = E["A"].copy()      # (Nq, 6, nT) 1/ps
logT = np.log(Tgrid)
# Many-body (GW) correction of the electron-phonon coupling (Lazzeri, Attaccalite, Wirtz,
# Mauri, PRB 78, 081406 (2008)): <D^2_K> = 193 eV^2/A^2 (DFT 92), <D^2_Gamma> = 62 (DFT 45.6).
# Our tight-binding coupling reproduces DFT at K and 0.0444 vs 0.0405 eV^2 at Gamma, so the
# intervalley (near K, K') rates are scaled by 193/92 and intravalley ones by
# (62/45.6)*(0.0405/0.0444).  Set GW=0 to use DFT couplings.
GW = int(os.environ.get("GW", "1"))
if GW:
    qf = E["qfrac"][:, :2]
    qc = qf @ np.array([[1, 0], [0.5, np.sqrt(3) / 2]])   # reduced coords in hexagonal metric
    def dist(q0):
        best = np.inf * np.ones(len(qf))
        for s1 in (-1, 0, 1):
            for s2 in (-1, 0, 1):
                d = qf - np.array(q0) - np.array([s1, s2])
                cart = d @ np.array([[1, 0], [-0.5, np.sqrt(3) / 2]])
                best = np.minimum(best, np.linalg.norm(cart, axis=1))
        return best
    near_K = np.minimum(dist([1 / 3, 1 / 3]), dist([2 / 3, 2 / 3])) < dist([0, 0])
    scale = np.where(near_K, 193 / 92, (62 / 45.6) * (0.0405 / 0.0444))
    Atab *= scale[:, None, None]
C = float(P["C"]) * float(os.environ.get("PHPH_SCALE", "1"))
trip = P["trip"]; wt = P["wt"]; phi2 = P["phi2"]
valid = freq > 0.05


def gauss(x):
    return np.exp(-0.5 * (x / SIG_DYN) ** 2) / (np.sqrt(2 * np.pi) * SIG_DYN)


# sparse lists of three-phonon processes for each of the three energy conditions
print("building three-phonon process lists ...", flush=True)
t0 = time.time()
procs = []
for kind in range(3):
    Q, J, I1, J1, I2, J2, VAL = [], [], [], [], [], [], []
    for q in range(Nq):
        nt = int((wt[q] > 0).sum())
        q1 = trip[q, :nt, 0]; q2 = trip[q, :nt, 1]
        w0 = freq[q][None, :, None, None]
        w1 = freq[q1][:, None, :, None]; w2 = freq[q2][:, None, None, :]
        if kind == 0:
            dlt = w0 - w1 - w2
        elif kind == 1:
            dlt = w0 + w1 - w2
        else:
            dlt = w0 - w1 + w2
        ok = (np.abs(dlt) < 4 * SIG_DYN) & valid[q][None, :, None, None] & \
             valid[q1][:, None, :, None] & valid[q2][:, None, None, :]
        t_, j_, j1_, j2_ = np.nonzero(ok)
        val = C * wt[q, t_] * phi2[q, t_, j_, j1_, j2_] * gauss(dlt[t_, j_, j1_, j2_])
        keep = val > 0
        Q.append(np.full(keep.sum(), q)); J.append(j_[keep])
        I1.append(q1[t_[keep]]); J1.append(j1_[keep]); I2.append(q2[t_[keep]]); J2.append(j2_[keep])
        VAL.append(val[keep])
    cat = lambda L: np.concatenate(L)
    procs.append((cat(Q) * nb + cat(J), cat(I1) * nb + cat(J1), cat(I2) * nb + cat(J2), cat(VAL)))
    print(f"  kind {kind}: {len(procs[-1][3])} processes", flush=True)
print(f"done {time.time()-t0:.0f}s")

Etab_T = np.geomspace(T0, 60000, 400)
Etab_E = NLAYER * np.array([electron_energy(T, nk=300) for T in Etab_T])
E0 = Etab_E[0]


def Te_of(Ee):
    return np.interp(Ee, Etab_E - E0, Etab_T)


def power(t):
    s = FWHM / (2 * np.sqrt(2 * np.log(2)))
    return E_abs_cell * np.exp(-0.5 * (t / s) ** 2) / (np.sqrt(2 * np.pi) * s)


def phph_rate(n):
    out = np.zeros(Nq * nb)
    nf = n.reshape(-1)
    for kind, (a, b, c, v) in enumerate(procs):
        na, n1, n2 = nf[a], nf[b], nf[c]
        if kind == 0:
            r = (na + 1) * n1 * n2 - na * (n1 + 1) * (n2 + 1)
        elif kind == 1:
            r = (na + 1) * (n1 + 1) * n2 - na * n1 * (n2 + 1)
        else:
            r = (na + 1) * n1 * (n2 + 1) - na * (n1 + 1) * n2
        out += np.bincount(a, weights=v * r, minlength=Nq * nb)
    return out.reshape(Nq, nb)


def eph_rate(n, Te):
    x = np.clip((np.log(Te) - logT[0]) / (logT[1] - logT[0]), 0, len(Tgrid) - 1.000001)
    i = int(x); f = x - i
    A = (1 - f) * Atab[:, :, i] + f * Atab[:, :, i + 1]
    boltz = np.exp(-hw / (KB_eV * Te))
    r = A * ((n + 1) * boltz - n)
    return np.where(valid, r, 0.0)


def nB(T):
    return np.where(valid, 1 / np.expm1(np.maximum(hw / (KB_eV * T), 1e-12)), 0.0)


def rhs(t, y):
    Ee = y[0]; n = y[1:].reshape(Nq, nb)
    Te = max(Te_of(Ee), T0 * 0.5)
    r_ep = eph_rate(n, Te)
    r_pp = phph_rate(n)
    # exact energy conservation of the three-phonon term: the finite-width energy deltas
    # leave a small residual sum(hw * dn); remove it with a heat-capacity-weighted shift
    resid = np.sum(hw * r_pp)
    cw = np.where(valid, hw * n * (n + 1), 0.0)
    r_pp = r_pp - resid * cw / np.sum(hw * cw)
    dn = r_ep + r_pp
    dEe = power(t) - np.sum(hw * r_ep) / Nq
    return np.concatenate([[dEe], dn.reshape(-1)])


if __name__ == "__main__":
    n0 = nB(T0)
    y0 = np.concatenate([[0.0], n0.reshape(-1)])
    tout = np.concatenate([np.linspace(-0.1, 2.0, 106), np.linspace(2.05, 10, 160), np.linspace(10.5, 100, 180)])
    t0 = time.time()
    sol = solve_ivp(rhs, (-0.1, 100.0), y0, method="RK45", t_eval=tout, rtol=1e-4, atol=1e-7,
                    max_step=0.01 if False else np.inf)
    print("integrated", sol.status, sol.message, f"{time.time()-t0:.0f}s, nfev={sol.nfev}")
    Ee = sol.y[0]; n = sol.y[1:].T.reshape(len(sol.t), Nq, nb)
    Te = Te_of(Ee)
    Eph = (n * hw[None]).sum((1, 2)) / Nq
    np.savez(os.path.join(HERE, os.environ.get("OUT", "dynamics_out.npz")), t=sol.t, n=n.astype(np.float32), Te=Te,
             Ee=Ee, Eph=Eph, freq=freq, qfrac=E["qfrac"], E_abs_cell=E_abs_cell, F_abs=F_abs)
    print("energy check (eV/cell): absorbed %.4f  final electrons+phonon excess %.4f" %
          (E_abs_cell, Ee[-1] + Eph[-1] - Eph[0]))
