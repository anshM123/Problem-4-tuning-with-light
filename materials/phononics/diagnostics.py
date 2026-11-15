"""Controlled diagnostics: separating the method from the MACE potential.

Part 1 (polarization reversal, Mankowsky et al. 2017).  The same Maxwell + coherent-mode +
incoherent-bath equations (drive.py: FDTD light, DFPT mode charges, MACE three-phonon decay
and bath, unchanged) are rerun with LITERATURE energy surfaces for the TO4-TO1 (Q_IR-Q_P)
system instead of MACE:
  L1  Abalmasov, arXiv:1910.05020 (2019), thermodynamic potential with only the parameters he
      takes from first principles / spectroscopy: Landau double well for Q_P (nu_P = 7.5 THz,
      Q_P0 = 3.1 Å amu^1/2 from ab initio displacements), biquadratic coupling
      c2 Q_P^2 Q_IR^2, c2 = 34 meV amu^-2 Å^-4 (from DFT paraelectric vs ferroelectric TO4
      frequencies); harmonic Q_IR.  No parameter adjusted to the experiment.
  L2  L1 plus the DFT anharmonicity of Q_IR of von Hoegen et al. (2018), Methods
      (a3 = 1568, a4 = 901, a5 = 7.1 meV): an all-literature first-principles surface.
  L3  Abalmasov's full model: L1 plus c1 Q_P Q_IR^3 (10 meV) and c3 Q_P^3 Q_IR (17 meV) and a
      depolarizing field E_d (1 - Q_P/Q_P0), E_d = 2.8 MV/cm.  These three were ADJUSTED TO
      EXPERIMENT B in the source.  Note: with c3 the stationary point of his potential moves to
      Q_P = 4.89, Q_IR = -0.90 Å amu^1/2 (he starts from Q_P0 = 3.1, which is not stationary), so
      L3 cannot be carried over consistently; it is reported for completeness only.
  (Subedi, PRB 92, 214303 (2015) states that LiNbO3 has a large g Q_P Q_IR^2 coupling but gives
  no LiNbO3 coefficients; the Mankowsky 2017 supplement uses a dimensionless model; Hardy et
  al., arXiv:2312.04732, give no coefficients and only relative fields.)
  The TO4 harmonic frequency is kept at the MACE value in every surface (Omega_IR^2 =
  w4^2 - 2 c2 Q_P0^2), so the linear optics is identical to the default run.

Part 2 (frequency renormalization, von Hoegen et al. 2018, Ext. Data Fig. 6).  The field scan at
17.5 THz is repeated with the electric field recorded to 18 um, and omega' is extracted like for
like from both measured observables: the polarization-rotation signal (E integrated along the
800 nm light line, group index 2.30, as in the paper's FDTD analysis) and the SH signal
(dQ_IR/dt along the 400 nm light line over 2 um, group index 3.03), FFT over the full measured
window (-0.3 to 1.2 ps) and over the free-oscillation window (0.25 to 1.2 ps), for the
polariton (12-17 THz), the drive and the TO (17.6-21 THz) components; the fit f = sqrt(1+aE^2)
of the paper is applied to the model points.

Part 3 (method comparison): the minimal single-cell model of the literature, in which the
incident vacuum field drives the modes directly (no reflection, no screening), with each surface.

Part 4: experiment B with the pump on the MACE TO4 resonance (same relative detuning as 19 THz
from the measured TO4), the largest internal drive the Maxwell treatment allows.

Writes diag_part1..4.pkl.  Usage: python diagnostics.py [1] [2] [3] [4]  (all by default;
about 20, 12, 1 and 3 min).
"""
import os, sys, time, pickle
import numpy as np
import drive as d
from params import EXP_A, fluence_to_field, C0

HERE = os.path.dirname(os.path.abspath(__file__))
eVconv = 1.602176634e-19 / (1e-20 * 1.66053906660e-27)       # eV/(Å^2 amu) -> s^-2


def w2_eV(nu_THz):
    return (2 * np.pi * nu_THz * 1e12) ** 2 / eVconv


W4_2 = w2_eV(d.FA1[3])
NU_P, QP0, C2 = 7.5, 3.1, 0.034


def make_surface(kind):
    wP2 = w2_eV(NU_P)
    cP = wP2 / (2 * QP0 ** 2)
    Om2 = W4_2 - 2 * C2 * QP0 ** 2
    c1, c3 = (0.010, 0.017) if kind == "L3" else (0.0, 0.0)
    a3, a4, a5 = (1.56765, 0.9008, 0.0071) if kind == "L2" else (0.0, 0.0, 0.0)

    def V(P, I):
        return (-wP2 / 4 * P ** 2 + cP / 4 * P ** 4 + 0.5 * Om2 * I ** 2 + a3 / 3 * I ** 3 + a4 / 4 * I ** 4
                + a5 / 5 * I ** 5 + C2 * P ** 2 * I ** 2 + c1 * P * I ** 3 + c3 * P ** 3 * I)
    from scipy.optimize import minimize
    r = minimize(lambda x: V(x[0], x[1]), [QP0, 0.0], method="Nelder-Mead",
                 options=dict(xatol=1e-12, fatol=1e-14, maxiter=4000))
    P0, I0 = r.x
    V0 = V(P0, I0)
    return (lambda q1, q4: V(P0 + q1, I0 + q4) - V0), P0, I0


def part1():
    out = {}
    sgnP = float(d.an["sgnP"])
    for kind in ("L1", "L2", "L3"):
        Vlit, P0, I0 = make_surface(kind)
        d.set_surface("literature", Vlit, Ed_MVcm=2.8 if kind == "L3" else 0.0, Q0=P0)
        print(f"{kind}: equilibrium Q_P = {P0:.3f}, Q_IR = {I0:+.3f} Å amu^1/2", flush=True)
        rows = []
        for nu, Fl in [(19.0, [6.7, 12.8, 20.7, 28.7, 37.8, 51.8, 65.9, 94.5]), (21.0, [24.4, 52.4, 90.2]),
                       (22.0, [49.4, 103.7, 173.8])]:
            for F in Fl:
                r = d.run(fluence_to_field(F), nu, quiet=True)
                ratio, tt = d.polarization_peak(r, sgnP=sgnP)
                q1 = r["Qs"][:, 0]
                rows.append(dict(nu=nu, F=F, Pmin=ratio.min(), q1min=q1.min(),
                                 tmin=tt[np.argmin(ratio)], ratio=ratio if F in (12.8, 51.8, 94.5) else None,
                                 t=tt if F in (12.8, 51.8, 94.5) else None))
                print(f"  {kind} nu={nu} F={F:6.1f}: P/P0 min {ratio.min():+.3f}  q1 min {q1.min():+.2f}", flush=True)
        out[kind] = rows
    d.set_surface("mace")
    return out


def pr_signal(r, ng=2.30):
    t, z, E = r["t"], r["zE"], r["E"]
    S = np.zeros_like(t)
    for k, zz in enumerate(z):
        S += np.interp(t + zz * ng / C0, t, E[:, k], left=0, right=0)
    return S * (z[1] - z[0])


def peaks(t, s):
    res = {}
    for wname, (a, b) in (("full", (-0.3e-12, 1.2e-12)), ("free", (0.25e-12, 1.2e-12))):
        f, F = d.spectrum(t, s, a, b, pad=32)
        res[wname] = dict(dom=d.peak(f, F, 10, 24)[0], pol=d.peak(f, F, 12, 17.0)[0],
                          to=d.peak(f, F, 17.6, 21)[0])
    return res


def part2():
    out = {}
    for label, nu, pot in (("mace_17.5", EXP_A["nu"], "mace"), ("mace_16.8", 16.8, "mace"),
                           ("dftIR_17.5", EXP_A["nu"], "dft")):
        d.set_surface("mace")
        if pot == "dft":
            d.set_ir_potential("dft")
        rows = []
        for E0 in (2.0, 5.7, 8.0, 11.6, 14.8, 17.6, 20.9):
            r = d.run(E0, nu, quiet=True, record_E=True)
            spr = pr_signal(r); ssh = d.sh_signal(r)
            rows.append(dict(E0=E0, PR=peaks(r["t"], spr), SH=peaks(r["t"], ssh),
                             tr=dict(t=r["t"], PR=spr, SH=ssh) if E0 in (2.0, 20.9) else None))
            p = rows[-1]
            print(f"  {label} E0={E0:5.1f}: PR full dom/pol/TO {p['PR']['full']['dom']:.3f}/{p['PR']['full']['pol']:.3f}/"
                  f"{p['PR']['full']['to']:.3f}  free dom {p['PR']['free']['dom']:.3f} | SH full dom/pol/TO "
                  f"{p['SH']['full']['dom']:.3f}/{p['SH']['full']['pol']:.3f}/{p['SH']['full']['to']:.3f} "
                  f"free dom {p['SH']['free']['dom']:.3f}", flush=True)
        out[label] = rows
    d.set_surface("mace")
    return out


def local_run(E0_MVcm, nu, fwhm_fs=150.0, t0=-0.5e-12, t1=1.4e-12, dt=0.05e-15):
    """Single-cell model (no Maxwell): the incident vacuum field acts directly on the modes, as in
    the minimal models of the literature.  Same surface, mode charges and 300 K linewidths."""
    from params import QSI
    zs = d.ZSI
    gam = 2 * np.pi * 1e12 * np.array(list(d.gam_lo) + [float(d.ph["gam_to4"][0])])
    Q = np.zeros((4, 1)); V = np.zeros((4, 1))
    n = int((t1 - t0) / dt); ts = t0 + dt * np.arange(n)
    qs = np.zeros((n, 4)); om = 2 * np.pi * nu * 1e12; E0 = E0_MVcm * 1e8
    for i, t in enumerate(ts):
        E = E0 * np.sin(om * t) * np.exp(-2 * np.log(2) * t ** 2 / (fwhm_fs * 1e-15) ** 2)
        if d.EDEP["Ed"] != 0.0:
            E = E - d.EDEP["Ed"] * Q[0, 0] / d.EDEP["Q0"]
        F = d.forces(Q)[:, 0] * d.FSI + zs * E
        VSI = V[:, 0] * QSI
        VSI = ((1 - gam * dt) * VSI + dt * F) / (1 + gam * dt)
        V[:, 0] = VSI / QSI
        Q[:, 0] += dt * V[:, 0]
        qs[i] = Q[:, 0]
    dP = (zs[None, :] * qs * QSI).sum(1) / d.OMEGA
    s = 35e-15 / (2 * np.sqrt(2 * np.log(2)))
    k = np.exp(-0.5 * ((np.arange(-3000, 3001) * dt) / s) ** 2); k /= k.sum()
    ratio = 1 + float(d.an["sgnP"]) * np.convolve(dP, k, mode="same") / d.PS0
    return ratio, qs, ts


def part3():
    """Unscreened-field single-cell model with each surface (method comparison)."""
    out = {}
    for kind in ("MACE", "L1", "L2"):
        if kind == "MACE":
            d.set_surface("mace")
        else:
            Vlit, P0, I0 = make_surface(kind)
            d.set_surface("literature", Vlit, Ed_MVcm=2.8 if kind == "L3" else 0.0, Q0=P0)
        rows = []
        for F in (6.7, 12.8, 20.7, 28.7, 37.8, 51.8, 65.9, 94.5):
            ratio, qs, ts = local_run(fluence_to_field(F), 19.0)
            m = slice(1000, -1000)
            rows.append(dict(F=F, Pmin=ratio[m].min(), q4max=np.abs(qs[:, 3]).max(), q1min=qs[:, 0].min()))
            print(f"  local {kind} F={F:6.1f}: P/P0 min {ratio[m].min():+.3f}  |Q_IR| max {np.abs(qs[:,3]).max():.2f}  "
                  f"q1 min {qs[:,0].min():+.2f}", flush=True)
        out[kind] = rows
    # frequency of the driven mode in the same single-cell model, 17.5 THz pump (experiment A)
    for pot in ("mace", "dft"):
        d.set_surface("mace")
        if pot == "dft":
            d.set_ir_potential("dft")
        rows = []
        for E0 in (2.0, 5.7, 11.6, 17.6, 20.9):
            ratio, qs, ts = local_run(E0, EXP_A["nu"])
            v = np.gradient(qs[:, 3], ts)
            pk = [d.peak(*d.spectrum(ts, v, a, b, pad=8), 10, 24)[0]
                  for a, b in ((-0.3e-12, 1.2e-12), (0.25e-12, 1.2e-12))]
            rows.append(dict(E0=E0, qmax=np.abs(qs[:, 3]).max(), full=pk[0], free=pk[1]))
            print(f"  local freq {pot} E0={E0:5.1f}: |Q_IR| max {rows[-1]['qmax']:.2f}  peak full {pk[0]:.3f} free {pk[1]:.3f}", flush=True)
        out[f"freq_{pot}"] = rows
    d.set_surface("mace")
    return out


def part4():
    """Experiment B with the pump at the same relative detuning from the MACE TO4 as 19 THz is
    from the measured TO4 (i.e. exactly on the MACE resonance, 18.24 THz): upper bound on the
    internal drive within the Maxwell treatment."""
    out = {}
    sgnP = float(d.an["sgnP"])
    nu = 19.0 * d.FA1[3] / 19.0
    for kind in ("MACE", "L1"):
        if kind == "MACE":
            d.set_surface("mace")
        else:
            Vlit, P0, I0 = make_surface(kind)
            d.set_surface("literature", Vlit, Q0=P0)
        rows = []
        for F in (12.8, 51.8, 94.5):
            r = d.run(fluence_to_field(F), nu, quiet=True)
            ratio, tt = d.polarization_peak(r, sgnP=sgnP)
            rows.append(dict(F=F, Pmin=ratio.min(), q4=np.abs(r["Q4s"]).max()))
            print(f"  resonant {kind} nu={nu:.2f} F={F:6.1f}: P/P0 min {ratio.min():+.3f} |Q_IR| {np.abs(r['Q4s']).max():.2f}", flush=True)
        out[kind] = rows
    d.set_surface("mace")
    return out


if __name__ == "__main__":
    parts = sys.argv[1:] or ["1", "2", "3", "4"]
    t0 = time.time()
    if "1" in parts:
        pickle.dump(part1(), open(os.path.join(HERE, "diag_part1.pkl"), "wb"))
    if "4" in parts:
        pickle.dump(part4(), open(os.path.join(HERE, "diag_part4.pkl"), "wb"))
    if "3" in parts:
        pickle.dump(part3(), open(os.path.join(HERE, "diag_part3.pkl"), "wb"))
    if "2" in parts:
        pickle.dump(part2(), open(os.path.join(HERE, "diag_part2.pkl"), "wb"))
    print(f"done ({time.time()-t0:.0f}s)")
