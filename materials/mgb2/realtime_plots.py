"""Figure and summary table for the real-time MgB2 runs (figures/mgb2_realtime.png, rt_out/rt_summary.txt).
Reads rt_out/*.npz written by realtime_runs.py; missing runs are skipped."""
import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import realtime as rt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "rt_out")
FIG = os.path.join(HERE, "..", "..", "figures")
THZ = rt.THZ

# measured (refs/experimental_targets.md, Paper C)
S3B = {21: 0.2, 47: 1.0, 60: 1.4, 68: 2.1, 79: 3.8, 87: 5.4, 116: 6.6, 150: 6.6}          # max dE_probe (a.u.)
S3C = {"150": [(0.5, 1.6), (0.8, 2.55), (1.1, 1.9), (1.4, 1.5), (1.6, 1.0), (1.85, 0.8), (2.0, 1.0)],
       "79": [(0.5, 0.7), (0.75, 1.2), (1.2, 0.9), (1.7, 0.75), (2.0, 1.1)],
       "68": [(0.5, 0.25), (0.75, 0.65), (1.15, 0.4), (1.85, 0.65)],
       "nb21": [(0.5, 0.0), (0.78, 0.25), (1.0, 0.0), (1.9, 0.6)]}
S3D = {"150": [(0.5, 2.6), (0.65, 1.4), (1.0, 1.55), (1.8, 1.0), (2.0, 0.5)],
       "79": [(0.5, 3.7), (0.75, 3.0), (1.0, 2.5), (1.8, 1.0)],
       "68": [(0.5, 4.7)], "nb21": [(0.5, 5.6)]}
F3B = {12: (5e-7, 2.8e-7), 15: (2e-6, 1.2e-6), 18: (1e-5, 3.5e-6), 21: (3.5e-5, 9e-6), 25: (2.1e-4, 1.6e-5), 30: (4.5e-4, 2.2e-5)}
OLD = {"s1_08": 2.23e4, "s2_05": 1.98e4, "s2_10": 1.15e4}      # bottleneck.py (Owen-Scalapino state)


def load(name):
    p = os.path.join(OUT, name)
    return np.load(p) if os.path.exists(p) else None


def at(f, s, x):
    return s[np.argmin(abs(f - x))]


def harm_I(d):
    t, ENL = d["t"], d["ENL"]
    dt = t[1] - t[0]
    F = np.fft.rfft(ENL * np.hanning(len(ENL)), 8 * len(ENL))
    fr = np.fft.rfftfreq(8 * len(ENL), dt)
    P = np.abs(F) ** 2
    return P[(fr > 0.3) & (fr < 1.0)].sum(), P[(fr > 1.6) & (fr < 2.2)].sum()


def main():
    lines = []
    tags = [t for t in ("nb21", "10", "20", "30", "47", "58", "68", "79", "90", "100", "116", "130", "150") if load(f"rt_pump_{t}.npz") is not None]
    runs = {t: load(f"rt_pump_{t}.npz") for t in tags}
    val = load("rt_validate.npz")
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    # (a) gaps vs time
    a = ax[0, 0]
    for t, d in runs.items():
        tt = d["t"] - (d["tpr"] - rt_tau(t))
        l, = a.plot(tt, 2 * abs(d["D"][:, 0]) / THZ, label=f"{t} kV/cm" if t != "nb21" else "21 kV/cm (narrowband)")
        a.plot(tt, 2 * abs(d["D"][:, 1]) / THZ, ":", color=l.get_color())
        kt = d["kin_t"] - (d["tpr"] - rt_tau(t))
        a.plot(kt, 2 * abs(d["kin_D"][:, 0]) / THZ, color=l.get_color(), alpha=0.5)
        a.plot(kt, 2 * abs(d["kin_D"][:, 1]) / THZ, ":", color=l.get_color(), alpha=0.5)
    a.axvline(2.0, color="k", lw=0.5)
    a.set_xscale("symlog", linthresh=5); a.set_xlim(-2, 100)
    a.set_xlabel("time after pump peak (ps)"); a.set_ylabel(r"$2\Delta_\pi$ (solid), $2\Delta_\sigma$ (dotted)  (THz)")
    a.set_title("real-time gaps (TDBdG + e-ph kinetics)"); a.legend(fontsize=6)
    # (b) sigma at 2 ps, 150 kV/cm
    a = ax[0, 1]
    if "150" in runs:
        d = runs["150"]
        m = (d["f"] > 0.35) & (d["f"] < 1.3)
        a.plot(d["f"][m], d["sig"][m].real / 100, "r-", label=r"$\sigma_1$ real-time probe")
        a.plot(d["f"][m], d["sig"][m].imag / 100, "b-", label=r"$\sigma_2$ real-time probe")
        a.plot(d["fmb"], d["sigmb"].real / 100, "r--", label=r"$\sigma_1$ MB with computed $f$, $\Delta$")
        a.plot(d["fmb"], d["sigmb"].imag / 100, "b--", label=r"$\sigma_2$ MB with computed $f$, $\Delta$")
    a.plot(*np.array(S3C["150"]).T * [[1], [1e4]], "ro", label=r"measured $\sigma_1$")
    a.plot(*np.array(S3D["150"]).T * [[1], [1e4]], "bs", label=r"measured $\sigma_2$")
    a.plot([0.8, 0.5, 1.0], [OLD["s1_08"], OLD["s2_05"], OLD["s2_10"]], "k*", ms=9, label="previous bottleneck state")
    a.set_xlabel("f (THz)"); a.set_ylabel(r"$\Omega^{-1}$cm$^{-1}$"); a.set_ylim(-0.5e4, 6e4)
    a.set_title("2 ps after 150 kV/cm, 7 K"); a.legend(fontsize=6)
    # (c) field dependence of sigma at 2 ps
    a = ax[0, 2]
    Es, s108, s205 = [], [], []
    for t, d in runs.items():
        E = 21.0 if t == "nb21" else float(t)
        Es.append(E); s108.append(at(d["f"], d["sig"], 0.8).real / 100); s205.append(at(d["f"], d["sig"], 0.5).imag / 100)
    o = np.argsort(Es)
    a.plot(np.array(Es)[o], np.array(s108)[o], "r-o", label=r"model $\sigma_1(0.8)$")
    a.plot(np.array(Es)[o], np.array(s205)[o], "b-o", label=r"model $\sigma_2(0.5)$")
    me = [(21, 0.25, 5.6), (68, 0.65, 4.7), (79, 1.2, 3.7), (150, 2.55, 2.6)]
    a.plot([x[0] for x in me], [x[1] * 1e4 for x in me], "rs", mfc="none", label=r"measured $\sigma_1(0.8)$")
    a.plot([x[0] for x in me], [x[2] * 1e4 for x in me], "bs", mfc="none", label=r"measured $\sigma_2(0.5)$")
    _, bp = beam_profile_report(E0s=(30, 47, 60, 68, 79, 87, 100, 116, 130, 150))
    if bp is not None:
        fb, rb = bp
        Eb = sorted(rb)
        a.plot(Eb, [at(fb, rb[e][0], 0.8).real / 100 for e in Eb], "r--", label=r"model $\sigma_1(0.8)$, beam-profile average")
        a.plot(Eb, [at(fb, rb[e][0], 0.5).imag / 100 for e in Eb], "b--", label=r"model $\sigma_2(0.5)$, beam-profile average")
    a.set_xlabel("pump peak field (kV/cm)"); a.set_ylabel(r"$\Omega^{-1}$cm$^{-1}$ at 2 ps"); a.legend(fontsize=6)
    a.set_title("field dependence of the pumped state (S3c,d)")
    # (d) dE_probe vs field
    a = ax[1, 0]
    dEm = None
    if val is not None:
        # rms change of the transmitted probe waveform over the probe window (same timing as validate run)
        dEm = []
        for t in [x for x in np.array(list(runs))[o]]:
            d = runs[t]
            tv = val["tp"] - 3.6
            dd = np.interp(tv + d["tpr"], d["tp"], d["dE"])
            dEm.append(np.sqrt(np.mean((dd - val["dE"]) ** 2)))
        Eo = np.array(Es)[o]
        dEm = np.array(dEm)
        a.plot(Eo, dEm / dEm[-1] * S3B[150], "k-o", label=r"model rms $\delta E_{probe}$ (scaled to 150 kV/cm)")
        lines.append("pump-probe signal (rms change of transmitted probe, normalized to 150 kV/cm): "
                     + ", ".join(f"{e:g}: {x / dEm[-1]:.2f}" for e, x in zip(Eo, dEm)))
        lines.append("measured max dE_probe / value at 150 kV/cm: " + ", ".join(f"{k}: {v / 6.6:.2f}" for k, v in S3B.items()))
    if bp is not None:
        a.plot(Eb, [rb[e][1] / rb[150][1] * S3B[150] for e in Eb], "k--", label="model, beam-profile average")
    a.plot(list(S3B), list(S3B.values()), "rs", label="measured max $\\delta E_{probe}$ (Fig. S3b)")
    a.set_xlabel("pump peak field (kV/cm)"); a.set_ylabel("a.u."); a.legend(fontsize=6)
    a.set_title("pump-probe signal vs pump field")
    # (e) FH / TH vs field
    a = ax[1, 1]
    import re
    hs = [p for p in glob.glob(os.path.join(OUT, "rt_harm_*.npz")) if re.fullmatch(r"rt_harm_[0-9.]+\.npz", os.path.basename(p))]
    hs = sorted(hs, key=lambda p: float(os.path.basename(p)[8:-4]))
    if hs:
        Eh = np.array([float(p.split("_")[-1][:-4]) for p in hs])
        I = np.array([harm_I(np.load(p)) for p in hs])
        k12 = np.argmin(abs(Eh - 12))
        a.loglog(Eh, I[:, 0] / I[k12, 0] * F3B[12][0], "b-o", label="FH model (scaled at 12 kV/cm)")
        a.loglog(Eh, I[:, 1] / I[k12, 1] * F3B[12][1], "r-o", label="TH model (scaled at 12 kV/cm)")
        lines.append("narrowband harmonics (I normalized to 12 kV/cm):")
        for e, (x, y) in zip(Eh, I):
            lines.append(f"  E = {e:4.0f} kV/cm: FH {x / I[k12, 0]:8.2f}  TH {y / I[k12, 1]:8.2f}   E^6 {(e/12)**6:8.2f}"
                         + (f"   measured FH {F3B[int(e)][0]/F3B[12][0]:7.1f} TH {F3B[int(e)][1]/F3B[12][1]:7.1f}" if int(e) in F3B else ""))
    Em = np.array(list(F3B)); a.loglog(Em, [F3B[e][0] for e in Em], "bs", mfc="none", label="FH measured (Fig. 3b)")
    a.loglog(Em, [F3B[e][1] for e in Em], "rs", mfc="none", label="TH measured")
    a.loglog(Em, 5e-7 * (Em / 12) ** 6, "k:", label=r"$\propto E^6$")
    a.set_xlabel("peak field per pulse (kV/cm)"); a.set_ylabel(r"$I_{NL}$"); a.legend(fontsize=6)
    a.set_title("narrowband 0.63 THz, 7 K")
    # (f) distributions and phonons at 2 ps (150 kV/cm)
    a = ax[1, 2]
    for t, c in (("150", "r"), ("79", "g"), ("47", "b")):
        if t not in runs:
            continue
        d = runs[t]
        o2 = np.argsort(d["snap_E0"])
        a.semilogy(d["snap_E0"][o2], np.clip(d["snap_f0"][o2], 1e-6, 1), c + "-", lw=0.8, label=rf"$f_\pi(E)$, {t} kV/cm")
        o3 = np.argsort(d["snap_E1"])
        a.semilogy(d["snap_E1"][o3], np.clip(d["snap_f1"][o3], 1e-6, 1), c + ":", lw=0.8, label=rf"$f_\sigma(E)$")
        a.semilogy(d["phW"], np.clip(d["snap_n"], 1e-6, 10), c + "--", lw=0.6, label="phonons n(Ω)")
    a.semilogy(np.linspace(0, 20, 50), rt.fermi(np.linspace(0, 20, 50), 7.0), "k:", lw=0.5, label="7 K")
    a.set_xlim(0, 20); a.set_ylim(1e-5, 2); a.set_xlabel("E, Ω (meV)"); a.legend(fontsize=5)
    a.set_title("distributions 2 ps after the pump")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "mgb2_realtime.png"), dpi=130)
    # summary
    lines.insert(0, "pumped state at tau = 2 ps (Ohm^-1 cm^-1); probe = real-time weak probe, MB = quasi-stationary Mattis-Bardeen")
    for t, d in runs.items():
        sp = [at(d["f"], d["sig"], x) / 100 for x in (0.5, 0.8, 1.0)]
        sm = [at(d["fmb"], d["sigmb"], x) / 100 for x in (0.5, 0.8, 1.0)]
        lines.append(f"  {t:>5}: U_abs = {float(d['Uabs_end'])*17.46e-6:6.3f} J/mol, 2Dpi(2ps) = {2*abs(d['snap_D'][0])/THZ:.3f}, "
                     f"2Dsig(2ps) = {2*abs(d['snap_D'][1])/THZ:.3f} THz | probe s1(0.8) {sp[1].real:.2e} s2(0.5) {sp[0].imag:.2e} "
                     f"s2(1.0) {sp[2].imag:.2e} | MB s1(0.8) {sm[1].real:.2e} s2(0.5) {sm[0].imag:.2e} s2(1.0) {sm[2].imag:.2e}")
    if val is not None:
        lines.append("equilibrium 7 K, real-time probe vs Mattis-Bardeen (twoband.sigma_total):")
        import twoband as tb
        D = tb.gaps(7.0, rt.WC)
        for x in (0.3, 0.5, 0.8, 1.0, 1.5):
            s = at(val["f"], val["sig"], x) / 100; r = tb.sigma_total(x, D, 7.0) / 100
            lines.append(f"  {x:.1f} THz: s1 {s.real:.2e} ({r.real:.2e})  s2 {s.imag:.2e} ({r.imag:.2e})")
    open(os.path.join(OUT, "rt_summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def rt_tau(tag):
    return 2.0


def bottleneck_report(tag="47", times=(7.0, 20.0, 50.0, 98.0)):
    """Rothwarf-Taylor diagnostic from the stored continuation: ratio of phonon pair breaking to
    recombination for each band (1 = quasiparticles and pair-breaking phonons in mutual equilibrium),
    and the fraction of the excess energy held by phonons."""
    d = load(f"rt_pump_{tag}.npz")
    if d is None:
        return []
    import realtime_runs as rr
    bands = [rt.Band(0, rr.NB["N"][0], rr.NB["W"][0], rr.NB["seed"]), rt.Band(1, rr.NB["N"][1], rr.NB["W"][1], rr.NB["seed"] + 1)]
    ph = rt.Phonons(7.0)
    out = [f"Rothwarf-Taylor diagnostic, {tag} kV/cm (P/R = pair breaking / recombination by phonons):"]
    t0 = d["tpr"] - 2.0
    for tt in times:
        k = np.argmin(abs(d["kin_t"] - t0 - tt))
        D = d["kin_D"][k]; fs = [d["kin_f0"][k], d["kin_f1"][k]]; ph.n = d["kin_n"][k]
        PR = []
        for A in bands:
            a = A.a
            P = R = 0.0
            for B in bands:
                b = B.a
                Ea = np.sqrt(A.eps ** 2 + abs(D[a]) ** 2); Eb = np.sqrt(B.eps ** 2 + abs(D[b]) ** 2)
                Om = Ea[:, None] + Eb[None, :]
                g = rt._a2f(a, b, Om) * 0.5 * (1 - (A.eps / Ea)[:, None] * (B.eps / Eb)[None, :] + (abs(D[a]) / Ea)[:, None] * (abs(D[b]) / Eb)[None, :])
                n = ph.at(Om); fa = fs[a][:, None]; fb = fs[b][None, :]
                P += np.sum(g * (1 - fa) * (1 - fb) * n); R += np.sum(g * fa * fb * (n + 1))
            PR.append(P / R)
        Ue = sum(2 * rt.NUD[B.a] * B.de * np.sum(np.sqrt(B.eps ** 2 + abs(D[B.a]) ** 2) * fs[B.a]) for B in bands) * rt.MEV
        Up = ph.energy()
        out.append(f"  t = {tt:5.1f} ps: 2Dpi = {2*abs(D[0])/THZ:.3f}, 2Dsig = {2*abs(D[1])/THZ:.3f} THz; P/R (pi, sigma) = {PR[0]:.2f}, {PR[1]:.2f}; "
                   f"phonon excess energy {Up:.2e} J/m^3")
    return out


# ------------------------------------------------------------------ pump beam profile
def beam_profile_report(E0s=(47, 60, 68, 79, 87, 116, 150)):
    """Gaussian pump and probe beams of equal size at the sample.  The detected probe field is the
    probe-mode-weighted average of the local transmission; the local pump peak field is E0 u with
    u = exp(-r^2/w^2) and weight 2u du on u in [0, 1].  Local responses are interpolated between the
    computed pump fields (0 = equilibrium run).  Returns lines with sigma_eff(2 ps) and the rms
    probe change versus the nominal peak field E0."""
    val = load("rt_validate.npz")
    tags = [t for t in ("10", "20", "30", "47", "58", "68", "79", "90", "100", "116", "130", "150") if load(f"rt_pump_{t}.npz") is not None]
    if val is None or not tags:
        return [], None
    f = val["f"]
    Es = np.array([0.0] + [float(t) for t in tags])
    sig = [val["sig"]] + [np.interp(f, load(f"rt_pump_{t}.npz")["f"], load(f"rt_pump_{t}.npz")["sig"]) for t in tags]
    tr = np.array([2 / (1 + rt.NS + rt.Z0 * s * rt.DFILM) for s in sig])          # local transmission
    tv = val["tp"] - 3.6
    dEs = [val["dE"]] + [np.interp(tv + load(f"rt_pump_{t}.npz")["tpr"], load(f"rt_pump_{t}.npz")["tp"], load(f"rt_pump_{t}.npz")["dE"]) for t in tags]
    dEs = np.array(dEs)
    u = np.linspace(0, 1, 401)[1:]; w = 2 * u; w /= w.sum()
    out = ["pump beam profile (equal Gaussian pump and probe spots), state at 2 ps vs nominal peak field:"]
    res = {}
    for E0 in E0s:
        El = E0 * u
        ta = sum(wi * np.array([np.interp(e, Es, tr[:, k]) for k in range(len(f))]) for wi, e in zip(w[::8] * 8, El[::8]))
        sg = (2 / ta - 1 - rt.NS) / (rt.Z0 * rt.DFILM)
        dd = sum(wi * np.array([np.interp(e, Es, dEs[:, k]) for k in range(dEs.shape[1])]) for wi, e in zip(w[::8] * 8, El[::8]))
        rms = np.sqrt(np.mean((dd - val["dE"]) ** 2))
        res[E0] = (sg, rms)
        out.append(f"  E0 = {E0:4.0f} kV/cm: s1(0.8) {at(f, sg, 0.8).real/100:.2e}  s2(0.5) {at(f, sg, 0.5).imag/100:.2e}  "
                   f"s2(1.0) {at(f, sg, 1.0).imag/100:.2e}")
    r150 = res[max(E0s)][1]
    out.append("  rms probe change normalized at 150 kV/cm: " + ", ".join(f"{e}: {res[e][1]/r150:.2f}" for e in E0s))
    return out, (f, res)


def harm_beam_report():
    """Narrowband harmonics with the beam profile: the detected nonlinear field is the
    detection-mode-weighted average of the local E_NL, E_NL,det(t; E0) = int_0^1 E_NL(t; E0 u) du
    (Gaussian beams of equal size; local fields between computed points interpolated, E^3 below
    the lowest computed field)."""
    import re
    hs = [p for p in glob.glob(os.path.join(OUT, "rt_harm_*.npz")) if re.fullmatch(r"rt_harm_[0-9.]+\.npz", os.path.basename(p))]
    if len(hs) < 3:
        return []
    Eh = np.array(sorted(float(os.path.basename(p)[8:-4]) for p in hs))
    tr = {e: np.load(os.path.join(OUT, f"rt_harm_{e:g}.npz")) for e in Eh}
    t = tr[Eh[0]]["t"]
    X = np.array([tr[e]["ENL"] for e in Eh])

    def local(E):
        if E <= Eh[0]:
            return X[0] * (E / Eh[0]) ** 3
        k = np.searchsorted(Eh, E) - 1
        k = min(k, len(Eh) - 2)
        w = (E - Eh[k]) / (Eh[k + 1] - Eh[k])
        return (1 - w) * X[k] + w * X[k + 1]
    u = (np.arange(200) + 0.5) / 200
    res = {}
    for E0 in [e for e in (12, 15, 18, 21, 25, 30) if e <= Eh[-1]]:
        EN = sum(local(E0 * x) for x in u) / len(u)
        res[E0] = harm_I({"t": t, "ENL": EN})
    out = ["narrowband harmonics with the beam profile (I normalized to 12 kV/cm):"]
    for e, (a, b) in res.items():
        out.append(f"  E0 = {e:4.0f}: FH {a / res[12][0]:7.2f}  TH {b / res[12][1]:7.2f}  E^6 {(e / 12) ** 6:7.2f}"
                   + (f"  measured FH {F3B[e][0] / F3B[12][0]:6.1f} TH {F3B[e][1] / F3B[12][1]:6.1f}" if e in F3B else ""))
    return out


if __name__ == "__main__":
    main()
    extra, _ = beam_profile_report()
    extra += harm_beam_report()
    for tg in ("47", "68", "150"):
        extra += bottleneck_report(tg)
    print("\n".join(extra))
    with open(os.path.join(OUT, "rt_summary.txt"), "a") as fh:
        fh.write("\n".join(extra) + "\n")
