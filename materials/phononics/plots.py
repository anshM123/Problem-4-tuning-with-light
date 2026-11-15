"""Figure figures/phononics.png and summary.txt (model vs measured) for LiNbO3."""
import os, pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, brentq
import expdata as X
from params import fluence_to_field, ZMODE, EPS_INF, PS0

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
an = np.load(os.path.join(HERE, "anharmonic.npz"))
ph = np.load(os.path.join(HERE, "phph.npz"))
D = pickle.load(open(os.path.join(HERE, "drive_out.pkl"), "rb"))
runsA, runsB, sens, tr = D["runsA"], D["runsB"], D["sens"], D["traces"]
p8 = an["p8"]
VM = lambda q: np.polyval(p8, q)
CM = 33.35641

# ---------------------------------------------------------------- potential-shape comparison
# the published points are scaled onto the DFT curve by one factor B (Q -> B Q, U -> B^2 U);
# refit that single factor against the MACE potential
# (the published procedure: f(Q) = U(B Q) / B^2 fitted to the reconstructed points)
def cost(B, V):
    return np.sum((X.POT_U - V(X.POT_Q * B) / B ** 2) ** 2)


rB = minimize_scalar(lambda b: cost(b, VM), bounds=(0.3, 5.0), method="bounded")
B_m = rB.x
rms_m = np.sqrt(cost(B_m, VM) / len(X.POT_Q))
rms_dft = np.sqrt(cost(1.0, X.dft_potential) / len(X.POT_Q))
rms_m1 = np.sqrt(cost(1.0, VM) / len(X.POT_Q))
rB_d = minimize_scalar(lambda b: cost(b, X.dft_potential), bounds=(0.3, 5.0), method="bounded")
# excursions at 0.7 eV
def turning(V, U):
    return brentq(lambda q: V(q) - U, -3.0, -1e-3), brentq(lambda q: V(q) - U, 1e-3, 3.0)
qm_m, qp_m = turning(VM, 0.7)
qm_d, qp_d = turning(X.dft_potential, 0.7)

# ---------------------------------------------------------------- field scan (experiment A)
E = np.array([r["E0"] for r in runsA])
f1 = np.array([r["f1"] for r in runsA])
ren = f1 / f1[0]
A1 = np.array([r["A1"] for r in runsA]); A2 = np.array([r["A2"] for r in runsA]); A3 = np.array([r["A3"] for r in runsA])
m = E >= 5.0
slopes = [np.polyfit(np.log(E[m]), np.log(A[m]), 1)[0] for A in (A1, A2, A3)]
i20 = int(np.argmin(abs(E - 20.0)))
r20 = runsA[i20]
# frequency renormalization at the measured fields
ren_model = np.interp(X.REN_E, E, ren)

# ---------------------------------------------------------------- fluence scan (experiment B)
B19 = [r for r in runsB if r["nu"] == 19.0]
F19 = np.array([r["F"] for r in B19]); P19 = np.array([r["Pmin"] for r in B19])

# ---------------------------------------------------------------- TO4 linewidth
fwhm300 = 2 * ph["gam_to4"][0] * CM
fwhm_mesh = [2 * ph[f"gam_mesh{k}"][0] * CM for k in (6, 8, 10)]

# ---------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 8})
fig, ax = plt.subplots(3, 4, figsize=(15, 10.8))
a = ax[0, 0]
q = np.linspace(-3.2, 2.0, 400)
a.plot(q, VM(q), "C0", lw=2, label="MACE (frozen phonon)")
a.plot(an["qad"], an["Vad"], "C0--", lw=1, label="MACE, other modes relaxed")
a.plot(q, X.dft_potential(q), color="0.5", lw=1.5, label="DFT (von Hoegen et al.)")
w2 = an["w2fit"]
a.plot(q, 0.5 * w2 * q ** 2, ":", color="0.6", label="harmonic (MACE)")
a.plot(X.POT_Q, X.POT_U, "o", color="C3", ms=4, label="measured (scaled onto DFT)")
a.plot(X.POT_Q * B_m, X.POT_U * B_m ** 2, "s", mfc="none", color="C3", ms=3,
       label=f"measured, rescaled onto MACE (B={B_m:.2f})")
a.set_xlim(-3.2, 2.0); a.set_ylim(0, 3.8)
a.set_xlabel("Q_IR (Å amu$^{1/2}$)"); a.set_ylabel("U (eV / cell)")
a.set_title("(a) A1(TO4) potential"); a.legend(fontsize=6)

a = ax[0, 1]
for key, lab, c in (("A_20", "20 MV/cm", "C3"), ("A_2", "2 MV/cm (×10)", "C0")):
    if key in tr:
        tt = tr[key]["t"] * 1e12
        s = tr[key]["S"] * (10 if key == "A_2" else 1)
        a.plot(tt, s / np.abs(tr["A_20"]["S"]).max(), c, lw=0.8, label=lab)
a.set_xlim(-0.4, 1.2); a.set_xlabel("delay (ps)"); a.set_ylabel("SH signal ∝ ∫dQ/dt (norm.)")
a.set_title("(b) simulated SH probe signal, 17.5 THz pump"); a.legend(fontsize=7)

a = ax[0, 2]
for key, lab, c in (("A_20", "20 MV/cm", "C3"), ("A_2", "2 MV/cm", "C0")):
    if key in tr:
        F = tr[key]["F"]; f = tr[key]["f"]
        a.semilogy(f, F / F.max(), c, lw=0.8, label=lab)
for n in range(1, 6):
    a.axvline(n * X.NU_POL_SH, color="0.8", lw=0.5)
a.set_xlim(0, 90); a.set_ylim(1e-4, 2); a.set_xlabel("frequency (THz)"); a.set_ylabel("|FFT| (norm.)")
a.set_title("(c) harmonics (grey: n × 16.2 THz measured)"); a.legend(fontsize=7)

a = ax[0, 3]
a.loglog(E, A1 / A1[i20], "o-", label=f"n=1 (slope {slopes[0]:.2f})")
a.loglog(E, A2 / A2[i20], "s-", label=f"n=2 (slope {slopes[1]:.2f})")
a.loglog(E, A3 / A3[i20], "^-", label=f"n=3 (slope {slopes[2]:.2f})")
for k in (1, 2, 3):
    a.loglog(E, (E / 20.0) ** k, ":", color="0.6")
a.set_xlabel("peak field (MV/cm)"); a.set_ylabel("harmonic amplitude (norm. at 20 MV/cm)")
a.set_title("(d) harmonic scaling (measured: 1, 2, 3)"); a.legend(fontsize=7)

a = ax[1, 0]
a.errorbar(X.REN_E, X.REN_W, yerr=X.REN_ERR, fmt="o", color="C3", label="measured (Ext. Data Fig. 6)")
a.plot(E, ren, "C0o-", label="model (fundamental of SH signal)")
rd = D.get("runsA_dft", [])
if rd:
    a.plot([r["E0"] for r in rd], np.array([r["f1"] for r in rd]) / rd[0]["f1"], "s--", color="0.5",
           label="same, Q_IR potential of the published DFT fit")
a.axhline(1, color="0.7", lw=0.5)
a.set_xlabel("peak field (MV/cm)"); a.set_ylabel("ω'/ω (low field)")
a.set_title("(e) frequency renormalization"); a.legend(fontsize=7, loc="lower left")

a = ax[1, 1]
if "A_20" in tr:
    tt = tr["A_20"]["t"] * 1e12
    Ec = tr["A_20"]["Ecoh"]
    a.plot(tt, Ec[:, 0], "C3", lw=0.8, label="coherent modes, surface cell (eV/cell)")
    a2 = a.twinx()
    a2.plot(tt, tr["A_20"]["Ts"][:, 0], "C1", label="bath T*, surface")
    a2.plot(tt, tr["A_20"]["Ts"][:, 2], "C1--", label="bath T*, 1.2 µm")
    a2.set_ylabel("bath temperature T* (K)")
    a2.legend(fontsize=6, loc="lower right")
a.set_xlabel("time (ps)"); a.set_ylabel("energy (eV / cell)"); a.set_xlim(-0.4, 1.4)
a.set_title("(f) coherent → incoherent, 20 MV/cm"); a.legend(fontsize=6, loc="upper left")

a = ax[1, 2]
cols = {19.0: "C3", 21.0: "C1", 22.0: "k"}
for nu, (Fx, Px) in X.REV.items():
    a.plot(Fx, Px, "o", color=cols[nu], ms=4, label=f"measured {nu:.0f} THz")
    rr = [r for r in runsB if r["nu"] == nu]
    a.plot([r["F"] for r in rr], [r["Pmin"] for r in rr], "-", color=cols[nu], label=f"model {nu:.0f} THz")
rbd = D.get("runsB_dft", [])
if rbd:
    a.plot([r["F"] for r in rbd], [r["Pmin"] for r in rbd], "x--", color="0.5", label="model 19 THz, DFT Q_IR potential")
a.axhline(0, color="0.7", lw=0.5)
a.set_xlabel("fluence (mJ/cm²)"); a.set_ylabel("(P_S/P_S0) at peak")
a.set_title("(g) rectified polarization (Mankowsky et al.)"); a.legend(fontsize=6)

a = ax[1, 3]
# linear check: c-axis reflectivity from MACE frequencies + three-phonon widths + DFPT polarities
from params import E_CH, AMU, EPS0
OMv = abs(np.linalg.det(np.load(os.path.join(HERE, "gamma_modes.npz"))["cell"])) * 1e-30
nu = np.linspace(5, 32, 1200); w = 2 * np.pi * nu * 1e12
gam_A1 = np.array([ph["bath_gamma"][0, 0, int(b)] for b in an["ia1"][:3]] + [ph["gam_to4"][0]])  # THz half widths
eps = EPS_INF + 0j
for j in range(4):
    S = (ZMODE[j] * E_CH) ** 2 / AMU / (EPS0 * OMv)
    wj = 2 * np.pi * an["freqA1"][j] * 1e12
    eps = eps + S / (wj ** 2 - w ** 2 - 1j * w * 2 * (2 * np.pi * gam_A1[j] * 1e12))
R_mod = np.abs((np.sqrt(eps) - 1) / (np.sqrt(eps) + 1)) ** 2
osc = [(249.3, 922.8, 27.7), (271.6, 384.1, 20.0), (632.0, 955.9, 33.5), (696.7, 352.5, 76.2)]
wc = nu * CM
eps_x = 4.4054 + sum(S ** 2 / (w0 ** 2 - wc ** 2 - 1j * g * wc) for w0, S, g in osc)
R_fit = np.abs((np.sqrt(eps_x) - 1) / (np.sqrt(eps_x) + 1)) ** 2
a.plot(nu, R_fit, color="C3", label="measured (4-oscillator fit, Ext. Data Table 1)")
a.plot(nu, R_mod, "C0", label="model: MACE ω, Γ; DFPT polarities; ε∞ = 4.6")
a.axvline(17.5, color="0.7", ls=":", lw=0.8); a.axvline(19.0, color="0.7", ls=":", lw=0.8)
a.set_xlabel("frequency (THz)"); a.set_ylabel("reflectivity (E ∥ c)")
a.set_title("(h) linear reflectivity, E ∥ c")
a.text(0.03, 0.04, f"TO4 FWHM 300 K: MACE {fwhm300:.0f} cm⁻¹; measured {X.TO4_FWHM_STOICH_293:.0f} (stoich.), {X.TO4_FWHM_CM:.0f} (congr.)",
       transform=a.transAxes, fontsize=6)
a.legend(fontsize=6)
# ---------------------------------------------------------------- controlled diagnostics
P1 = os.path.join(HERE, "diag_part1.pkl"); P2 = os.path.join(HERE, "diag_part2.pkl")
P3 = os.path.join(HERE, "diag_part3.pkl")
D3 = pickle.load(open(P3, "rb")) if os.path.exists(P3) else None
D1 = pickle.load(open(P1, "rb")) if os.path.exists(P1) else None
D2 = pickle.load(open(P2, "rb")) if os.path.exists(P2) else None
LAB = {"L1": "L1: Abalmasov first-principles part (c₂, double well)",
       "L2": "L2: L1 + DFT Q_IR anharmonicity (von Hoegen)",
       "L3": "L3: Abalmasov full (c₁, c₃, E_d fitted to this data)"}
COLL = {"L1": "C2", "L2": "C4", "L3": "C5"}
a = ax[2, 0]
Fx, Px = X.REV[19.0]
a.plot(Fx, Px, "o", color="C3", ms=4, label="measured 19 THz")
a.plot(F19, P19, "C0-", label="model, MACE surface (default)")
if D1:
    for k in ("L1", "L2"):
        rr = [r for r in D1[k] if r["nu"] == 19.0]
        a.plot([r["F"] for r in rr], [r["Pmin"] for r in rr], "s--", color=COLL[k], ms=3, label=LAB[k] + ", Maxwell")
if D3:
    for k, c in (("MACE", "C0"), ("L1", COLL["L1"]), ("L2", COLL["L2"])):
        rr = D3[k]
        a.plot([r["F"] for r in rr], [max(r["Pmin"], -1.2) for r in rr], ":", marker="v", ms=3, color=c,
               label=f"{k}, vacuum field on the mode (no Maxwell)")
a.axhline(0, color="0.7", lw=0.5)
a.set_xlabel("fluence (mJ/cm²)"); a.set_ylabel("(P_S/P_S0) at peak, 19 THz")
a.set_ylim(-1.25, 1.05)
a.set_title("(i) surfaces × field treatment, 19 THz"); a.legend(fontsize=5)

a = ax[2, 1]
def slope(rows, nu):
    rr = [r for r in rows if r["nu"] == nu and r["F"] <= 55]
    if not rr:
        rr = [r for r in rows if r["nu"] == nu][:1]
    Fm = np.array([r["F"] for r in rr]); Pm = np.array([r["Pmin"] for r in rr])
    return -np.sum(Fm * (Pm - 1)) / np.sum(Fm ** 2)
sets = [("measured", [dict(nu=nu, F=f, Pmin=p) for nu in X.REV for f, p in zip(*X.REV[nu])], "C3"),
        ("MACE", runsB, "C0")] + ([(k, D1[k], COLL[k]) for k in ("L1", "L2")] if D1 else [])
xs = np.arange(3)
for i, (name, rows, c) in enumerate(sets):
    sl = np.array([slope(rows, nu) for nu in (19.0, 21.0, 22.0)])
    a.bar(xs + i * 0.15 - 0.3, sl / sl[0], width=0.14, color=c, label=f"{name} (19 THz: {sl[0]:.4f}/mJcm⁻²)")
a.set_xticks(xs); a.set_xticklabels(["19 THz", "21 THz", "22 THz"])
a.set_ylabel("susceptibility, normalized to 19 THz"); a.set_title("(j) pump-frequency ordering")
a.legend(fontsize=5.5)

a = ax[2, 2]
if D2:
    for lab, c in (("mace_17.5", "C0"),):
        for row, ls in ((D2[lab][0], "-"), (D2[lab][-1], "-")):
            tr_ = row["tr"]
            f, F = __import__("drive").spectrum(tr_["t"], tr_["PR"], -0.3e-12, 1.2e-12, pad=32)
            a.semilogy(f, F / F.max(), color="C0" if row["E0"] < 5 else "C3", lw=0.8,
                       label=f"PR proxy, {row['E0']:.1f} MV/cm")
    a.axvline(15.3, color="C3", ls=":", lw=0.8, label="measured ν_p (PR) 15.3 THz")
    a.axvline(17.5, color="0.6", ls="--", lw=0.6, label="pump 17.5 THz")
    a.set_xlim(5, 60); a.set_ylim(1e-4, 2)
    a.set_xlabel("frequency (THz)"); a.set_ylabel("|FFT| (norm.)")
    a.set_title("(k) polarization-rotation signal (800 nm light line)"); a.legend(fontsize=6)

a = ax[2, 3]
a.errorbar(X.REN_E, X.REN_W, yerr=X.REN_ERR, fmt="o", color="C3", label="measured")
if D2:
    for lab, key, win, sty, txt in (("mace_17.5", "PR", "full", "C0o-", "PR polariton, full window"),
                                    ("mace_17.5", "SH", "full", "C1s-", "SH dominant, full window"),
                                    ("mace_17.5", "SH", "free", "C2^-", "SH, free window 0.25–1.2 ps"),
                                    ("mace_16.8", "PR", "full", "C0o:", "PR polariton, measured detuning"),
                                    ("dftIR_17.5", "SH", "full", "s--", "SH dominant, DFT Q_IR potential")):
        if lab not in D2:
            continue
        rows = D2[lab]
        E_ = np.array([r["E0"] for r in rows]); w_ = np.array([r[key][win]["dom"] for r in rows])
        kw = dict(color="0.5") if sty.startswith("s--") else {}
        a.plot(E_, w_ / w_[0], sty if not kw else "s--", ms=3, label=txt, **kw)
if D3 and "freq_mace" in D3:
    rr = D3["freq_mace"]
    a.plot([r["E0"] for r in rr], np.array([r["free"] for r in rr]) / rr[0]["free"], "v:", color="C0", ms=3,
           label="vacuum field on the mode, MACE, free window")
a.axhline(1, color="0.7", lw=0.5); a.set_ylim(0.9, 1.045)
a.set_xlabel("peak field (MV/cm)"); a.set_ylabel("ω'/ω (low field)")
a.set_title("(l) ω' extracted like for like"); a.legend(fontsize=5.5, loc="lower left")
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "figures", "phononics.png"), dpi=130)

# ---------------------------------------------------------------- summary
L = []
p = L.append
p("LiNbO3: light resonantly driving the A1(TO4) phonon (nonlinear phononics)")
p("All model numbers are predictions: MACE medium-omat-0 (phonons, anharmonic surface, fc3),")
p("DFPT mode polarities (Veithen, Gonze & Ghosez 2004), eps_inf = 4.6 (Barker & Loudon 1967).")
p("No parameter was adjusted to the measured nonlinear data.")
p("")
p("Linear properties")
p(f"  A1 TO frequencies (cm-1): MACE {', '.join(f'{x*CM:.0f}' for x in an['freqA1'])}; measured 252, 275, 332, 632")
p(f"  TO4 linewidth FWHM at 300 K: {fwhm300:.1f} cm-1 (meshes 6/8/10: {', '.join(f'{x:.1f}' for x in fwhm_mesh)}); "
  f"measured: {X.TO4_FWHM_STOICH_293} cm-1 (stoichiometric, Raman, Schwarz & Maier 1997), "
  f"{X.TO4_FWHM_CM} cm-1 (congruent, reflectivity fit of experiment A)")
p(f"  TO4 FWHM at 10 K: MACE {2*float(ph['gam_to4_10K'])*CM:.1f} cm-1; measured low-temperature limit (stoichiometric) {X.TO4_FWHM_STOICH_0} cm-1")
p(f"  TO4 bubble self-energy shift: {ph['shift'][0]:+.3f} THz at 300 K, {ph['shift'][1]:+.3f} THz at 500 K")
p("")
p("Anharmonic potential along Q_IR (eV, Å amu^1/2)")
p(f"  MACE: a3 = {an['a3']*1e3:.0f}, a4 = {an['a4']*1e3:.0f}, a5 = {an['a5']*1e3:.0f} meV; DFT (von Hoegen et al.): 1568, 901, 7.1 meV")
p(f"  excursions at U = 0.7 eV: MACE {qm_m:.2f} / +{qp_m:.2f}; DFT {qm_d:.2f} / +{qp_d:.2f}; measured (scaled onto DFT) -1.33 / +0.75 (U = 0.70, 0.66 eV)")
p(f"  measured points as published vs MACE: rms {rms_m1*1e3:.0f} meV (vs DFT curve: {rms_dft*1e3:.0f} meV)")
p(f"  one free scale factor (published procedure) onto MACE: B = {B_m:.2f}, rms {rms_m*1e3:.0f} meV; onto DFT: B = {rB_d.x:.2f}")
p(f"    i.e. the measured shape matches MACE only if all amplitudes are {B_m:.1f}x and energies {B_m**2:.1f}x those quoted")
for fn, keys in (("anharmonic.log", ("fc3 check", "g_TO")), ("model_check.log", ("nu ",))):
    fp = os.path.join(HERE, fn)
    if os.path.exists(fp):
        for line in open(fp, encoding="utf-8"):
            if any(k in line for k in keys):
                p("  " + line.strip())
p("")
p("Experiment A (von Hoegen et al., Nature 555, 79 (2018)): 17.5 THz, 150 fs")
p(f"  at 20 MV/cm: Q_IR at the surface spans [{r20['qmin']:.2f}, {r20['qmax']:.2f}] Å amu^1/2 (measured max excursion 1.4 on the soft side)")
p(f"  peak coherent energy, surface cell: {r20['Ecoh_max_surf']:.2f} eV/cell (measured 0.7 eV at max excursion; 0.6 eV deposited, estimate)")
p(f"  peak coherent energy, mean over first 1.3 µm: {r20['Ecoh_max_layer']:.2f} eV/cell")
p(f"  bath temperature T* at the surface, max: {r20['Tmax']:.0f} K; TO4 decay rate Gamma max {r20['G4max']:.3f} THz (300 K: {ph['gam_to4'][0]:.3f})")
p(f"  harmonic amplitudes vs field (5-21 MV/cm), power-law exponents n=1,2,3: {slopes[0]:.2f}, {slopes[1]:.2f}, {slopes[2]:.2f}; measured 1, 2, 3")
qr = r20["Q1rect"]
p(f"  rectified (cycle-averaged) displacements at the surface, 20 MV/cm (Å amu^1/2): TO1 {qr[0]:+.3f}, TO2 {qr[1]:+.3f}, "
  f"TO3 {qr[2]:+.3f}, TO4 {qr[3]:+.3f}; P_S change at peak (SH layer) {r20['Pmin']-1:+.3f} P_S0")
p(f"  fundamental of the SH signal at low field: {f1[0]:.2f} THz (measured polariton 16.2 THz, TO 19 THz; MACE TO 18.24 THz)")
p("  frequency renormalization omega'/omega:")
for e_, w_, mm in zip(X.REN_E, X.REN_W, ren_model):
    p(f"    {e_:5.1f} MV/cm: model {mm:.4f}   measured {w_:.3f} +- {X.REN_ERR}")
p("  sensitivity at 20 MV/cm (f1 THz, Q_IR range, coherent energy eV):")
for k, v in sens.items():
    p(f"    {k:20s} f1 {v['f1']:.3f}  Q [{v['qmin']:.2f}, {v['qmax']:.2f}]  E {v['Ecoh']:.2f}")
p("")
p("Experiment B (Mankowsky et al., PRL 118, 197601 (2017)): (P_S/P_S0) at peak")
for nu in (19.0, 21.0, 22.0):
    rr = [r for r in runsB if r["nu"] == nu]
    Fx, Px = X.REV[nu]
    for r in rr:
        meas = np.interp(r["F"], Fx, Px)
        p(f"  {nu:.0f} THz  {r['F']:6.1f} mJ/cm2 ({r['E0']:.1f} MV/cm): model {r['Pmin']:.3f}  measured {meas:.2f}")
slopes_B = {}
for nu in (19.0, 21.0, 22.0):
    rr = [r for r in runsB if r["nu"] == nu]
    Fm = np.array([r["F"] for r in rr]); Pm = np.array([r["Pmin"] for r in rr])
    Fx, Px = X.REV[nu]
    mm = Fx <= 55
    slopes_B[nu] = (-np.polyfit(Fm, Pm - 1, 1)[0], -np.polyfit(Fx[mm], Px[mm], 1)[0])
p("  susceptibility -d(P/P0)/dF (per mJ/cm^2; measured from points <= 55 mJ/cm^2): " + "; ".join(
    f"{nu:.0f} THz model {v[0]:.5f} measured {v[1]:.4f}" for nu, v in slopes_B.items()))
p("  normalized to 19 THz: model " + ", ".join(f"{slopes_B[n][0]/slopes_B[19.0][0]:.2f}" for n in (19.0, 21.0, 22.0))
  + "; measured " + ", ".join(f"{slopes_B[n][1]/slopes_B[19.0][1]:.2f}" for n in (19.0, 21.0, 22.0)))
if rbd:
    p("  with the published DFT Q_IR potential (19 THz): " + ", ".join(f"{r['F']:.1f} mJ/cm2 -> {r['Pmin']:.3f}" for r in rbd))
if rd:
    p("")
    p("Diagnostic: experiment A with the Q_IR self-potential of the published DFT fit (MACE couplings, MACE fc3):")
    for r in rd:
        p(f"  {r['E0']:5.1f} MV/cm: Q_IR [{r['qmin']:.2f}, {r['qmax']:.2f}], coherent energy {r['Ecoh_max_surf']:.2f} eV, "
          f"omega'/omega {r['f1']/rd[0]['f1']:.4f}")
p("")
p("Controlled diagnostics (diagnostics.py): same method, other energy surfaces / like-for-like readout")
def a_fit(E_, w_):
    return np.sum(E_ ** 2 * (w_ ** 2 - 1)) / np.sum(E_ ** 4)
p(f"  measured: fit sqrt(1 + a E^2) to Ext. Data Fig. 6 gives a = {a_fit(X.REN_E, X.REN_W):+.2e} (MV/cm)^-2")
if D1:
    p("  Part 1, (P_S/P_S0) at peak, 19 THz, fluence 12.8 / 51.8 / 94.5 mJ/cm2 (measured 0.79 / 0.17 / -0.35):")
    p(f"    MACE (default)  {' / '.join(f'{np.interp(F, F19, P19):.3f}' for F in (12.8, 51.8, 94.5))}")
    for k in ("L1", "L2", "L3"):
        rr = [r for r in D1[k] if r["nu"] == 19.0]
        Fk = np.array([r["F"] for r in rr]); Pk = np.array([r["Pmin"] for r in rr])
        thr = Fk[np.argmax(Pk < 0)] if np.any(Pk < 0) else None
        sl = [slope(D1[k], nu) for nu in (19.0, 21.0, 22.0)]
        p(f"    {LAB[k]}: {' / '.join(f'{np.interp(F, Fk, Pk):.3f}' for F in (12.8, 51.8, 94.5))}; "
          f"first sign change at {thr if thr is not None else 'none up to 94.5'} mJ/cm2; "
          f"susceptibility 19/21/22 THz = 1 : {sl[1]/sl[0]:.2f} : {sl[2]/sl[0]:.2f} (19 THz {sl[0]:.5f})")
    slm = [slope([dict(nu=nu, F=f, Pmin=pp) for nu in X.REV for f, pp in zip(*X.REV[nu])], nu) for nu in (19.0, 21.0, 22.0)]
    p(f"    measured: first sign change between 51.8 and 65.9 mJ/cm2 (threshold 60); susceptibility 1 : "
      f"{slm[1]/slm[0]:.2f} : {slm[2]/slm[0]:.2f} (19 THz {slm[0]:.4f})")
if D1:
    p("    note: L3 (Abalmasov's c1, c3, E_d, adjusted to experiment B in the source) is not a consistent surface: with c3 its")
    p("    stationary point moves to Q_P = 4.89, Q_IR = -0.90 (the source starts from the non-stationary Q_P0 = 3.1); listed only.")
if D3:
    p("  Part 3, field treatment: vacuum peak field applied directly to the modes (single cell, no reflection or screening,")
    p("  as in the minimal models of the literature); (P_S/P_S0) min at 12.8 / 51.8 / 65.9 / 94.5 mJ/cm2, max |Q_IR|:")
    for k in ("MACE", "L1", "L2"):
        rr = D3[k]; Fk = np.array([r["F"] for r in rr]); Pk = np.array([r["Pmin"] for r in rr])
        p(f"    {k:5s} {' / '.join(f'{np.interp(F, Fk, Pk):+.3f}' for F in (12.8, 51.8, 65.9, 94.5))}   |Q_IR| max at 12.8 mJ/cm2 "
          f"{rr[1]['q4max']:.2f} (Maxwell, surface cell: {max(abs(B19[2]['qmin']), B19[2]['qmax']):.2f})")
    for pot in ("mace", "dft"):
        rr = D3.get(f"freq_{pot}")
        if rr:
            p(f"    driven-mode frequency, vacuum field, {pot} Q_IR potential, free window: "
              + ", ".join(f"{r['E0']:.1f}: {r['free']/rr[0]['free']:.3f}" for r in rr) + f"  (|Q_IR| max {rr[-1]['qmax']:.2f} at 20.9)")
P4 = os.path.join(HERE, "diag_part4.pkl")
if os.path.exists(P4):
    D4 = pickle.load(open(P4, "rb"))
    p("  Part 4, Maxwell, pump on the MACE TO4 resonance (18.24 THz, same relative detuning as 19 THz in the experiment):")
    for k in ("MACE", "L1"):
        p(f"    {k:5s} (P_S/P_S0) min at 12.8 / 51.8 / 94.5 mJ/cm2: " + " / ".join(f"{r['Pmin']:.3f}" for r in D4[k])
          + f"; |Q_IR| max at 94.5: {D4[k][-1]['q4']:.2f}")
if D2:
    p("  Part 2, omega'/omega at 5.7 / 11.6 / 17.6 / 20.9 MV/cm (measured 1.004 / 1.012 / 1.020 / 1.036):")
    for lab, key, win, txt in (("mace_17.5", "PR", "full", "polarization rotation, polariton peak, full window"),
                               ("mace_17.5", "PR", "free", "polarization rotation, free window"),
                               ("mace_17.5", "SH", "full", "SH signal, dominant peak, full window"),
                               ("mace_17.5", "SH", "free", "SH signal, free window"),
                               ("mace_16.8", "PR", "full", "PR, pump at the measured detuning (16.8 THz)"),
                               ("mace_16.8", "SH", "full", "SH, pump at the measured detuning"),
                               ("dftIR_17.5", "PR", "full", "PR, DFT Q_IR potential"),
                               ("dftIR_17.5", "SH", "full", "SH, DFT Q_IR potential")):
        if lab not in D2:
            continue
        rows = D2[lab]
        E_ = np.array([r["E0"] for r in rows]); w_ = np.array([r[key][win]["dom"] for r in rows])
        w_ = w_ / w_[0]
        p(f"    {txt:52s} {' / '.join(f'{np.interp(e, E_, w_):.4f}' for e in X.REN_E)}   a = {a_fit(E_[1:], w_[1:]):+.2e}")
    r0 = D2["mace_17.5"][0]
    p(f"  low-field peaks: PR polariton {r0['PR']['full']['dom']:.2f} THz (measured 15.3); SH dominant "
      f"{r0['SH']['full']['dom']:.2f} THz, free window {r0['SH']['free']['dom']:.2f} THz (measured 16.2)")
fp = os.path.join(HERE, "isotope_check.log")
if os.path.exists(fp):
    p("  Linewidth: isotope (mass-disorder) scattering of TO4, natural abundances: " + open(fp).read().strip().splitlines()[-1])
p("")
bud = r20["budget"]
bal = (bud["Uref"] + bud["Uabs"] + bud["Ufield"] + bud["Umech"] + bud["Udiss"]) / bud["Uin"]
p(f"Energy budget at 20 MV/cm (fractions of incident): reflected {bud['Uref']/bud['Uin']:.3f}, coherent {bud['Umech']/bud['Uin']:.4f}, "
  f"dissipated to bath {bud['Udiss']/bud['Uin']:.4f}, absorbed deep {bud['Uabs']/bud['Uin']:.3f}, field {bud['Ufield']/bud['Uin']:.3f}; sum {bal:.4f}")
open(os.path.join(HERE, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
