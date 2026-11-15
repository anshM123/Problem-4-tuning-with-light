"""figures/unified_limits.png: the unified solver against exact dynamics (a) and its three limits
(b: Level A, c: Level B, d: Level C), from the data written by the scripts in this folder."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "data")
FIG = os.path.join(HERE, "..", "..", "figures")
INK, BLUE, ORANGE, AQUA, VIOLET = "#0b0b0b", "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#52514e", "axes.labelcolor": INK, "xtick.color": "#52514e",
                     "ytick.color": "#52514e", "legend.frameon": False, "lines.linewidth": 1.4})
fig, ax = plt.subplots(2, 4, figsize=(13, 6.2))

# (a) full level vs exact (Hubbard-Holstein dimer, cubic phonons, light pulse at t = 8)
F = np.load(os.path.join(D, "full_vs_ed.npz"))
t = F["t"]
for j, (key, q, yl, ttl) in enumerate((("U0_g0.2", 1, r"$E_{opt}$ (J)", "(a) U=0, g=0.2: optical phonon energy"),
                                        ("U2_g0.1", 0, "double occupancy D", "(a) U=2J, g=0.1: all interactions"),
                                        ("U2_g0.2", 1, r"$E_{opt}$ (J)", "(a) U=2J, g=0.2: closure drifts after the pulse"))):
    a = ax[0, j]
    a.plot(t, F[key + "_ed"][q], color=INK, lw=2.0, label="exact (ED)")
    a.plot(t, F[key + "_full"][q], color=BLUE, label="unified, all correlations")
    a.plot(t, F[key + "_ehrenfest"][q], color=ORANGE, ls="--", label="same e-e, Ehrenfest phonons")
    a.plot(t, F[key + "_meanfield"][q], color=AQUA, ls=":", label="TDHF + Ehrenfest")
    a.set_title(ttl, fontsize=8, loc="left")
    a.set_xlabel(r"t ($\hbar$/J)"); a.set_ylabel(yl)
    a.axvspan(4, 12, color="#e9e8e4", lw=0)
fig.legend(*ax[0, 0].get_legend_handles_labels(), loc="lower center", ncol=4, fontsize=7.5,
           bbox_to_anchor=(0.37, 0.0))

# (b) Level A limit
A = np.load(os.path.join(D, "levelA_limit.npz"))
a = ax[0, 3]
for mode, c, lab in (("factorized", ORANGE, "T-matrix ladder"), ("sic", BLUE, "self-interaction-free")):
    a.plot(A[mode + "_t"], A[mode + "_I_levelA"], color=c, lw=3, alpha=0.35)
    a.plot(A[mode + "_t"], A[mode + "_I_unified"], color=c, lw=1, label=f"{lab}: unified (thin) = Level A (thick)")
    dev = np.abs(A[mode + "_I_unified"] - A[mode + "_I_levelA"]).max()
    a.text(0.98, 0.95 if mode == "sic" else 0.86, f"{lab}: max dev {dev:.0e}", transform=a.transAxes,
           ha="right", fontsize=7, color=INK)
a.set_title("(b) Level A limit: 11-site CDW, U=5J, no phonons", fontsize=8, loc="left")
a.set_xlabel(r"t ($\hbar$/J)"); a.set_ylabel("density imbalance I")
a.legend(fontsize=6.5, loc="upper right", bbox_to_anchor=(1.0, 0.82))

# (c) Level B limit
B = np.load(os.path.join(D, "levelB_limit_12.npz"))
a = ax[1, 0]
tt = B["p2_t"]
a.plot(tt, B["p2_graphite"][1], color=INK, lw=2.5, alpha=0.5, label="Level B code (graphite rate functions)")
a.plot(tt, B["p2_unified"][1], color=BLUE, label="unified, Markov + random-phase limit")
a.plot(tt, B["p2_graphite"][2:].sum(0) / 10, color=INK, lw=2.5, alpha=0.5, ls="--")
a.plot(tt, B["p2_unified"][2:].sum(0) / 10, color=BLUE, ls="--", label="acoustic total / 10")
a.set_title("(c) Level B code vs kinetic limit", fontsize=8, loc="left")
a.set_xlabel("t (1/w0)"); a.set_ylabel("phonon occupation"); a.legend(fontsize=6.5)
P3 = np.load(os.path.join(D, "levelB_limit_3.npz"))
for j, kind, lab in ((1, "pp", "anharmonic decay, 24 modes"), (2, "ep", "heating by hot electrons")):
    a = ax[1, j]
    show = {"pp": ("0.3", "0.2", "0.14"), "ep": ("0.5", "0.35", "0.25")}[kind]   # all couplings: data/levelB_part3.log
    cols = [ORANGE, AQUA, VIOLET, BLUE]
    for c, sc in zip(cols, show):
        p = f"p3{kind}_{sc}_"
        r = float(P3[p + "rate"]); ts = P3[p + "t"] * r
        a.plot(ts, P3[p + "full"][:, 0], color=c, lw=1, label=f"full, coupling x{sc}")
        a.plot(ts, P3[p + "kin"][:, 0], color=c, lw=2.2, alpha=0.45, ls="--")
    a.plot([], [], color="#888", lw=2.2, alpha=0.6, ls="--", label="kinetic limit (same couplings)")
    a.set_title(f"(c) full vs kinetic: {lab}", fontsize=8, loc="left")
    a.axhline(0, color="#c3c2b7", lw=0.6)
    a.set_xlabel("t x initial kinetic rate"); a.set_ylabel("optical occupation N0"); a.legend(fontsize=6.5)

# (d) Level C limit
C = np.load(os.path.join(D, "levelC_limit.npz"))
a = ax[1, 3]
w = C["w"]; win = w <= 15.0                       # probe-pulse window
for b, c, lab in ((0, BLUE, r"$\pi$"), (1, ORANGE, r"$\sigma$")):
    a.plot(w[win], np.abs(C[f"RlevelC{b}"][b][win]), color=c, lw=3, alpha=0.35)
    a.plot(w[win], np.abs(C[f"Rsolver{b}"][b][win]), color=c, lw=1, label=f"{lab} band: unified (thin), Level C (thick)")
for g in C["gap_solver"]:
    a.axvline(2 * g, color="#888", lw=0.8, ls=":")
a.set_yscale("log")
a.set_title("(d) Level C: amplitude response, MgB2, 5 K", fontsize=8, loc="left")
a.set_xlabel(r"$\hbar\omega$ (meV)   (dotted: $2\Delta_\pi$, $2\Delta_\sigma$ from the solver)")
a.set_ylabel(r"|$\delta\Delta_a/h_a$|"); a.legend(fontsize=6.5, loc="lower center")
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(os.path.join(FIG, "unified_limits.png"), dpi=150)
print("written figures/unified_limits.png")
