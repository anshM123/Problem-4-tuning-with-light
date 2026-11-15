"""WSe2 figure and summary: valley-resolved exciton populations vs the measured traces."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from excitons import exciton_table
from dynamics import PARAMS, run, run_above_gap

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "..", "figures"); os.makedirs(FIG, exist_ok=True)
tab = exciton_table()
P = dict(PARAMS)
P["dE_SK"] = P["dCB_SK"] - (tab["Eb_SK"] - tab["Eb_KK"])
P["M_KK"], P["M_SK"] = tab["M_KK"], tab["M_SK"]
exp_t = np.array([-0.5, -0.2, 0.09, 0.4, 0.7, 1.0, 1.3, 1.6, 2.1, 2.6, 3.1, 3.6])
exp_K = np.array([0.1, 0.4, 7.5, 6.8, 5.1, 4.2, 3.5, 3.2, 2.5, 1.75, 1.65, 1.35])
exp_S = np.array([0.1, 0.25, 6.55, 9.8, 9.55, 8.95, 7.8, 6.55, 5.25, 4.45, 3.35, 2.75])

res = run(P, "resonant")
hot = run_above_gap(P, tau_form=0.5e-12)
t = res["t"] * 1e12
K = res["KK"] + res["KKp"]; S = res["SK"]
scale = exp_K.max() / K.max()
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
ax[0].plot(t, K * scale, "b-", label="K valley (model)")
ax[0].plot(t, S * scale, "r-", label="Sigma valley (model)")
ax[0].plot(exp_t, exp_K, "bo", ms=4, label="K (Madeo et al.)")
ax[0].plot(exp_t, exp_S, "ro", ms=4, label="Sigma (Madeo et al.)")
ax[0].set_xlabel("delay (ps)"); ax[0].set_ylabel("population (a.u., one common scale)")
ax[0].set_title("resonant 1.72 eV pump, 90 K"); ax[0].legend(fontsize=7)
ax[1].plot(t, K / np.maximum(S, 1e-30), "k-", label="K / Sigma (model)")
ax[1].plot(exp_t[2:], exp_K[2:] / exp_S[2:], "ko", ms=4, label="measured")
ax[1].set_ylim(0, 3.5); ax[1].set_xlim(0, 4); ax[1].set_xlabel("delay (ps)"); ax[1].set_title("bright/dark ratio")
ax[1].legend(fontsize=7)
th = hot["t"] * 1e12
ax[2].plot(th, hot["K"] / np.maximum(hot["S"], 1e-30), "k-", label="K / Sigma (model: carriers + excitons)")
ax[2].plot([0.1, 3.5], [0.15, 0.32], "ks", label="measured (start, late)")
ax[2].set_ylim(0, 1); ax[2].set_xlim(0, 4); ax[2].set_xlabel("delay (ps)"); ax[2].set_title("above-gap 2.48 eV pump")
ax[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "wse2_valleys.png"), dpi=150)

# summary numbers
iS = np.argmax(S); iK = np.argmax(K)
def t_at(frac, y):
    k = np.where(y >= frac * y.max())[0][0]; return t[k]
rise_S = t_at(0.5, S) - t_at(0.5, K)
late = (t > 1.0) & (t < 3.6)
ratio_late = np.mean(K[late] / S[late])
slope = np.polyfit(t[late], np.log(S[late]), 1)[0]
with open(os.path.join(HERE, "summary.txt"), "w") as fh:
    fh.write(f"E(KK) reference = 1.73 eV (A resonance); E(SK)-E(KK) = {1e3*P['dE_SK']:.1f} meV "
             f"(binding energies {1e3*tab['Eb_KK']:.0f} / {1e3*tab['Eb_SK']:.0f} meV)\n")
    fh.write(f"Sigma peak at {t[iS]:.2f} ps; K peak at {t[iK]:.2f} ps (measured maxima at the samples "
             f"0.4 and 0.09 ps)\n")
    # timing against the measured points themselves (the data are sampled at -0.2, 0.09, 0.4,
    # 0.7 ps, so rise and peak times between samples are not resolved by the measurement)
    fh.write(f"Sigma at 0.09 ps / Sigma peak: {np.interp(0.09, t, S)/S.max():.2f} "
             f"(measured {exp_S[2]/exp_S.max():.2f})\n")
    fh.write(f"K/Sigma at 0.09 ps: {np.interp(0.09, t, K)/np.interp(0.09, t, S):.2f} (measured {exp_K[2]/exp_S[2]:.2f}); "
             f"at 0.4 ps: {np.interp(0.4, t, K)/np.interp(0.4, t, S):.2f} (measured {exp_K[3]/exp_S[3]:.2f})\n")
    mm = exp_t > 0
    a = np.r_[np.interp(exp_t[mm], t, K), np.interp(exp_t[mm], t, S)]; b = np.r_[exp_K[mm], exp_S[mm]]
    c = (a @ b) / (a @ a)
    fh.write(f"rms deviation from all measured K and Sigma points (one common scale): "
             f"{np.sqrt(np.mean((c * a - b) ** 2)) / exp_S.max() * 100:.0f} % of the Sigma peak\n")
    fh.write(f"half-rise delay Sigma vs K: {rise_S*1e3:.0f} fs (not resolved by the measured sampling)\n")
    # homogeneous linewidth implied by the couplings (out-scattering rate of Q = 0 bright excitons)
    from dynamics import build
    PP = dict(P); PP["T"] = 50.0; mdl = build(PP)
    G = sum(r for (a_, b_, r, sh) in mdl["W"] if a_ == 0 and sh >= 0)
    fh.write(f"phonon linewidth (FWHM) of the bright exciton at 50 K: {G * 1e-12 * 0.6582:.1f} meV "
             f"(measured phonon part 6.0 meV, Moody et al., Nat. Commun. 6, 8315 (2015))\n")
    fh.write(f"mean K/Sigma over 1-3.6 ps: {ratio_late:.2f} (measured ~0.47)\n")
    fh.write(f"late decay time of Sigma: {-1/slope:.2f} ps (measured ~2.2-2.5 ps)\n")
    fh.write(f"peak Sigma / peak K: {S.max()/K.max():.2f} (measured 1.3)\n")
    rr = hot["K"] / np.maximum(hot["S"], 1)
    fh.write("above-gap K/Sigma at 0.1, 1, 3.5 ps: " + ", ".join(f"{np.interp(x, th, rr):.2f}" for x in (0.1, 1.0, 3.5))
             + " (measured 0.15 -> 0.32)\n")
print(open(os.path.join(HERE, "summary.txt")).read())
