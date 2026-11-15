"""WSe2 figure and summary: coherent + incoherent exciton kinetics against the measured traces.

Resonant case: excitonic Bloch equations (bloch.py), coherent polarization at Q = 0 plus
energy-resolved incoherent occupations, DFPT exciton-phonon couplings.  The simulated signal is
passed through the measurement operator of Madeo et al. (sampling at their delays, their
definition of the timing) before comparison.  Above-gap case: free-carrier stage (dynamics.py).
"""
import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from excitons import exciton_table
from dynamics import PARAMS, run_above_gap
import bloch as B

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "..", "figures"); os.makedirs(FIG, exist_ok=True)
t_start = time.time()
tab = exciton_table()
P = dict(PARAMS)
P["dE_SK"] = P["dCB_SK"] - (tab["Eb_SK"] - tab["Eb_KK"])
P["M_KK"], P["M_SK"] = tab["M_KK"], tab["M_SK"]
exp_t, exp_K, exp_S = B.EXP_T, B.EXP_K, B.EXP_S

res = B.run_bloch(P)                       # coherent polarization + incoherent excitons
inc = B.run_bloch(P, coherent=False)       # rate-equation limit (incoherent injection at Q = 0)
hot = run_above_gap(P, tau_form=0.5e-12)
t = res["t"] * 1e12
K, S, P2 = res["K"], res["S"], res["P2"]
meas = B.as_measured(t, K, S)
pos = exp_t > 0
a = np.r_[meas["Ks"][pos], meas["Ss"][pos]]; b = np.r_[exp_K[pos], exp_S[pos]]
c = (a @ b) / (a @ a)                      # one common scale for both valleys
rms = np.sqrt(np.mean((c * a - b) ** 2)) / exp_S.max()
fit_d = B.fit_two_level(exp_t, exp_K, exp_S, sigma=0.5)
fit_m = B.fit_two_level(exp_t, c * meas["Ks"], c * meas["Ss"])
modes = B.relaxation_modes(P)
rt50 = B.rates_table(P, 50.0); rt90 = B.rates_table(P, 90.0)

fig, ax = plt.subplots(1, 4, figsize=(17, 3.7))
ax[0].plot(t, K * c, "b-", label="K valleys (model)")
ax[0].plot(t, S * c, "r-", label="Sigma valleys (model)")
ax[0].plot(t, P2 * c, "b:", lw=1, label="coherent part |p|^2")
ax[0].plot(exp_t, exp_K, "bo", ms=4, label="K (Madeo et al.)")
ax[0].plot(exp_t, exp_S, "ro", ms=4, label="Sigma (Madeo et al.)")
ax[0].set_xlabel("delay (ps)"); ax[0].set_ylabel("population (a.u., one common scale)")
ax[0].set_title("resonant 1.72 eV pump, 90 K"); ax[0].legend(fontsize=7)
for y, col, d in ((meas["Ks"], "b", exp_K), (meas["Ss"], "r", exp_S)):
    ax[1].plot(exp_t, y / y.max(), col + "-", lw=1)
    ax[1].plot(exp_t, d / d.max(), col + "o", ms=4)
ax[1].axvline(meas["t_Kmax"], color="b", ls=":"); ax[1].axvline(meas["t_Smax"], color="r", ls=":")
ax[1].set_xlim(-0.5, 1.05); ax[1].set_xlabel("delay (ps)"); ax[1].set_ylabel("normalized")
ax[1].set_title("as measured (Fig. S11): lines model, dots data", fontsize=9)
ax[2].plot(t, K / np.maximum(S, 1e-30), "k-", label="K / Sigma (model)")
ax[2].plot(exp_t[2:], exp_K[2:] / exp_S[2:], "ko", ms=4, label="measured")
ax[2].set_ylim(0, 3.5); ax[2].set_xlim(0, 4); ax[2].set_xlabel("delay (ps)"); ax[2].set_title("bright/dark ratio")
ax[2].legend(fontsize=7)
th = hot["t"] * 1e12
ax[3].plot(th, hot["K"] / np.maximum(hot["S"], 1e-30), "k-", label="K / Sigma (model: carriers + excitons)")
ax[3].plot([0.1, 3.5], [0.15, 0.32], "ks", label="measured (start, late)")
ax[3].set_ylim(0, 1); ax[3].set_xlim(0, 4); ax[3].set_xlabel("delay (ps)"); ax[3].set_title("above-gap 2.48 eV pump")
ax[3].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "wse2_valleys.png"), dpi=150)

iS = np.argmax(S); iK = np.argmax(K)
late = (t > 1.0) & (t < 3.6)
ratio_late = np.mean(K[late] / S[late])
slope = np.polyfit(t[late], np.log(S[late]), 1)[0]
def half(y):
    return t[np.where(y >= 0.5 * y.max())[0][0]]
slow = np.isfinite(modes["tau"]) & (modes["tau"] > 0.1e-12) & (modes["tau"] < 1e-9)
tau_slow = (modes["amp"][slow] * modes["tau"][slow]).sum() / modes["amp"][slow].sum()
ti = inc["t"] * 1e12; mi = B.as_measured(ti, inc["K"], inc["S"])
lines = [
    f"E(KK) reference = 1.73 eV (A resonance); E(SK)-E(KK) = {1e3*P['dE_SK']:.1f} meV "
    f"(binding energies {1e3*tab['Eb_KK']:.0f} / {1e3*tab['Eb_SK']:.0f} meV)",
    "-- timing as Madeo et al. define it (sampled at their delays; Fig. 3d,e and Fig. S11) --",
    f"Sigma maximum among the measured delays: {meas['t_Smax']:.2f} ps (measured 0.40 ps, the '~400 fs')",
    f"delay between sampled K and Sigma maxima: {meas['delay']:.2f} ps (measured 0.31 ps, the '~300 fs' of SI 9)",
    f"effective K->Sigma time of a two-level fit at the measured delays: model {1e3*fit_m['tau_f']:.0f} fs, "
    f"data {1e3*fit_d['tau_f']:.0f} fs (1-sigma {1e3*fit_d['tau_f_lo']:.0f}-{1e3*fit_d['tau_f_hi']:.0f} fs, 0.5 a.u. reading error)",
    "-- continuous-time model quantities (not resolved by the 0.3 ps delay step) --",
    f"Sigma peak at {t[iS]:.2f} ps; K peak at {t[iK]:.2f} ps; half-rise delay Sigma vs K {1e3*(half(S)-half(K)):.0f} fs",
    f"Q = 0 bright exciton out-scattering at 90 K: {rt90['G_Q0']*1e-12:.1f} /ps ({1e3/(rt90['G_Q0']*1e-12):.0f} fs), "
    f"{100*rt90['G_Q0_to_S']/rt90['G_Q0']:.0f} % into Sigma-K (phonon absorption)",
    f"thermal K->Sigma and Sigma->K rates at 90 K: {rt90['thermal_K_to_S']*1e-12:.0f} and {rt90['thermal_S_to_K']*1e-12:.0f} /ps",
    f"slow relaxation modes of the K/Sigma imbalance (energy relaxation of the cold resonantly created gas): "
    f"amplitude-weighted {1e15*tau_slow:.0f} fs",
    f"coherent excitons: peak |p|^2 / peak K signal = {P2.max()/K.max():.2f}; radiated coherently "
    f"{res['n_rad_coh']/P['n_exc']*100:.0f} % of the absorbed density",
    f"rate-equation limit (no coherent stage): sampled maxima {mi['t_Kmax']:.2f} / {mi['t_Smax']:.2f} ps, "
    f"K/Sigma(0.09 ps) {mi['Ks'][2]/mi['Ss'][2]:.2f}",
    "-- comparison with the measured points --",
    f"Sigma at 0.09 ps / Sigma at 0.4 ps: {meas['Ss'][2]/meas['Ss'][3]:.2f} (measured {exp_S[2]/exp_S[3]:.2f})",
    f"K/Sigma at 0.09 ps: {meas['Ks'][2]/meas['Ss'][2]:.2f} (measured {exp_K[2]/exp_S[2]:.2f}); "
    f"at 0.4 ps: {meas['Ks'][3]/meas['Ss'][3]:.2f} (measured {exp_K[3]/exp_S[3]:.2f})",
    f"rms deviation from all measured K and Sigma points (one common scale): {100*rms:.1f} % of the Sigma peak",
    f"phonon linewidth (FWHM) of the bright exciton at 50 K: {rt50['G_Q0']*0.6582e-12:.1f} meV "
    f"(measured phonon part 6.0 meV, Moody et al., Nat. Commun. 6, 8315 (2015)); at 90 K: {rt90['G_Q0']*0.6582e-12:.1f} meV",
    f"mean K/Sigma over 1-3.6 ps: {ratio_late:.2f} (measured ~0.47)",
    f"late decay time of Sigma (1-3.6 ps): {-1/slope:.2f} ps (measured ~2.2-2.5 ps); two-level fit decay: "
    f"model {fit_m['tau_d']:.2f} ps, data {fit_d['tau_d']:.2f} ps",
    f"peak Sigma / peak K: {S.max()/K.max():.2f} (measured 1.3)",
]
rr = hot["K"] / np.maximum(hot["S"], 1)
lines.append("above-gap K/Sigma at 0.1, 1, 3.5 ps: " + ", ".join(f"{np.interp(x, th, rr):.2f}" for x in (0.1, 1.0, 3.5))
             + " (measured 0.15 -> 0.32)")
with open(os.path.join(HERE, "summary.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
print(open(os.path.join(HERE, "summary.txt")).read())
print(f"runtime {time.time()-t_start:.0f} s")
