"""MgB2 figures: gaps, linear conductivity, nonlinear response vs temperature, pumped state."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from twoband import fit_wc, gaps, sigma_total, THZ

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "..", "figures"); os.makedirs(FIG, exist_ok=True)
wc = fit_wc()
fig, ax = plt.subplots(2, 3, figsize=(14, 7))
# (a) gaps
T = np.linspace(2, 37.9, 60)
G = np.array([gaps(t, wc) for t in T])
ax[0, 0].plot(T, 2 * G[:, 0] / THZ, "b-", label=r"$2\Delta_\pi$ (two-band BCS)")
ax[0, 0].plot(T, 2 * G[:, 1] / THZ, "r-", label=r"$2\Delta_\sigma$")
dT = [5, 12, 20, 25, 30, 32, 34, 36]; dv = [1.10, 1.05, 0.87, 0.76, 0.68, 0.57, 0.45, 0.35]
ax[0, 0].plot(dT, dv, "bo", label=r"$2\Delta_\pi$ from $\sigma(\omega)$ (Katsumi)")
ax[0, 0].plot([5, 20, 28, 32, 37], [3.4, 3.3, 2.9, 2.4, 1.0], "r^", label=r"$2\Delta_\sigma$ (STM)")
ax[0, 0].set_xlabel("T (K)"); ax[0, 0].set_ylabel("THz"); ax[0, 0].legend(fontsize=7)
# (b) sigma2 at 5 K and sigma1 at 38 K
f = np.linspace(0.25, 3, 60)
s5 = np.array([sigma_total(x, gaps(5, wc), 5) for x in f]) / 100
ax[0, 1].plot(f, s5.imag, "b-", label=r"$\sigma_2$ 5 K (model)")
ax[0, 1].plot(f, s5.real, "b--", label=r"$\sigma_1$ 5 K (model)")
ax[0, 1].plot([0.3, 0.5, 1.0, 2.0, 3.0], [8e4, 6e4, 2e4, 1e4, 1e4], "bo", label=r"$\sigma_2$ 5 K (meas.)")
ax[0, 1].set_yscale("log"); ax[0, 1].set_ylim(1e2, 2e5); ax[0, 1].set_xlabel("f (THz)")
ax[0, 1].set_ylabel(r"$\Omega^{-1}$cm$^{-1}$"); ax[0, 1].legend(fontsize=7)
# (c) broadband peak
d = np.load(os.path.join(HERE, "nonlinear_out.npz"))
ax[0, 2].plot(d["Tb"], d["peaks"], "rs-", label="model: peak of nonlinear spectrum")
ax[0, 2].plot([7, 12, 17, 22, 27, 30, 32, 36, 37], [1.05, 1.03, 0.95, 0.9, 0.8, 0.7, 0.65, 0.6, 0.55], "ko",
              label="measured (Fig. 2a)")
ax[0, 2].plot(T, 2 * G[:, 0] / THZ, "b:", label=r"$2\Delta_\pi$")
ax[0, 2].set_ylim(0, 1.4); ax[0, 2].set_xlabel("T (K)"); ax[0, 2].set_ylabel("THz"); ax[0, 2].legend(fontsize=7)
# (d) FH, (e) TH normalized
eFH = ([7.5, 12, 17, 20, 23, 27.5, 31, 35], [1.0, 0.83, 0.68, 0.65, 0.5, 0.32, 0.2, 0.1])
eTH = ([8, 12, 17, 20, 22, 25, 28], [1.0, 0.9, 0.5, 0.3, 0.2, 0.1, 0.0])
for k, (key, keyc, e, name) in enumerate((("FH", "FHc", eFH, "first harmonic"), ("TH", "THc", eTH, "third harmonic"))):
    a = ax[1, k]
    a.plot(d["Ts"], d[keyc], "rs-", label="model, Raman-type (clean) vertex")
    a.plot(d["Ts"], d[key], "g^--", label="model, superfluid-weight (dirty) vertex")
    a.plot(*e, "ko", label="measured, 12 kV/cm (Fig. 4)")
    a.set_xlabel("T (K)"); a.set_ylabel(r"$I_{NL}/E_A^6$ (norm.)"); a.set_title(name); a.legend(fontsize=7)
    a.set_ylim(0, 2.6)
# (f) pumped state: quasiparticle-bottleneck (Owen-Scalapino) state vs thermal state
import bottleneck as bn, pump as pm
from scipy.optimize import brentq
U_abs = pm.absorbed(normal=False) * pm.VM
D0 = bn.gaps_mu(0.0); E0 = bn.qp_energy(0.0, D0)
mu = brentq(lambda m: bn.qp_energy(m, bn.gaps_mu(m)) - E0 - U_abs, 0.0, 6.0)
Dmu = bn.gaps_mu(mu)
ff = np.linspace(0.3, 2.0, 30)
sb = np.array([bn.sigma_mu(x, Dmu, mu) for x in ff]) / 100
Tst = 20.5
st = np.array([sigma_total(x, gaps(Tst, wc), Tst) for x in ff]) / 100
ax[1, 2].plot(ff, sb.real, "r-", label=r"$\sigma_1$, bottleneck state ($\mu^*$=%.1f meV)" % mu)
ax[1, 2].plot(ff, sb.imag, "b-", label=r"$\sigma_2$, bottleneck state")
ax[1, 2].plot(ff, st.real, "r:", label=r"$\sigma_1$, thermal T*=20.5 K")
ax[1, 2].plot(ff, st.imag, "b:", label=r"$\sigma_2$, thermal T*")
ax[1, 2].plot([0.8], [2.55e4], "ro", label=r"measured $\sigma_1$ (2 ps, 150 kV/cm)")
ax[1, 2].plot([0.5, 1.0], [2.6e4, 1.5e4], "bo", label=r"measured $\sigma_2$")
ax[1, 2].set_xlabel("f (THz)"); ax[1, 2].set_ylabel(r"$\Omega^{-1}$cm$^{-1}$"); ax[1, 2].legend(fontsize=6)
ax[1, 2].set_title("2 ps after 150 kV/cm pump, 7 K")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "mgb2.png"), dpi=140)
print("ok")
