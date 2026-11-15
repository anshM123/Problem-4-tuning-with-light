"""Self-check of run_material.py on graphite: compares the general pipeline (universal tight-binding
electron-phonon coupling, no material-specific calibration) with the dedicated graphite model of
materials/graphite (DFT-calibrated pi-band coupling with GW enhancement). Identical post-processing
is applied to both dynamics files.

    python materials/blind/graphite/compare_with_repo.py
"""
import json, os
import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..", "..", "graphite")
THz2eV = 4.135667696e-3
KB = 8.617333262e-5


def load(path):
    d = np.load(path)
    return d["t"], d["n"].astype(float), d["freq"], d["qfrac"], d["Te"]


def a1k(t, n, freq, qf):
    iK = [i for i, x in enumerate(qf) if np.allclose(np.asarray(x)[:2], [1 / 3, 1 / 3])][0]
    y = n[:, iK, 5] - n[0, iK, 5]
    pk = np.argmax(y)
    rise = t[np.where(y >= (1 - np.exp(-1)) * y[pk])[0][0]]
    rel = (y[pk:] - y[-1]) / (y[pk] - y[-1])
    dec = t[pk + np.where(rel <= np.exp(-1))[0][0]] - t[pk]
    # 1/e decay of the excess itself (repo definition: decays to 1/e of its peak excess)
    k = np.where(y[pk:] <= np.exp(-1) * y[pk])[0]
    dec_abs = t[pk + k[0]] - t[pk] if len(k) else np.nan
    m10 = (t >= -0.1) & (t <= 10)
    fK = lambda x, A, tr, td, C: np.where(x > 0, A * (1 - np.exp(-x / tr)) * np.exp(-x / td) + C * (1 - np.exp(-x / td)), 0)
    p, _ = curve_fit(fK, t[m10], y[m10], p0=[y.max(), 0.2, 2, y[-1]], maxfev=20000)
    return dict(rise63_ps=float(rise), peak_ps=float(t[pk]), decay_to_1e_of_peak_ps=float(dec_abs),
                decay_1e_above_final_ps=float(dec), fit_rise_ps=float(p[1]), fit_decay_ps=float(p[2]),
                peak_excess=float(y[pk]))


def branch_energy(t, n, freq, times):
    hw = freq * THz2eV
    E = (n * hw[None]).sum(1) / n.shape[1]
    E = (E - E[0]) * 1e3
    return {f"{x:g}ps": dict(acoustic_meV=float(E[np.argmin(abs(t - x)), :3].sum()),
                             optical_meV=float(E[np.argmin(abs(t - x)), 3:].sum())) for x in times}


out = {}
for lab, path in (("dedicated_graphite_model", os.path.join(REPO, "dynamics_out.npz")),
                  ("general_pipeline", os.path.join(HERE, "dynamics.npz"))):
    t, n, freq, qf, Te = load(path)
    out[lab] = dict(Te_peak_K=float(Te.max()), A1K=a1k(t, n, freq, qf),
                    energy=branch_energy(t, n, freq, [0.5, 1.5, 5, 50]))
R = json.load(open(os.path.join(HERE, "results.json")))
out["general_pipeline"]["absorbed_meV_per_cell"] = R["optics"]["absorbed_energy_meV_per_cell"]
out["general_pipeline"]["diffuse_K_fit_rise_decay"] = R["diffuse"]["K"]["fit_rise_decay"]
fits = np.load(os.path.join(REPO, "diffuse_fits.npz"))
out["dedicated_graphite_model"]["diffuse_K_fit"] = dict(tau_rise_ps=float(fits["K_rise"]), tau_decay_ps=float(fits["K_decay"]))
out["measured_Stern2018"] = dict(K_rise_ps=0.28, K_decay_ps=1.7, A1K_peak_ps="0.6-0.7")
json.dump(out, open(os.path.join(HERE, "comparison.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
