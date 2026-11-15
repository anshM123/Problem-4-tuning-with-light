"""Compare the end-to-end unified run (unified_out.npz) with the Level B results
(materials/graphite/dynamics_out.npz: thermal electrons; dynamics_neq.npz: non-thermal carriers).
The same code computes, for every run:
  * A1'(K) excess population: 63 % rise time, peak time, 1/e decay time after the peak;
  * one-phonon diffuse intensity at K, K/2, M/2, M and the experiment-style fits of
    materials/graphite/diffuse.py (K: rise x decay to a plateau; others: single-exponential rise);
  * excess energy per branch character (TO, LO, LA, TA, ZO, ZA) at 0.5, 1.5, 5, 100 ps; Te.
Writes comparison.txt and graphite_unified.png (this folder).
Usage: python analysis.py [unified_out*.npz ...]
"""
import os, sys
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

HERE = os.path.dirname(os.path.abspath(__file__))
GDIR = os.path.join(HERE, "..", "graphite")
RUNS = {"Level B (thermal e)": os.path.join(GDIR, "dynamics_out.npz"),
        "Level B (non-thermal e)": os.path.join(GDIR, "dynamics_neq.npz")}
NAMES = {"unified_out.npz": "unified", "unified_out_matched.npz": "unified, field x1.09",
         "unified_out_localU.npz": "unified, local-U e-e"}
for fn in (sys.argv[1:] or ["unified_out.npz", "unified_out_matched.npz", "unified_out_localU.npz"]):
    RUNS[NAMES.get(fn, fn)] = os.path.join(HERE, fn)
kB = 8.617333e-5
d0 = np.load(RUNS["Level B (thermal e)"])
freq, qf = d0["freq"], d0["qfrac"]
hw = freq * 4.135667e-3
Nq = len(qf)
a_ = 2.464
cell = [[a_, 0, 0], [-a_ / 2, a_ * np.sqrt(3) / 2, 0], [0, 0, 15.0]]
ph = Phonopy(PhonopyAtoms(symbols=["C", "C"], cell=cell, scaled_positions=[[1/3, 2/3, .5], [2/3, 1/3, .5]]),
             supercell_matrix=[8, 8, 1])
ph.force_constants = h5py.File(os.path.join(GDIR, "fc2.hdf5"))["force_constants"][:]
ph.run_qpoints(np.c_[qf, np.zeros(Nq)], with_eigenvectors=True)
ev = ph.qpoints.eigenvectors
B = 2 * np.pi * np.linalg.inv(np.array(cell)[:2, :2]).T
qc = qf @ B
# branch characters (as materials/graphite/plots.py)
label = np.empty((Nq, 6), object)
for i in range(Nq):
    qhat = qc[i] / (np.linalg.norm(qc[i]) + 1e-12)
    zs = [j for j in range(6) if np.sum(abs(ev[i][:, j].reshape(2, 3)[:, 2]) ** 2) > 0.5]
    ins = [j for j in range(6) if j not in zs]
    zs = sorted(zs, key=lambda j: freq[i, j]); ins = sorted(ins, key=lambda j: freq[i, j])
    if len(zs) == 2:
        label[i, zs[0]], label[i, zs[1]] = "ZA", "ZO"
    for pair, suf in ((ins[:2], "A"), (ins[2:], "O")):
        if len(pair) < 2:
            for j in pair:
                label[i, j] = "L" + suf
            continue
        lon = [np.sum(abs(ev[i][:, j].reshape(2, 3)[:, :2] @ qhat) ** 2) for j in pair]
        jl = pair[int(np.argmax(lon))]; jt = pair[1 - int(np.argmax(lon))]
        label[i, jl], label[i, jt] = "L" + suf, "T" + suf
iK = [i for i, x in enumerate(qf) if np.allclose(x, [1/3, 1/3])][0]
# diffuse-scattering weights (as materials/graphite/diffuse.py)
DT_a = [0.7307, 1.1951, 0.4563, 0.1247]; DT_b = [36.9951, 11.2966, 2.8139, 0.3456]
fC = lambda Q: sum(a * np.exp(-b * (Q / (4 * np.pi)) ** 2) for a, b in zip(DT_a, DT_b))
rfrac = np.array([[1/3, 2/3], [2/3, 1/3]])
pts, lab_q, lab_G = [], [], []
for g1 in range(-4, 5):
    for g2 in range(-4, 5):
        G = np.array([g1, g2])
        for iq in range(Nq):
            Q = (qf[iq] + G) @ B
            if 0.3 < np.linalg.norm(Q) < 7.0:
                pts.append(Q); lab_q.append(iq); lab_G.append(G)
pts = np.array(pts); lab_q = np.array(lab_q); lab_G = np.array(lab_G); absQ = np.linalg.norm(pts, axis=1)
Wd = np.zeros((len(pts), 6))
for k in range(len(pts)):
    iq = lab_q[k]; Q = pts[k]; G = lab_G[k]
    for j in range(6):
        if freq[iq, j] < 0.1:
            continue
        e = ev[iq][:, j].reshape(2, 3)[:, :2]
        F = sum(fC(absQ[k]) * (Q @ e[s]) * np.exp(-2j * np.pi * (G @ rfrac[s])) for s in range(2))
        Wd[k, j] = abs(F) ** 2 / freq[iq, j]
fK = lambda x, A, tr_, td, C: np.where(x > 0, A * (1 - np.exp(-x / tr_)) * np.exp(-x / td) + C * (1 - np.exp(-x / td)), 0)
f1 = lambda x, A, tau, t0: np.where(x > t0, A * (1 - np.exp(-(x - t0) / tau)), 0)


def metrics(fn):
    d = np.load(fn)
    t, n = d["t"], d["n"].astype(float)
    out = {}
    x = n[:, iK, 5] - n[0, iK, 5]
    ip = int(np.argmax(x)); pk = x[ip]
    out["A1' peak (ps)"] = t[ip]
    out["A1' rise 63% (ps)"] = t[np.nonzero(x >= (1 - np.exp(-1)) * pk)[0][0]]
    after = np.nonzero((t > t[ip]) & (x <= pk / np.e))[0]
    out["A1' 1/e decay after peak (ps)"] = (t[after[0]] - t[ip]) if len(after) else np.nan
    out["A1' peak excess"] = pk
    m10 = (t >= -0.1) & (t <= 10)
    def tr(qq):
        sel = [k for k in range(len(pts)) if np.allclose(qf[lab_q[k]], qq)]
        return np.array([np.sum((n[it][lab_q[sel]] + 0.5) * Wd[sel]) for it in range(len(t))])
    y = tr([1/3, 1/3]); y = y - y[0]
    try:
        p, _ = curve_fit(fK, t[m10], y[m10], p0=[y.max(), 0.2, 2, y[-1]], maxfev=20000)
        out["diffuse K rise (fit, ps)"], out["diffuse K decay (fit, ps)"] = p[1], p[2]
    except Exception:
        out["diffuse K rise (fit, ps)"] = out["diffuse K decay (fit, ps)"] = np.nan
    for name, qq in (("K/2", [1/6, 1/6]), ("M/2", [0.25, 0.0]), ("M", [0.5, 0.0])):
        y = tr(qq); y = y - y[0]
        try:
            p, _ = curve_fit(f1, t[m10], y[m10], p0=[y[m10].max(), 1, 0], bounds=([0, 0.01, -0.5], [np.inf, 1e3, 1]), maxfev=20000)
            out[f"diffuse {name} rise (fit, ps)"] = p[1]
        except Exception:
            out[f"diffuse {name} rise (fit, ps)"] = np.nan
    for tk in (0.5, 1.5, 5.0, 100.0):
        it = np.argmin(abs(t - tk))
        for ch in ("TO", "LO", "LA", "TA", "ZO", "ZA"):
            msk = label == ch
            out[f"E_{ch} at {tk:g} ps (meV/cell)"] = 1e3 * np.sum(((n[it] - n[0]) * hw)[msk]) / Nq
        out[f"Te at {tk:g} ps (K)"] = float(d["Te"][it]) if "Te" in d.files else np.nan
    out["total phonon excess at 100 ps (meV/cell)"] = 1e3 * np.sum((n[-1] - n[0]) * hw) / Nq
    out["absorbed energy (meV/cell)"] = 1e3 * float(d["E_abs_cell"]) if "E_abs_cell" in d.files else np.nan
    return out, t, n, d


if __name__ == "__main__":
    res = {}; series = {}
    for name, fn in RUNS.items():
        if os.path.exists(fn):
            res[name], t, n, d = metrics(fn)
            series[name] = (t, n[:, iK, 5] - n[0, iK, 5], d["Te"] if "Te" in d.files else None, n)
    keys = list(next(iter(res.values())).keys())
    lines = ["quantity".ljust(42) + "".join(k.rjust(25) for k in res)]
    for k in keys:
        lines.append(k.ljust(42) + "".join(f"{res[r][k]:25.3f}" for r in res))
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(HERE, "comparison.txt"), "w").write(txt + "\n")
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    cols = {"Level B (thermal e)": "#52514e", "Level B (non-thermal e)": "#eb6834", "unified": "#2a78d6",
            "unified, field x1.09": "#1baf7a", "unified, local-U e-e": "#6250d6"}
    for name, (t, x, Te, n) in series.items():
        ax[0].plot(t, x, color=cols[name], label=name)
        if Te is not None:
            ax[1].plot(t, Te, color=cols[name], label=name)
        Eb = [1e3 * np.sum(((n[i] - n[0]) * hw)[(label == "ZA") | (label == "TA") | (label == "LA")]) / Nq for i in range(len(t))]
        ax[2].plot(t, Eb, color=cols[name])
    ax[0].set_xlim(-0.2, 10); ax[0].set_xlabel("time (ps)"); ax[0].set_ylabel("A1'(K) excess population"); ax[0].legend(fontsize=7)
    ax[1].set_xscale("symlog", linthresh=1); ax[1].set_xlabel("time (ps)"); ax[1].set_ylabel("electron temperature (K)")
    ax[2].set_xscale("symlog", linthresh=1); ax[2].set_xlabel("time (ps)"); ax[2].set_ylabel("acoustic excess energy (meV/cell)")
    for a in ax:
        a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "graphite_unified.png"), dpi=140)
