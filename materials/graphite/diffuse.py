"""One-phonon ultrafast electron diffuse scattering (UEDS) of photoexcited graphite from the
simulated mode- and wavevector-resolved phonon populations.

  I1(Q, t) ~ sum_lambda (n_lambda(q,t) + 1/2) / omega_lambda(q) |F_lambda(Q)|^2,
  F_lambda(Q) = sum_s f_C(Q) / sqrt(M) (Q . e_{s,lambda}(q)) exp(i G . r_s),   Q = q + G,
with the Doyle-Turner electron scattering factor of carbon and the Debye-Waller factor
neglected (it changes slowly and only rescales maps).  Differential maps
dI(Q,t) = I1(Q,t) - I1(Q,0) are computed in the in-plane reciprocal plane up to |Q| = 7 A^-1
(the range of the experiment), together with time traces at K, M and near-Gamma points and
their rise/decay times, which are compared with the published UEDS time constants
(Stern et al., PRB 97, 165416 (2018); Rene de Cotret et al., PRB 100, 214115 (2019)).
"""
import os, sys
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "..", "figures")
d = np.load(os.path.join(HERE, sys.argv[1] if len(sys.argv) > 1 else "dynamics_out.npz"))
t, n, freq, qf = d["t"], d["n"].astype(float), d["freq"], d["qfrac"]
Nq = len(qf); N1 = int(round(np.sqrt(Nq)))
a_ = 2.464
cell = [[a_, 0, 0], [-a_ / 2, a_ * np.sqrt(3) / 2, 0], [0, 0, 15.0]]
ph = Phonopy(PhonopyAtoms(symbols=["C", "C"], cell=cell, scaled_positions=[[1/3, 2/3, .5], [2/3, 1/3, .5]]),
             supercell_matrix=[8, 8, 1])
ph.force_constants = h5py.File(os.path.join(HERE, "fc2.hdf5"))["force_constants"][:]
ph.run_qpoints(np.c_[qf, np.zeros(Nq)], with_eigenvectors=True)
ev = ph.qpoints.eigenvectors                                   # (Nq, 6, 6)
B = 2 * np.pi * np.linalg.inv(np.array(cell)[:2, :2]).T        # rows b1, b2 (1/A)
rfrac = np.array([[1/3, 2/3], [2/3, 1/3]])
# Doyle-Turner electron scattering factor of carbon, s = Q / (4 pi)
DT_a = [0.7307, 1.1951, 0.4563, 0.1247]; DT_b = [36.9951, 11.2966, 2.8139, 0.3456]
def fC(Q):
    s2 = (Q / (4 * np.pi)) ** 2
    return sum(a * np.exp(-b * s2) for a, b in zip(DT_a, DT_b))

# reciprocal-plane sampling: all grid q plus G vectors with |Q| < 7 A^-1
pts, lab_q, lab_G = [], [], []
for g1 in range(-4, 5):
    for g2 in range(-4, 5):
        G = np.array([g1, g2])
        for iq in range(Nq):
            Qf = qf[iq] + G
            Q = Qf @ B
            if 0.3 < np.linalg.norm(Q) < 7.0:
                pts.append(Q); lab_q.append(iq); lab_G.append(G)
pts = np.array(pts); lab_q = np.array(lab_q); lab_G = np.array(lab_G)
absQ = np.linalg.norm(pts, axis=1)
# |F|^2 / omega for every point and branch
W = np.zeros((len(pts), 6))
for k in range(len(pts)):
    iq = lab_q[k]; Q = pts[k]; G = lab_G[k]
    for j in range(6):
        if freq[iq, j] < 0.1:
            continue
        e = ev[iq][:, j].reshape(2, 3)[:, :2]
        F = sum(fC(absQ[k]) * (Q @ e[s]) * np.exp(-2j * np.pi * (G @ rfrac[s])) for s in range(2))
        W[k, j] = abs(F) ** 2 / freq[iq, j]
I = lambda it: np.sum((n[it][lab_q] + 0.5) * W, axis=1)
I0 = I(0)
times = [0.5, 1.5, 5.0, 100.0]
fig, ax = plt.subplots(1, 4, figsize=(17, 4.3))
for k, tk in enumerate(times):
    it = np.argmin(abs(t - tk))
    dI = I(it) - I0
    rel = dI / np.maximum(I0, 1e-30)
    sc = ax[k].scatter(pts[:, 0], pts[:, 1], c=rel, s=3, cmap="inferno", vmin=0, vmax=np.percentile(rel, 99))
    ax[k].set_aspect("equal"); ax[k].set_title(f"dI/I0 at {tk:g} ps"); ax[k].set_xticks([]); ax[k].set_yticks([])
    fig.colorbar(sc, ax=ax[k], fraction=0.046)
fig.suptitle("Simulated one-phonon diffuse scattering of graphite, differential maps (in-plane, |Q| < 7 1/A)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "graphite_diffuse_maps.png"), dpi=120)

# time traces at reduced wavevectors K, M and near Gamma (averaged over all zones)
def trace(qtarget):
    sel = [k for k in range(len(pts)) if np.allclose(qf[lab_q[k]], qtarget)]
    return np.array([np.sum((n[it][lab_q[sel]] + 0.5) * W[sel]) for it in range(len(t))])
def fit_rise(tr):
    y = (tr - tr[0]) / (tr.max() - tr[0]); k63 = np.where(y > 1 - np.exp(-1))[0][0]
    return t[k63], t[np.argmax(tr)]
out = {}
for name, qq in (("K", [1/3, 1/3]), ("M", [0.5, 0.0]), ("Gamma+1", [1 / N1, 0.0])):
    tr = trace(qq); out[name] = tr
    r, p = fit_rise(tr)
    print(f"{name:8s}: rise (63%) {r:.2f} ps, maximum at {p:.1f} ps")
fig, ax = plt.subplots(figsize=(6, 3.6))
for name, tr in out.items():
    ax.plot(t, (tr - tr[0]) / (tr.max() - tr[0]), label=name)
ax.set_xscale("symlog", linthresh=1); ax.set_xlabel("time (ps)"); ax.set_ylabel("dI (normalized)"); ax.legend()
ax.set_title("diffuse intensity at reduced wavevectors")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "graphite_diffuse_traces.png"), dpi=140)
print("measured (Stern 2018): K rise 0.28 ps (A1'), K-LO 0.73 ps, M rise 2.1 ps, acoustic thermalization 34 ps")

# Time constants extracted the way the experiment reports them (fits of the diffuse intensity,
# not 63 % times of a population):
#  K: rise x decay to a plateau, dI = A (1 - e^{-t/tr}) e^{-t/td} + C (1 - e^{-t/td}), 0..10 ps
#  mid-zone (K/2, M/2) and M: single-exponential rise fitted to the early dynamics, 0..10 ps
from scipy.optimize import curve_fit
def tr_at(qq):
    sel = [k for k in range(len(pts)) if np.allclose(qf[lab_q[k]], qq)]
    return np.array([np.sum((n[it][lab_q[sel]] + 0.5) * W[sel]) for it in range(len(t))])
m10 = (t >= -0.1) & (t <= 10)
y = tr_at([1/3, 1/3]); y = y - y[0]
fK = lambda x, A, tr_, td, C: np.where(x > 0, A * (1 - np.exp(-x / tr_)) * np.exp(-x / td) + C * (1 - np.exp(-x / td)), 0)
pK, _ = curve_fit(fK, t[m10], y[m10], p0=[y.max(), 0.2, 2, y[-1]], maxfev=20000)
fits = {"K_rise": pK[1], "K_decay": pK[2]}
f1 = lambda x, A, tau, t0: np.where(x > t0, A * (1 - np.exp(-(x - t0) / tau)), 0)
for name, qq in (("K/2", [1/6, 1/6]), ("M/2", [0.25, 0.0]), ("M", [0.5, 0.0])):
    y = tr_at(qq); y = y - y[0]
    p, _ = curve_fit(f1, t[m10], y[m10], p0=[y[m10].max(), 1, 0], bounds=([0, 0.01, -0.5], [np.inf, 1e3, 1]), maxfev=20000)
    fits[name + "_rise"] = p[1]
print("experiment-style fits: K rise %.2f ps, K decay %.2f ps (measured 0.28, 1.7); K/2 rise %.2f, M/2 rise %.2f ps "
      "(measured 1.0-1.5); M rise %.1f ps (measured 2.1)" % (fits["K_rise"], fits["K_decay"], fits["K/2_rise"],
                                                          fits["M/2_rise"], fits["M_rise"]))
np.savez(os.path.join(HERE, "diffuse_fits.npz"), **fits)
