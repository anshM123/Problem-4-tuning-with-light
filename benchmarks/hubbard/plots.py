"""Hubbard benchmark figures: CDW melting (exact vs fast) and the U sweep; ladder spin
correlations."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "data")
FIG = os.path.join(HERE, "..", "..", "figures"); os.makedirs(FIG, exist_ok=True)

fig, ax = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
for c, U in enumerate([5, 15]):
    e = np.load(os.path.join(D, f"ed_cdw_U{U}.npz")); f = np.load(os.path.join(D, f"fast_cdw_U{U}_factorized.npz"))
    h = os.path.join(D, f"fast_cdw_U{U}_hf.npz")
    ax[0, c].plot(e["t"], e["I"], "k-", lw=1.2, label="exact (ED)")
    ax[0, c].plot(f["t"], f["I"], "r-", lw=1, label="factorized pp ladder (this work)")
    if os.path.exists(h):
        hh = np.load(h); ax[0, c].plot(hh["t"], hh["I"], "b:", lw=1, label="TDHF")
    ax[0, c].set_title(f"U/J = {U}"); ax[0, c].set_ylabel("density imbalance I")
    ax[1, c].plot(e["t"], e["dcorr"].sum(1), "k-", lw=1.2)
    ax[1, c].plot(f["t"], f["dcorr"].sum(1), "r-", lw=1)
    ax[1, c].set_ylabel(r"$\sum_i d^{corr}_i$"); ax[1, c].set_xlabel("t J / hbar")
ax[0, 0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "hubbard_cdw.png"), dpi=150)

S = np.loadtxt(os.path.join(D, "cdw_sweep_summary.txt"))
fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
ax[0].plot(S[:, 0], S[:, 1], "ko-", label="exact"); ax[0].plot(S[:, 0], S[:, 2], "rs-", label="fast")
ax[0].set_xlabel("U/J"); ax[0].set_ylabel(r"$\langle I\rangle_{t\in[20,50]}$"); ax[0].legend()
ax[1].plot(S[:, 0], S[:, 3], "ko-"); ax[1].plot(S[:, 0], S[:, 4], "rs-")
ax[1].set_xlabel("U/J"); ax[1].set_ylabel("first zero crossing of I (hbar/J)")
ax[2].plot(S[:, 0], S[:, 5], "ko-"); ax[2].plot(S[:, 0], S[:, 6], "rs-")
ax[2].set_xlabel("U/J"); ax[2].set_ylabel(r"$\langle\sum d^{corr}\rangle_{t\in[20,50]}$")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "hubbard_cdw_sweep.png"), dpi=150)

ed = np.loadtxt(os.path.join(D, "ladder_gs_ed.txt"))
fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
others = [("ladder_gs_ppladder.txt", "pp ladder (factorized, this work)", "bs-"),
          ("ladder_gs_fast_gwx.txt", "ph fluctuations (RPA kernel)", "r^-")]
for k, name in enumerate(["C_loc", "C_leg", "C_rung"]):
    ax[k].plot(ed[:, 0], ed[:, k + 1], "ko-", label="exact (ED)")
    for fn, lab, sty in others:
        p = os.path.join(D, fn)
        if os.path.exists(p):
            d = np.atleast_2d(np.loadtxt(p)); ax[k].plot(d[:, 0], d[:, k + 1], sty, label=lab)
    ax[k].set_xlabel("U/J"); ax[k].set_title(name)
ax[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "hubbard_ladder.png"), dpi=150)

# ---- laser-driven ladder ---------------------------------------------------------
pe = os.path.join(D, "ladder_laser_ed.npz")
if os.path.exists(pe):
    fig, ax = plt.subplots(2, 3, figsize=(13, 6))
    sets = [("ladder_laser_ed.npz", "exact (ED)", "k-"), ("ladder_laser_fast_gwx.npz", "ph fluctuations", "r-"),
            ("ladder_laser_pp.npz", "pp ladder", "b-")]
    w = np.linspace(2.5, 5.0, 300)
    for fn, lab, sty in sets:
        d = np.load(os.path.join(D, fn)); t = d["t"]; dC = d["C"] - d["C"][0]
        m = (t > 15) & (t <= 100)
        for k, name in enumerate(["C_loc", "C_leg", "C_rung"]):
            ax[0, k].plot(t, dC[:, k], sty, lw=1, label=lab)
            x = dC[m, k] - dC[m, k].mean()
            win = np.hanning(m.sum())
            spec = np.abs(np.sum((x * win)[:, None] * np.exp(1j * w[None, :] * t[m][:, None]), 0))
            ax[1, k].plot(w, spec / spec.max(), sty, lw=1, label=lab)
            ax[0, k].set_title(r"$\Delta$" + name); ax[1, k].set_xlabel(r"$\hbar\omega/J$")
    ax[0, 0].legend(fontsize=7); ax[0, 0].set_xlabel("t J/hbar"); ax[1, 0].set_ylabel("|FT| (norm.)")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "hubbard_ladder_laser.png"), dpi=140)

# ---- optical conductivity of the chain -----------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(13, 3.4))
for k, U in enumerate([3, 6, 12]):
    p = os.path.join(D, f"chain_sigma_U{U}_gwx_c0.npz")
    if not os.path.exists(p):
        continue
    d = np.load(p)
    re = d["sigma"].real
    ax[k].plot(d["w"], re / re.max(), "r-", label="TDHF/RPA linear response (64 sites)")
    ax[k].axvline(float(d["gap"]), color="orange", ls="--", label="exact Bethe-ansatz gap")
    ax[k].set_title(f"U/J = {U}"); ax[k].set_xlabel(r"$\hbar\omega/J$"); ax[k].set_xlim(0, 1.5 * U + 6)
ax[0].legend(fontsize=7); ax[0].set_ylabel(r"Re $\sigma(\omega)$ (norm.)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "hubbard_chain_sigma.png"), dpi=140)
print("ok")
