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
    sc = os.path.join(D, f"fast_cdw_U{U}_sic.npz")
    h = os.path.join(D, f"fast_cdw_U{U}_hf.npz")
    ax[0, c].plot(e["t"], e["I"], "k-", lw=1.2, label="exact (ED)")
    ax[0, c].plot(f["t"], f["I"], "r-", lw=1, label="T-matrix pp ladder (factorized)")
    if os.path.exists(sc):
        ss = np.load(sc); ax[0, c].plot(ss["t"], ss["I"], "-", color="tab:blue", lw=1,
                                        label="self-interaction-free pp ladder (this work)")
        ax[1, c].plot(ss["t"], ss["dcorr"].sum(1), "-", color="tab:blue", lw=1)
    if os.path.exists(h):
        hh = np.load(h); ax[0, c].plot(hh["t"], hh["I"], "g:", lw=1, label="TDHF")
    ax[0, c].set_title(f"U/J = {U}"); ax[0, c].set_ylabel("density imbalance I")
    ax[1, c].plot(e["t"], e["dcorr"].sum(1), "k-", lw=1.2)
    ax[1, c].plot(f["t"], f["dcorr"].sum(1), "r-", lw=1)
    ax[1, c].set_ylabel(r"$\sum_i d^{corr}_i$"); ax[1, c].set_xlabel("t J / hbar")
ax[0, 0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "hubbard_cdw.png"), dpi=150)

# U sweep: the quantities of Fig. 3 of the reference plus the melting time t(I < 0.5)
Us = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 15, 16, 18, 20]
def metrics(t, I, dc):
    k = np.where(I < 0)[0]; z = t[k[0]] if len(k) else np.nan
    k = np.where(I < 0.5)[0]; m = t[k[0]] if len(k) else np.nan
    w = (t >= 20) & (t <= 50)
    return [I[w].mean(), z, m, dc.sum(1)[w].mean()]
rows = []
for U in Us:
    r = [U]
    for fn in (f"ed_cdw_U{U}.npz", f"fast_cdw_U{U}_factorized.npz", f"fast_cdw_U{U}_sic.npz"):
        p = os.path.join(D, fn)
        if os.path.exists(p):
            z = np.load(p); r += metrics(z["t"], z["I"], z["dcorr"])
        else:
            r += [np.nan] * 4
    rows.append(r)
S = np.array(rows)
np.savetxt(os.path.join(D, "cdw_sweep_metrics.txt"), S, fmt="%.4f",
           header="U | ED: <I>, t_zero, t_melt(I<0.5), <dcorr> | T-matrix: same | self-interaction-free: same")
fig, ax = plt.subplots(1, 4, figsize=(15, 3.4))
lab = ["exact", "T-matrix pp ladder", "self-interaction-free pp ladder"]
sty = ["ko-", "rs-", "^-"]
ylab = [r"$\langle I\rangle_{t\in[20,50]}$", "first zero crossing of I (hbar/J)",
        "melting time, I < 0.5 (hbar/J)", r"$\langle\sum d^{corr}\rangle_{t\in[20,50]}$"]
for k in range(4):
    for m in range(3):
        kw = dict(color="tab:blue") if m == 2 else {}
        ax[k].plot(S[:, 0], S[:, 1 + 4 * m + k], sty[m], ms=4, label=lab[m], **kw)
    ax[k].set_xlabel("U/J"); ax[k].set_ylabel(ylab[k])
ax[1].set_ylim(0, 30)
ax[0].legend(fontsize=7)
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
# Top: like-for-like on the open 14-site chain (exact diagonalization vs fast methods, all
# with Lorentzian eta = 0.1 J).  Bottom: long chains and the 14-site ring vs the exact
# Bethe-ansatz gap.  Numbers: data/chain_sigma_summary.txt (chain_sigma_compare.py).
def _sig(fn):
    p = os.path.join(D, fn)
    if not os.path.exists(p):
        return None
    d = np.load(p)
    return d["w"], np.real(d["sigma"])
fig, ax = plt.subplots(2, 3, figsize=(14, 6.6))
for k, U in enumerate([3, 6, 12]):
    pe = os.path.join(D, f"chain_sigma_ed_obc_L14_U{U}.npz")
    gap = float(np.load(pe)["gap"]) if os.path.exists(pe) else None
    ymax = [0.0, 0.0]
    for fn, lab, sty in ((f"chain_sigma_ed_obc_L14_U{U}.npz", "exact (ED), open 14 sites", "k-"),
                         (f"chain_sigma_dh_sw_L14_U{U}.npz", "doublon-holon pairs, HF + spin waves (this work)", "b-"),
                         (f"chain_sigma_dh_L14_U{U}.npz", "doublon-holon pairs, HF ground state", "c--"),
                         (f"chain_sigma_U{U}_gwx_c0_L14.npz", "TDHF/RPA", "r--"),
                         (f"chain_sigma_U{U}_gwx_c1_L14.npz", "ph fluctuations + collisions", "g:")):
        z = _sig(fn)
        if z is not None:
            ax[0, k].plot(z[0], z[1], sty, lw=1, label=lab)
            if "ed_obc" in fn or "dh_sw" in fn:
                ymax[0] = max(ymax[0], z[1].max())
    for fn, lab, sty in ((f"chain_sigma_ed_pbc_L14_U{U}.npz", "exact (ED), ring 14 sites", "k-"),
                         (f"chain_sigma_dh_sw_ring_L128_U{U}.npz", "doublon-holon pairs, HF + spin waves, ring 128", "b-"),
                         (f"chain_sigma_dh_ring_L128_U{U}.npz", "doublon-holon pairs, HF, ring 128", "c--"),
                         (f"chain_sigma_U{U}_gwx_c0.npz", "TDHF/RPA, 64 sites", "r--")):
        z = _sig(fn)
        if z is not None:
            ax[1, k].plot(z[0], z[1], sty, lw=1, label=lab)
            ymax[1] = max(ymax[1], z[1].max())
    for r in range(2):
        if gap is not None:
            ax[r, k].axvline(gap, color="orange", ls="--", lw=1, label="Bethe-ansatz gap (L = inf)")
        ax[r, k].set_xlim(0, 1.5 * U + 6); ax[r, k].set_ylim(0, 1.25 * ymax[r] if ymax[r] > 0 else None)
        ax[r, k].set_xlabel(r"$\hbar\omega/J$")
    ax[0, k].set_title(f"U/J = {U}")
ax[0, 0].set_ylabel(r"Re $\sigma(\omega)$, open 14-site chain"); ax[1, 0].set_ylabel(r"Re $\sigma(\omega)$, long chains")
ax[0, 0].legend(fontsize=7); ax[1, 0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "hubbard_chain_sigma.png"), dpi=140)
print("ok")
