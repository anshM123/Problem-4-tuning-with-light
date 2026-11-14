"""Figures for the graphite simulation: A1'(K) population vs time and branch-resolved
excess phonon populations over the Brillouin zone at 0.5, 1.5, 5 and 100 ps."""
import os
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "..", "figures"); os.makedirs(FIG, exist_ok=True)
d = np.load(os.path.join(HERE, "dynamics_out.npz"))
t, n, Te, f, q = d["t"], d["n"], d["Te"], d["freq"], d["qfrac"]
hw = f * 4.135667e-3; kB = 8.617333e-5
Nq = len(q); N1 = int(round(np.sqrt(Nq)))

# branch character from eigenvectors: out-of-plane (Z), longitudinal (L), transverse (T);
# acoustic vs optical from the relative phase of the two atoms
a_ = 2.464
cell = [[a_, 0, 0], [-a_ / 2, a_ * np.sqrt(3) / 2, 0], [0, 0, 15.0]]
ph = Phonopy(PhonopyAtoms(symbols=["C", "C"], cell=cell, scaled_positions=[[1/3, 2/3, .5], [2/3, 1/3, .5]]),
             supercell_matrix=[8, 8, 1])
ph.force_constants = h5py.File(os.path.join(HERE, "fc2.hdf5"))["force_constants"][:]
ph.run_qpoints(np.c_[q, np.zeros(Nq)], with_eigenvectors=True)
ev = ph.qpoints.eigenvectors
B = 2 * np.pi * np.linalg.inv(np.array(cell)[:2, :2]).T
qc = q @ B
label = np.empty((Nq, 6), object)
for i in range(Nq):
    qhat = qc[i] / (np.linalg.norm(qc[i]) + 1e-12)
    zs = [j for j in range(6) if np.sum(abs(ev[i][:, j].reshape(2, 3)[:, 2]) ** 2) > 0.5]
    ins = [j for j in range(6) if j not in zs]
    zs = sorted(zs, key=lambda j: f[i, j]); ins = sorted(ins, key=lambda j: f[i, j])
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
iK = [i for i, x in enumerate(q) if np.allclose(x, [1/3, 1/3])][0]
Teff = lambda nn, w: w / (kB * np.log1p(1 / np.maximum(nn, 1e-12)))

# ---- Fig: A1'(K) population and electron temperature ----------------------------
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
nA = n[:, iK, 5]
ax[0].plot(t, Teff(nA, hw[iK, 5]), label="A1' (K), effective T (this work)")
stern_t = [1.0, 10.0]; stern_T = [1300, 450]
ax[0].plot(stern_t, stern_T, "ks", label="K-TO temperature, Stern et al. 2018")
ax[0].set_xscale("symlog", linthresh=1); ax[0].set_xlabel("time (ps)"); ax[0].set_ylabel("effective temperature (K)")
ax[0].legend(fontsize=7)
frac = (nA - nA[0]) / (nA.max() - nA[0])
ax[1].plot(t, frac, label="A1'(K) excess population (normalized)")
tt = np.linspace(0, 10, 200)
ax[1].plot(tt, np.where(tt > 0, (1 - np.exp(-tt / 0.28)) * np.exp(-tt / 1.7) / 0.70, 0), "k--",
           label="exp. rise 0.28 ps / decay 1.7 ps")
ax[1].set_xlim(-0.2, 10); ax[1].set_xlabel("time (ps)"); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "graphite_A1K.png"), dpi=150)

# ---- Fig: excess populations over the BZ per branch character ---------------------
times = [0.5, 1.5, 5.0, 100.0]
chars = ["TO", "LO", "LA", "TA", "ZO", "ZA"]
fig, axs = plt.subplots(len(chars), len(times), figsize=(11, 15))
# extended zone points for display
ext = []
for s1 in (-1, 0, 1):
    for s2 in (-1, 0, 1):
        ext.append((q + np.array([s1, s2])) @ B)
ext = np.concatenate(ext)
for c, ch in enumerate(chars):
    for k, tk in enumerate(times):
        i = np.argmin(abs(t - tk))
        dn = np.zeros(Nq)
        for qq in range(Nq):
            js = [j for j in range(6) if label[qq, j] == ch]
            dn[qq] = sum(n[i, qq, j] - n[0, qq, j] for j in js)
        ax = axs[c, k]
        sc = ax.scatter(ext[:, 0], ext[:, 1], c=np.tile(dn, 9), s=6, cmap="inferno",
                        vmin=0, vmax=max(dn.max(), 1e-6))
        ax.set_xlim(-2.0, 2.0); ax.set_ylim(-2.0, 2.0); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        if c == 0:
            ax.set_title(f"{tk:g} ps")
        if k == 0:
            ax.set_ylabel(ch)
        fig.colorbar(sc, ax=ax, fraction=0.046)
fig.suptitle("Excess phonon population n(q,t) - n(q,0) per branch (graphene-layer BZ)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "graphite_bz_maps.png"), dpi=120)

# ---- energy per branch character vs time -----------------------------------------
fig, ax = plt.subplots(figsize=(6, 3.6))
for ch in chars:
    Eb = []
    for i in range(len(t)):
        Eb.append(sum((n[i, qq, j] - n[0, qq, j]) * hw[qq, j] for qq in range(Nq) for j in range(6)
                      if label[qq, j] == ch) / Nq * 1e3)
    ax.plot(t, Eb, label=ch)
ax.set_xscale("symlog", linthresh=1); ax.set_xlabel("time (ps)"); ax.set_ylabel("excess energy (meV / cell)")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(FIG, "graphite_branch_energy.png"), dpi=150)
# text summary
with open(os.path.join(HERE, "summary.txt"), "w") as fh:
    for tk in times:
        i = np.argmin(abs(t - tk))
        parts = []
        for ch in chars:
            E = sum((n[i, qq, j] - n[0, qq, j]) * hw[qq, j] for qq in range(Nq) for j in range(6) if label[qq, j] == ch) / Nq * 1e3
            parts.append(f"{ch} {E:.2f}")
        fh.write(f"t={tk:g} ps  Te={Te[i]:.0f} K  excess energy (meV/cell): " + ", ".join(parts) + "\n")
print(open(os.path.join(HERE, "summary.txt")).read())
