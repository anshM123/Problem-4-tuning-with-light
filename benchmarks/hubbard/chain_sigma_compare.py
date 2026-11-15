"""Thresholds, peaks and line-shape distances of the chain optical conductivity: fast methods
vs exact diagonalization on the same open chain, and vs the Bethe-ansatz gap.

Definitions (identical for every method and for ED; all spectra Re sigma(w) with the same
Lorentzian broadening eta = 0.1 J on the grid w = linspace(0.02, 1.5 U + 6, 700)):
  threshold  w_th : first frequency at which Re sigma reaches 10 % of its maximum;
  peak       w_pk : position of the maximum;
  distance   D    : int |s_m/W_m - s_e/W_e| dw with W = int s dw  (0 = identical line shape,
                    2 = no overlap);
  weight     W_m/W_e (f-sum ratio over the window).
The thermodynamic-limit threshold of the doublon-holon closure is its bulk band edge
U - 4J - 4 e_b + K_b (Hartree-Fock bond energy and kinetic energy of the 128-site ring).
Usage: python chain_sigma_compare.py      (writes data/chain_sigma_summary.txt)
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chain_sigma_ed import bethe_gap, sigma_from_poles

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
Us = [3, 6, 12]
LREF = 14


def load(name):
    p = os.path.join(D, name)
    if not os.path.exists(p):
        return None
    d = np.load(p)
    return d["w"], np.real(d["sigma"])


def metrics(w, s):
    return w[np.argmax(s >= 0.1 * s.max())], w[np.argmax(s)]


def respec(name, eta, w):
    """Re-broadened spectrum (arbitrary scale) from the stored current j(t) or ED poles."""
    d = np.load(os.path.join(D, name))
    if "wn" in d:
        return sigma_from_poles(w, d["wn"], d["an"], int(d["L"]), eta)
    t, j = d["t"], d["j"]; dt = t[1] - t[0]
    return np.array([-np.sum(j * np.exp((1j * x - eta) * t)) * dt for x in w]).real


def distance(w, s, we, se):
    s = np.interp(we, w, s)
    dw = we[1] - we[0]
    Wm, We = s.sum() * dw, se.sum() * dw
    return np.abs(s / Wm - se / We).sum() * dw, Wm / We


def table():
    rows = []
    for U in Us:
        ed = load(f"chain_sigma_ed_obc_L{LREF}_U{U}.npz")
        pbc = load(f"chain_sigma_ed_pbc_L{LREF}_U{U}.npz")
        gap = bethe_gap(U)
        r = {"U": U, "bethe": gap}
        if ed:
            r["ed"] = metrics(*ed)
        if pbc:
            r["ed_pbc"] = metrics(*pbc)
        for key, fn in (("rpa", f"chain_sigma_U{U}_gwx_c0_L{LREF}.npz"),
                        ("corr", f"chain_sigma_U{U}_gwx_c1_L{LREF}.npz"),
                        ("dh1", f"chain_sigma_dh_o1_L{LREF}_U{U}.npz"),
                        ("dh", f"chain_sigma_dh_L{LREF}_U{U}.npz"),
                        ("dhsw", f"chain_sigma_dh_sw_L{LREF}_U{U}.npz")):
            m = load(fn)
            if m is not None and ed is not None:
                we = ed[0]
                env = distance(we, respec(fn, 0.5, we), we, respec(f"chain_sigma_ed_obc_L{LREF}_U{U}.npz", 0.5, we))[0]
                r[key] = metrics(*m) + distance(*m, *ed) + (env,)
        for key, fn in (("rpa64", f"chain_sigma_U{U}_gwx_c0.npz"), ("corr16", f"chain_sigma_U{U}_L16_corr.npz"),
                        ("dh128", f"chain_sigma_dh_ring_L128_U{U}.npz"),
                        ("dhsw128", f"chain_sigma_dh_sw_ring_L128_U{U}.npz")):
            m = load(fn)
            if m is not None:
                r[key] = metrics(*m)
        p = os.path.join(D, f"chain_sigma_dh_ring_L128_U{U}.npz")
        if os.path.exists(p):
            r["dh_edge"] = float(np.load(p)["edge"])
        p = os.path.join(D, f"chain_sigma_dh_sw_ring_L128_U{U}.npz")
        if os.path.exists(p):
            r["dhsw_edge"] = float(np.load(p)["edge"])
        rows.append(r)
    return rows


if __name__ == "__main__":
    rows = table()
    f = lambda x: "  -  " if x is None else f"{x:5.2f}"
    lines = ["# Optical conductivity of the half-filled Hubbard chain, eta = 0.1 J.",
             f"# Finite-chain comparisons on the same open {LREF}-site chain as the exact diagonalization.",
             "# th = 10%-of-maximum threshold, pk = peak, D = normalized line-shape distance (0..2) at eta = 0.1,",
             "# D(eta=0.5) = same distance after re-broadening both spectra to eta = 0.5 J (envelope), W = weight ratio",
             ""]
    for r in rows:
        lines.append(f"U/J = {r['U']}:  Bethe-ansatz gap {r['bethe']:.3f}")
        if "ed" in r:
            lines.append(f"  ED open L={LREF}:      th {r['ed'][0]:.2f}  pk {r['ed'][1]:.2f}")
        if "ed_pbc" in r:
            lines.append(f"  ED ring L={LREF}:      th {r['ed_pbc'][0]:.2f}  pk {r['ed_pbc'][1]:.2f}")
        for key, lab in (("rpa", "TDHF/RPA (gwx)"), ("corr", "ph fluct + collisions"),
                         ("dh1", "d-h pairs, leading order"), ("dh", "d-h pairs, HF"),
                         ("dhsw", "d-h, HF+spin waves (main)")):
            if key in r:
                th, pk, dist, W, env = r[key]
                lines.append(f"  {lab:25s} L={LREF}: th {th:.2f}  pk {pk:.2f}  D {dist:.2f}  "
                             f"D(eta=0.5) {env:.2f}  W {W:.2f}")
        for key, lab in (("rpa64", "TDHF/RPA 64 sites"), ("corr16", "ph fluct + coll. 16 sites"), ("dh128", "d-h pairs, ring 128 sites"),
                         ("dhsw128", "d-h + spin waves, ring 128")):
            if key in r:
                lines.append(f"  {lab:26s}: th {r[key][0]:.2f}  pk {r[key][1]:.2f}")
        if "dh_edge" in r:
            lines.append(f"  doublon-holon bulk band edge (L -> inf): {r['dh_edge']:.3f}  "
                         f"({100 * (r['dh_edge'] / r['bethe'] - 1):+.1f} % vs Bethe ansatz)")
        if "dhsw_edge" in r:
            lines.append(f"  d-h + spin waves bulk band edge (L -> inf): {r['dhsw_edge']:.3f}  "
                         f"({100 * (r['dhsw_edge'] / r['bethe'] - 1):+.1f} % vs Bethe ansatz)")
        lines.append("")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(D, "chain_sigma_summary.txt"), "w", encoding="utf8").write(txt + "\n")
