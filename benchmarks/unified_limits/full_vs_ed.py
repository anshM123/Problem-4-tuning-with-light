"""(a) Full level against exact diagonalization: Hubbard-Holstein dimer with anharmonic
phonons, driven by a Peierls light pulse; e-e, e-ph and ph-ph interactions all active.

Model (hbar = J = 1): two sites, one up + one down electron, U = 0 and U = 2.  Each site carries an
optical (Holstein) mode a_i (w = 1, g n_i (a_i + a_i^+), g = 0.05 ... 0.2) and a soft mode b_i
(w = 0.5); the cubic term lam x_a x_b^2 (lam = 0.1) makes a_i decay resonantly into pairs of b_i,
the analogue of optical-phonon decay into acoustic phonons.  Light: Peierls pulse
A(t) = 0.6 cos(2 (t - 8)) exp(-(t - 8)^2 / 8) on the bond.

Initial state: both electrons on site 0 (doublon), phonon vacua (uncorrelated).
Compared with exact diagonalization (phonon cutoffs 10,10,13,13; change against 8,8,10,10 printed):
  full       unified solver, sic pair ladder + e-ph (incl. phonon-exchange G2) + ph-ph correlations
  ehrenfest  same e-e treatment, phonons at mean-field (Ehrenfest) level
  meanfield  TDHF + Ehrenfest
Usage: python full_vs_ed.py [tmax]    -> data/full_vs_ed.npz, printed errors
"""
import os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from lightmat.unified import Model, UnifiedSolver, peierls_chain
from ed_eph import EPhED

J, U = 1.0, 2.0
WA, WB, G, LAM = 1.0, 0.5, 0.2, 0.1
AL, WL, TL, SL = 0.6, 2.0, 8.0, 2.0


def Afun(t):
    return AL * np.cos(WL * (t - TL)) * np.exp(-(t - TL) ** 2 / (2 * SL ** 2))


def build(g=G, lam=LAM, U=U):
    L, Q = 2, 4                         # modes: a_0, a_1, b_0, b_1
    omega = np.array([WA, WA, WB, WB])
    gamma = np.zeros((Q, L, L))
    gamma[0, 0, 0] = gamma[1, 1, 1] = np.sqrt(2) * g
    kappa = np.zeros((Q, Q, Q))
    for i in range(2):                  # lam x_a x_b^2  ->  kappa_abb = 2 lam (all orderings)
        a, b = i, 2 + i
        for p in ((a, b, b), (b, a, b), (b, b, a)):
            kappa[p] = 2 * lam
    return L, Q, omega, gamma, kappa


def initial(kind):
    if kind == "doublon":
        r = np.diag([1.0, 0.0]).astype(complex)
        return r, r, [0], [0]
    v = np.array([1, 1]) / np.sqrt(2)
    r = np.outer(v, v).astype(complex)
    return r, r, None, None


def run_unified(kind, times, dt=0.01, g=G, lam=LAM, U=U, **kw):
    L, Q, omega, gamma, kappa = build(g, lam, U)
    m = Model(peierls_chain(L, J, Afun), U=U, omega=omega, gamma=gamma, kappa=kappa)
    ru, rd, _, _ = initial(kind)
    s = UnifiedSolver(m, ru, rd, **kw)
    E0 = s.energy()
    out = []
    nsub = int(round((times[1] - times[0]) / dt))
    for k, T in enumerate(times):
        if k > 0:
            for _ in range(nsub):
                s.step(dt)
        nu, nd = s.densities()
        out.append(dict(d=s.double_occ(), n=nu + nd, N=s.phonon_number(), E=s.energy() - s.W - E0))
    return out


def run_ed(kind, times, ncut=(10, 10, 14, 14), g=G, lam=LAM, U=U):
    L, Q, omega, gamma, kappa = build(g, lam, U)
    ed = EPhED(L, 1, 1, [(0, 1)], J, U, omega, gamma, kappa, list(ncut), Afun=Afun)
    if kind == "doublon":
        psi = ed.initial([0], [0])
    else:   # bonding orbital for both spins: (|0> + |1>)_up (|0> + |1>)_dn / 2 x vacuum
        pe = sum(ed.el.product_state([i], [j]) for i in range(2) for j in range(2)) / 2
        v = np.zeros(ed.dph); v[0] = 1
        psi = np.kron(pe, v).astype(complex)
    return ed.run(psi, times, nsub=8)


def table(res):
    D = np.array([np.sum(r["d"]) for r in res])
    Na = np.array([r["N"][0] + r["N"][1] for r in res])
    Nb = np.array([r["N"][2] + r["N"][3] for r in res])
    return D, WA * Na, WB * Nb


VARIANTS = (("full", {}),                                                     # all correlations
            ("ehrenfest", dict(eph="meanfield", phph="meanfield")),            # same e-e, phonons mean field
            ("meanfield", dict(ee="hf", eph="meanfield", phph="meanfield")))   # TDHF + Ehrenfest


def compare(kind, U, g, lam, times, ncut=(8, 8, 10, 10), ncut2=(10, 10, 13, 13)):
    ex = np.array(table(run_ed(kind, times, ncut=ncut2, g=g, lam=lam, U=U)))
    ex1 = np.array(table(run_ed(kind, times, ncut=ncut, g=g, lam=lam, U=U)))
    res = {"ed": ex, "cut": np.abs(ex - ex1).max(1)}
    for name, kw in VARIANTS:
        kw = dict(kw)
        if U == 0 and kw.get("ee") is None:
            kw["ee"] = "hf"
        out = run_unified(kind, times, g=g, lam=lam, U=U, **kw)
        res[name] = np.array(table(out))
        res[name + "_E"] = np.array([r["E"] for r in out])
    return res


if __name__ == "__main__":
    tmax = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
    times = np.arange(0, tmax + 1e-9, 0.25)
    store = {"t": times}
    for U in (0.0, 2.0):
        for g in (0.05, 0.1, 0.2):
            t0 = time.time()
            r = compare("doublon", U, g, LAM, times)
            for k, v in r.items():
                store[f"U{U:g}_g{g:g}_{k}"] = v
            line = f"U={U:g} g={g:g} lam={LAM:g} (doublon, pulse): ED cutoff change {r['cut'].max():.1e};"
            for name, _ in VARIANTS:
                e = np.abs(r[name] - r["ed"]).max(1)
                line += f" {name}: max err D {e[0]:.3f} E_opt {e[1]:.3f} E_soft {e[2]:.3f};"
            line += f" |E-W-E0| max {np.abs(r['full_E']).max():.1e}; scale: max D {r['ed'][0].max():.2f}, E_opt {r['ed'][1].max():.3f}, E_soft {r['ed'][2].max():.3f} [{time.time()-t0:.0f}s]"
            print(line, flush=True)
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    np.savez(os.path.join(HERE, "data", "full_vs_ed.npz"), **store)
