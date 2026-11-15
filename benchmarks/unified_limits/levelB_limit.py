"""(c) Level B limit: the kinetic equations of materials/graphite as the Markov + random-phase
limit of the unified solver.

Model (levelB_model.py): L electron levels (spin degenerate) in a band of width 4, one optical
mode (w0 = 1) with random e-ph matrix elements, M acoustic modes (0.2..0.8) into which the optical
mode decays by cubic anharmonicity (kappa_{0bc}).

1. Collision operators.  At random (f, N) and at thermal electron states, the solver's
   kinetic limit (UnifiedSolver.kinetic_rhs: sources of eqs. 4-5, adiabatic solution with a
   Gaussian memory kernel of width sigma, random phases) is compared with
     - the closed-form golden-rule integrals (levelB_model.boltzmann), and
     - the Level B code itself: materials/graphite/dynamics.py's phph_rate (three kinds of
       three-phonon processes) and eph_rate (A(Te)[(n+1) exp(-w/Te) - n]), fed with process
       lists / emission tables built from the same couplings.
2. Dynamics, graphite protocol (thermalized hot electrons, electron energy bookkeeping):
   graphite's rate functions vs the solver's kinetic limit, N_a(t) compared.
3. Full (non-Markov, coherent) solver vs its kinetic limit for decreasing coupling:
   anharmonic decay of an optical phonon into M = 24 acoustic modes, and phonon heating by
   hot electrons in L = 40 levels.

Usage: python levelB_limit.py [part ...]   (default: 1 2 3)  -> data/levelB_limit.npz
"""
import os, sys, time
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, HERE)
from levelB_model import build, boltzmann, gauss, nambu_occ, UnifiedSolver

SIG = 0.1           # energy-conservation width (Gaussian memory kernel), model units


def import_graphite():
    """Import materials/graphite/dynamics.py (Level B) without running its main block."""
    gdir = os.path.join(ROOT, "materials", "graphite")
    sys.path.insert(0, gdir)
    argv = sys.argv; sys.argv = sys.argv[:1]
    import dynamics as D
    sys.argv = argv
    return D


def install_model(D, m, eps, Tgrid):
    """Point graphite's rate functions at the model: process lists and emission table.
    Energies of the model are read as eV (Te in K via KB_eV)."""
    Q, w = m.Q, m.omega
    k2 = m.kappa ** 2
    procs = []
    for kind in range(3):
        a, b, c = np.nonzero(k2 > 0)
        if kind == 0:
            dl = w[a] - w[b] - w[c]
        elif kind == 1:
            dl = w[a] + w[b] - w[c]
        else:
            dl = w[a] - w[b] + w[c]
        procs.append((a, b, c, np.pi / 8 * k2[a, b, c] * gauss(dl, SIG)))
    D.procs = procs
    D.Nq, D.nb = Q, 1
    D.hw = w[:, None].copy()
    D.valid = np.ones((Q, 1), bool)
    D.Tgrid = Tgrid; D.logT = np.log(Tgrid)
    A = np.zeros((Q, 1, len(Tgrid)))
    de = eps[:, None] - eps[None, :]
    for k, T in enumerate(Tgrid):
        f = 1 / (np.exp(eps / (D.KB_eV * T)) + 1)
        for q in range(Q):
            g2 = np.abs(m.gamma[q]) ** 2
            # absorption coefficient: j (lower) + phonon -> i (upper), both spins
            A[q, 0, k] = 2 * np.pi * np.sum(g2 * gauss(de - w[q], SIG) * f[None, :] * (1 - f[:, None]))
    D.Atab = A


def unified_kinetic(s, f, N):
    """kinetic_rhs in the original level basis (H_N diagonal): returns df (per spin), dN."""
    L = s.L
    H = s.hfb(0.0, s.G, np.zeros(s.Q))
    E, U = np.linalg.eigh(H)
    P = np.abs(U) ** 2
    df_e, dN, _, _ = s.kinetic_rhs(P.T @ nambu_occ(f), N, SIG, "gauss")
    dfn = P @ df_e
    return dfn[:L], dN


def part1(D):
    out = {}
    m, eps = build()
    s = UnifiedSolver(m, np.diag((eps < 0) * 1.0), np.diag((eps < 0) * 1.0), ee=None)
    rng = np.random.default_rng(7)
    f = rng.uniform(0, 1, m.L); N = rng.uniform(0, 2, m.Q)
    df, dN = unified_kinetic(s, f, N)
    bf, bN = boltzmann(m, eps, f, N, SIG)
    e1 = max(np.abs(df - bf).max() / np.abs(bf).max(), np.abs(dN - bN).max() / np.abs(bN).max())
    print(f"[1] random non-thermal state: solver kinetic limit vs golden rule, max rel. diff {e1:.1e}")
    Tgrid = np.geomspace(0.02, 5.0, 2000) / D.KB_eV
    install_model(D, m, eps, Tgrid)
    worst = 0.0
    for Te in (0.15, 0.3, 0.6, 1.0):
        fth = 1 / (np.exp(eps / Te) + 1)
        N = rng.uniform(0, 1, m.Q)
        df, dN = unified_kinetic(s, fth, N)
        s_pp = UnifiedSolver(m, s.G[:m.L, :m.L], s.G[:m.L, :m.L], ee=None, eph="meanfield")
        _, dN_pp = unified_kinetic(s_pp, fth, N)
        g_pp = D.phph_rate(N[:, None])[:, 0]
        g_ep = D.eph_rate(N[:, None], Te / D.KB_eV)[:, 0]
        r_pp = np.abs(dN_pp - g_pp).max() / np.abs(g_pp).max()
        r_ep = np.abs((dN - dN_pp) - g_ep).max() / np.abs(g_ep).max()
        worst = max(worst, r_pp)
        print(f"    thermal electrons Te = {Te}: phph rate vs graphite phph_rate {r_pp:.1e};"
              f" eph rate vs graphite eph_rate {r_ep:.1e} (detailed-balance form with broadened delta)")
        out[f"ep_err_{Te}"] = r_ep
    out["rand_err"] = e1; out["pp_err"] = worst
    return out


def part2(D, tmax=60.0, Te0=0.8, Tl=0.05):
    """Graphite protocol: thermal electrons (energy bookkeeping), phonon occupations."""
    m, eps = build()
    s = UnifiedSolver(m, np.diag((eps < 0) * 1.0), np.diag((eps < 0) * 1.0), ee=None)
    Tgrid = np.geomspace(0.02, 5.0, 2000) / D.KB_eV
    install_model(D, m, eps, Tgrid)
    Tt = np.geomspace(0.01, 5.0, 3000)
    Et = np.array([2 * np.sum(eps / (np.exp(eps / T) + 1)) for T in Tt])
    Te_of = lambda E: np.interp(E, Et, Tt)
    N0 = 1 / np.expm1(m.omega / Tl)
    y0 = np.concatenate([[np.interp(Te0, Tt, Et)], N0])
    ts = np.linspace(0, tmax, 121)

    def rhs_graphite(t, y):
        E, n = y[0], y[1:]
        Te = Te_of(E)
        r_ep = D.eph_rate(n[:, None], Te / D.KB_eV)[:, 0]
        r_pp = D.phph_rate(n[:, None])[:, 0]
        return np.concatenate([[-np.sum(m.omega * r_ep)], r_ep + r_pp])

    def rhs_unified(t, y):
        E, n = y[0], y[1:]
        f = 1 / (np.exp(eps / Te_of(E)) + 1)
        df, dN = unified_kinetic(s, f, n)
        return np.concatenate([[2 * np.sum(eps * df)], dN])
    t0 = time.time()
    a = solve_ivp(rhs_graphite, (0, tmax), y0, t_eval=ts, rtol=1e-8, atol=1e-10).y
    ta = time.time() - t0; t0 = time.time()
    b = solve_ivp(rhs_unified, (0, tmax), y0, t_eval=ts, rtol=1e-8, atol=1e-10).y
    tb = time.time() - t0
    Na, Nb = a[1:], b[1:]
    rel0 = np.abs(Na[0] - Nb[0]).max() / Na[0].max()
    relac = np.abs(Na[1:].sum(0) - Nb[1:].sum(0)).max() / Na[1:].sum(0).max()
    relT = np.abs(Te_of(a[0]) - Te_of(b[0])).max() / Te0
    print(f"[2] graphite protocol, t = 0..{tmax}: optical N0 peak {Na[0].max():.4f} (Level B code) vs "
          f"{Nb[0].max():.4f} (solver kinetic limit); max rel. diff N0 {rel0:.1e}, acoustic sum {relac:.1e}, "
          f"Te {relT:.1e}  [{ta:.1f}s, {tb:.1f}s]")
    return dict(p2_t=ts, p2_graphite=a, p2_unified=b, p2_Te_graphite=Te_of(a[0]), p2_Te_unified=Te_of(b[0]),
                p2_err=max(rel0, relac))


def _run_pair(m, f0, n0, tmax, dt=0.05, nout=81):
    """Full coherent solver and its kinetic limit from the same initial state."""
    L = m.L
    G0 = np.diag(f0).astype(complex)
    s = UnifiedSolver(m, G0, G0, ee=None, nph0=n0)
    ts = np.linspace(0, tmax, nout)

    def rhs(t, y):
        df, dN = unified_kinetic(s, y[:L], y[L:])
        return np.concatenate([df, dN])
    t0 = time.time()
    kin = solve_ivp(rhs, (0, tmax), np.concatenate([f0, n0]), t_eval=ts, rtol=1e-7, atol=1e-10).y.T
    tk = time.time() - t0; t0 = time.time()
    full = []
    nsub = int(round((ts[1] - ts[0]) / dt))
    for k in range(len(ts)):
        if k:
            for _ in range(nsub):
                s.step(ts[1] / nsub)
        full.append(s.phonon_number())
    return ts, kin[:, L:], np.array(full), s.energy() - s.E0, tk, time.time() - t0


def _smooth(y, ts, width):
    """running average over a window of the given width (removes the 2w coherent beats)."""
    k = max(1, int(round(width / (ts[1] - ts[0]))))
    pad = np.concatenate([y[:1].repeat(k // 2), y, y[-1:].repeat(k // 2)])
    return np.convolve(pad, np.ones(k) / k, mode="valid")[:len(y)]


def part3(scales_pp=(0.3, 0.2, 0.14, 0.1), scales_ep=(1.0, 0.7, 0.5, 0.35, 0.25, 0.18), M=24, L=40):
    """Full coherent solver vs its kinetic limit for decreasing coupling.
    pp: optical mode (N0 = 1) decaying into M acoustic modes (no electrons coupled);
    ep: phonon heating by hot electrons (Te = 1) in L levels (no anharmonicity).
    Time window: 4 / (initial kinetic rate); error measures: max |N_full - N_kin| / max N_kin,
    and the same after a running average over 2 pi / w0 (coherent beats removed)."""
    out = {}
    for kind, scales in (("pp", scales_pp), ("ep", scales_ep)):
        for sc in scales:
            if kind == "pp":
                m, eps = build(L=2, M=M, g=0.0, k0=0.1 * sc, jitter=0.8)
                f0 = np.zeros(2); n0 = np.zeros(m.Q); n0[0] = 1.0
                obs = lambda N: N[:, 0]
            else:
                m, eps = build(L=L, M=1, g=0.05 * sc, k0=0.0, jitter=0.8)
                f0 = 1 / (np.exp(eps / 1.0) + 1); n0 = np.zeros(m.Q)
                obs = lambda N: N[:, 0]
            s0 = UnifiedSolver(m, np.diag(f0), np.diag(f0), ee=None, nph0=n0)
            _, dN = unified_kinetic(s0, f0, n0)
            Nfin = 1.0 if kind == "pp" else 0.5
            rate = abs(dN[0]) / Nfin
            tmax = min(4.0 / rate, 600.0)
            ts, Nk, Nf, edrift, tk, tf = _run_pair(m, f0, n0, tmax)
            a, b = obs(Nk), obs(Nf)
            scale = np.abs(a - a[0]).max()
            e_max = np.abs(a - b).max() / scale
            w = 2 * np.pi / m.omega[0]
            msk = ts > min(w, ts[-1] / 4)
            e_avg = np.abs(_smooth(a, ts, w) - _smooth(b, ts, w))[msk].max() / scale
            print(f"[3{kind}] coupling x{sc}: rate {rate:.3f}, t = 0..{tmax:.0f}; N0 kinetic {a[-1]:.4f}, full {b[-1]:.4f};"
                  f" max dev {e_max:.3f}, beat-averaged {e_avg:.3f} (of the change in N0); "
                  f"min N {Nf.min():.3f}; energy drift {edrift:.1e} [{tk:.0f}s, {tf:.0f}s]", flush=True)
            out[f"p3{kind}_{sc}"] = dict(t=ts, kin=Nk, full=Nf, e_max=e_max, e_avg=e_avg, rate=rate)
    return out


if __name__ == "__main__":
    parts = sys.argv[1:] or ["1", "2", "3"]
    res = {}
    D = import_graphite() if ("1" in parts or "2" in parts) else None
    if "1" in parts:
        res.update(part1(D))
    if "2" in parts:
        res.update(part2(D))
    if "3" in parts:
        for k, v in part3().items():
            for kk, vv in v.items():
                res[f"{k}_{kk}"] = vv
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    fn = os.path.join(HERE, "data", "levelB_limit_" + "".join(parts) + ".npz")
    np.savez(fn, **res)
