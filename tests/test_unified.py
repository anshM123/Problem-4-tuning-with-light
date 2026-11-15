"""Tests of the unified electron-phonon solver (lightmat/unified.py) and its three limits.
Run: python -m pytest tests/test_unified.py -q   (about one minute)"""
import os, sys
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "benchmarks", "unified_limits"))
sys.path.insert(0, os.path.join(ROOT, "materials", "mgb2"))
from lightmat.unified import Model, UnifiedSolver, peierls_chain
from lightmat.fluct_hubbard import PPLadder


def _cdw(L=7):
    h = np.zeros((L, L))
    for i in range(L - 1):
        h[i, i + 1] = h[i + 1, i] = -1.0
    r0 = np.diag([float(i % 2) for i in range(L)])
    return h, r0


def test_level_A_limit_is_ppladder():
    """Without phonons the unified solver is Level A (both pair ladders) to rounding."""
    h, r0 = _cdw(7)
    for mode in ("factorized", "sic"):
        a = PPLadder(h, 5.0, r0, r0, mode=mode)
        b = UnifiedSolver(Model(h, U=5.0), r0, r0, ee=mode)
        t = 0.0
        for _ in range(150):
            a.step(t, 0.01); t += 0.01
            b.step(0.01)
        na, da = a.densities()
        nb, db = b.densities()
        assert np.abs(na - nb).max() < 1e-12 and np.abs(da - db).max() < 1e-12
        assert np.abs(a.dcorr() - b._dcorr(b.G, b.x, b.y)).max() < 1e-12


def _dimer(g, lam, U=2.0, A0=0.6):
    import full_vs_ed as F
    L, Q, omega, gamma, kappa = F.build(g, lam, U)
    A = lambda t: A0 * np.cos(2.0 * (t - 3.0)) * np.exp(-(t - 3.0) ** 2 / 2.0)
    return Model(peierls_chain(L, 1.0, A), U=U, omega=omega, gamma=gamma, kappa=kappa), A


def test_full_equations_conserve_energy():
    """Eqs. (1)-(6) with e-e, e-ph and ph-ph correlations conserve <H> up to the field's work."""
    m, _ = _dimer(0.2, 0.1)
    r = np.diag([1.0, 0.0]).astype(complex)
    s = UnifiedSolver(m, r, r)
    for _ in range(600):
        s.step(0.01)
    assert abs(s.energy() - s.W - s.E0) < 1e-7
    assert abs(s.W) > 1e-2          # the pulse did work


def _dimer_errors(U, g, lam, kws, tmax=8.0):
    import full_vs_ed as F
    from ed_eph import EPhED
    m, A = _dimer(g, lam, U=U)
    L, Q, omega, gamma, kappa = F.build(g, lam, U)
    ed = EPhED(L, 1, 1, [(0, 1)], 1.0, U, omega, gamma, kappa, [7, 7, 9, 9], Afun=A)
    times = np.arange(0, tmax + 0.01, 0.5)
    ex = ed.run(ed.initial([0], [0]), times, nsub=10)
    r = np.diag([1.0, 0.0]).astype(complex)
    out = []
    for kw in kws:
        s = UnifiedSolver(m, r, r, **kw)
        err = []
        for k, T in enumerate(times):
            if k:
                for _ in range(50):
                    s.step(0.01)
            e = ex[k]; N = s.phonon_number()
            err.append([abs(s.double_occ().sum() - sum(e["d"])), abs(N[0] + N[1] - e["N"][0] - e["N"][1]),
                        abs(N[2] + N[3] - e["N"][2] - e["N"][3])])
        out.append(np.array(err).max(0))
    return out


def test_full_level_exact_to_second_order_without_ee():
    """U = 0 dimer, Holstein + cubic phonons + light: errors are O(g^4) (Ehrenfest: O(g^2))."""
    full, ehr = _dimer_errors(0.0, 0.1, 0.1, [dict(ee="hf"), dict(ee="hf", eph="meanfield", phph="meanfield")])
    assert full.max() < 6e-3
    assert np.all(full < 0.3 * ehr)


def test_full_level_all_interactions_dimer():
    """U = 2, g = 0.1: e-e (exact pair ladder), e-ph, ph-ph correlations and light at once.
    Both phonon channels improve on Ehrenfest (soft, decay-product modes 4x); D stays within 0.08
    of exact.  (At g = 0.2 the closure drifts after t ~ 8, see benchmarks/unified_limits.)"""
    full, ehr = _dimer_errors(2.0, 0.1, 0.1, [{}, dict(eph="meanfield", phph="meanfield")])
    assert full[0] < 0.08
    assert full[1] < ehr[1] and full[2] < 0.5 * ehr[2]


def test_kinetic_limit_is_golden_rule():
    """Markov + random-phase limit of eqs. (4), (5) = Boltzmann e-ph and three-phonon integrals."""
    from levelB_model import build, boltzmann, nambu_occ
    m, eps = build(L=12, M=6)
    s = UnifiedSolver(m, np.diag((eps < 0) * 1.0), np.diag((eps < 0) * 1.0), ee=None)
    rng = np.random.default_rng(3)
    f = rng.uniform(0, 1, m.L); N = rng.uniform(0, 1.5, m.Q)
    H = s.hfb(0, s.G, np.zeros(m.Q)); E, U = np.linalg.eigh(H)
    P = np.abs(U) ** 2
    dfe, dN, _, _ = s.kinetic_rhs(P.T @ nambu_occ(f), N, 0.1, "gauss")
    df = (P @ dfe)[:m.L]
    bf, bN = boltzmann(m, eps, f, N, 0.1)
    assert np.abs(df - bf).max() < 1e-12 * max(1, np.abs(bf).max()) + 1e-14
    assert np.abs(dN - bN).max() < 1e-12


def test_level_C_gap_equation():
    """Nambu mean field with the MgB2 pairing interaction = Level C two-band gap equation."""
    import twoband as TB
    import levelC_limit as C
    wc = TB.fit_wc()
    n, T = 100, 5.0
    Dc = C.levelC_gaps_same_levels(n, wc, T)
    m, xi, band, nb = C.build(n, wc)
    E = np.sqrt(xi ** 2 + Dc[band] ** 2)
    kap = np.diag(-Dc[band] / (2 * E) * np.tanh(E / (2 * TB.KB * T)) * 0.9)     # perturbed start
    s = UnifiedSolver(m, np.diag((xi < 0) * 1.0), np.diag((xi < 0) * 1.0), kappa=kap, ee=None)
    D = np.abs(s.equilibrium(TB.KB * T, tol=1e-12, mix=1.0))
    Ds = np.array([D[band == a].mean() for a in range(2)])
    assert np.abs(Ds - Dc).max() < 1e-9
    assert np.abs(Ds - TB.gaps(T, wc)).max() < 0.01        # discretization error vs continuum


def test_level_C_amplitude_mode():
    """Linearized dynamics around the gap = Level C amplitude propagator (smaller grid, its own
    level-density ratio N_s/N_p = 15/20)."""
    import twoband as TB
    wc = TB.fit_wc()
    nb, T, eta = (20, 15), 5.0, 0.3
    nr = np.array([1.0, nb[1] / nb[0]])
    xi, band = [], []
    for a in range(2):
        x = (np.arange(nb[a]) + 0.5) * wc / nb[a]
        xi += list(-x[::-1]) + list(x); band += [a] * (2 * nb[a])
    xi, band = np.array(xi), np.array(band)
    Vlev = -TB.V[band][:, band] * wc / nb[0]
    # Level C gap equation with this level-density ratio (twoband._I on the same levels)
    Dc = np.array([2.3, 7.0])
    for _ in range(5000):
        Dc = 0.5 * Dc + 0.5 * TB.V @ (nr * np.array([TB._I(Dc[a], T, wc, n=nb[a]) for a in range(2)]) * Dc)
    E = np.sqrt(xi ** 2 + Dc[band] ** 2)
    kap = np.diag(-Dc[band] / (2 * E) * np.tanh(E / (2 * TB.KB * T)))
    h0, tp, sp = 1e-3, 1.0, 0.25
    sel = (band == 0).astype(float)
    pf = lambda t: np.diag(sel * h0 * np.exp(-(t - tp) ** 2 / (2 * sp ** 2)))
    m = Model(np.diag(xi), V=Vlev, pair_field=pf)
    s = UnifiedSolver(m, np.diag((xi < 0) * 1.0), np.diag((xi < 0) * 1.0), kappa=kap, ee=None)
    s.equilibrium(TB.KB * T, tol=1e-13, mix=1.0)
    base = np.array([np.real(s.gap()[band == a]).mean() for a in range(2)])
    dt, nt = 0.025, 1200
    ts = np.arange(nt + 1) * dt
    dD = []
    for k in range(nt + 1):
        g = s.gap()
        dD.append([np.real(g[band == a]).mean() - base[a] for a in range(2)])
        if k < nt:
            s.step(dt)
    dD = np.array(dD)
    w = np.linspace(2.0, 12.0, 101)          # probe-pulse spectrum > 1e-2 of its peak
    wt = np.exp(-eta * ts) * dt
    ft = lambda y: np.array([np.sum(y * wt * np.exp(1j * x * ts)) for x in w])
    hw = ft(h0 * np.exp(-(ts - tp) ** 2 / (2 * sp ** 2)))
    Rs = ft(dD[:, 0]) / hw
    Rc = []
    for x in w:
        P = np.array([TB.P_amp(x, abs(base[a]), T, wc, eta=eta, n=nb[a]) for a in range(2)])
        Dm = np.linalg.inv(np.linalg.inv(TB.V) - np.diag(nr * P))
        Rc.append((Dm @ np.linalg.inv(TB.V) - np.eye(2))[0, 0])
    Rc = np.array(Rc)
    assert np.abs(Rs - Rc).max() / np.abs(Rc).max() < 0.02
    assert abs(w[np.argmax(np.abs(Rs))] - w[np.argmax(np.abs(Rc))]) < 0.15


def test_material_scale_markov_operators():
    """Large-scale Markov operators used for graphite: e-e and e-ph conserve particle number
    and energy and vanish in equilibrium; the three-phonon operator on process lists equals the
    solver's own kinetic limit (kinetic_rhs)."""
    from lightmat.unified import markov_ee_local, EphKernel, markov_phph
    from levelB_model import build, gauss
    de = 0.05
    eg = np.arange(-3, 3 + 1e-9, de)
    rho = 0.1 + 0.05 * np.abs(eg)
    g = 2 * rho * de
    rng = np.random.default_rng(0)
    F = np.clip(1 / (np.exp(eg / 0.3) + 1) + 0.1 * rng.uniform(-1, 1, len(eg)) * np.exp(-eg ** 2), 0, 1)
    gin, gout = markov_ee_local(F, rho, 2.0, de)
    dF = gin * (1 - F) - gout * F
    assert abs(np.sum(g * dF)) < 1e-10 * np.sum(g * np.abs(dF))
    assert abs(np.sum(g * eg * dF)) < 1e-10 * np.sum(g * np.abs(eg * dF))
    Fth = 1 / (np.exp(eg / 0.3) + 1)
    gin, gout = markov_ee_local(Fth, rho, 2.0, de)
    assert np.abs(gin * (1 - Fth) - gout * Fth).max() < 1e-12 * gout.max()
    w = np.array([0.07, 0.123, 0.2])
    W = rng.uniform(0, 1, (3, len(eg)))
    k = EphKernel(W, w, de, g, 10)
    N = np.array([0.3, 1.0, 0.1])
    dN, Gin, Gout = k.rates(F, N)
    dF = Gin * (1 - F) - Gout * F
    assert abs(np.sum(g * dF)) < 1e-12 * np.abs(dN).sum()
    assert abs(np.sum(g * eg * dF) + np.sum(w * dN) / 10) < 1e-12 * np.abs(w * dN).sum()
    m, eps = build(L=2, M=6, g=0.0, k0=0.1)
    s = UnifiedSolver(m, np.zeros((2, 2)), np.zeros((2, 2)), ee=None)
    Nph = rng.uniform(0, 1.5, m.Q)
    _, dN_kin, _, _ = s.kinetic_rhs(np.concatenate([np.zeros(2), np.ones(2)]), Nph, 0.1, "gauss")
    k2, wq = m.kappa ** 2, m.omega
    a, b, c = np.nonzero(k2 > 0)
    procs = [(a, b, c, np.pi / 8 * k2[a, b, c] * gauss(dl, 0.1))
             for dl in (wq[a] - wq[b] - wq[c], wq[a] + wq[b] - wq[c], wq[a] - wq[b] + wq[c])]
    assert np.abs(markov_phph(Nph, procs) - dN_kin).max() < 1e-12 * np.abs(dN_kin).max()
