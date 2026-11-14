"""Fast consistency tests of the core methods (run: python -m pytest tests -q)."""
import os, sys
import numpy as np
ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, ROOT)
from lightmat.ed import Hubbard, chain_bonds, evolve
from lightmat.fluct_hubbard import PPLadder
from lightmat.fluct_ph import PHFluct


def _chain(L):
    h = np.zeros((L, L))
    for i in range(L - 1):
        h[i, i + 1] = h[i + 1, i] = -1.0
    return h


def test_ed_dimer_ground_state():
    m = Hubbard(2, 1, 1, [(0, 1)], U=4.0)
    E, _ = m.ground_state()
    assert abs(E - (4 - np.sqrt(16 + 16)) / 2) < 1e-10


def test_factorized_equals_g1g2_without_collisions():
    """The positive-factorized pair dynamics is algebraically identical to G1-G2 (T-matrix)."""
    L, U = 5, 3.0
    h = _chain(L); r0 = np.diag([0, 1, 0, 1, 0.0])
    out = []
    for mode in ("factorized", "g1g2"):
        s = PPLadder(h, U, r0, r0, mode=mode)
        orig = s.rhs
        def rhs(t, st, orig=orig, s=s):
            k = orig(t, st)
            ru, rd = st[0], st[1]
            hu = h + U * np.diag(np.real(np.diag(rd))); hd = h + U * np.diag(np.real(np.diag(ru)))
            k[0] = -1j * (hu @ ru - ru @ hu); k[1] = -1j * (hd @ rd - rd @ hd)
            return k
        s.rhs = rhs
        t = 0.0
        for _ in range(200):
            s.step(t, 0.005); t += 0.005
        out.append(s.dcorr())
    assert np.max(np.abs(out[0] - out[1])) < 1e-6


def test_weak_coupling_vs_exact():
    L, U = 6, 0.5
    occ = [0, 2, 4]
    m = Hubbard(L, 3, 3, chain_bonds(L), U)
    psi = m.product_state(occ, occ)
    ex = evolve(m, psi, np.array([0.0, 1.0]), observe=lambda p, t: m.double_occ(p))
    r0 = np.diag([1.0 if i in occ else 0.0 for i in range(L)])
    s = PPLadder(_chain(L), U, r0, r0)
    t = 0.0
    for _ in range(200):
        s.step(t, 0.005); t += 0.005
    nu, nd = s.densities()
    assert np.max(np.abs(nu * nd + s.dcorr() - ex[1])) < 0.02


def test_positivity_of_pair_matrices():
    """X = x x^+ and Y = y y^+ stay positive semidefinite at strong coupling."""
    L, U = 7, 15.0
    r0 = np.diag([1.0 if i % 2 else 0.0 for i in range(L)])
    s = PPLadder(_chain(L), U, r0, r0)
    t = 0.0
    for _ in range(400):
        s.step(t, 0.002); t += 0.002
    for f in (s.x, s.y):
        w = np.linalg.eigvalsh(f @ f.conj().T)
        assert w.min() > -1e-10
    assert np.all(np.isfinite(s.dcorr()))


def test_sparse_ph_kernel_matches_dense():
    L = 4
    bonds = [(i, i + 1, 1.0, -1.0) for i in range(L - 1)]
    rng = np.random.default_rng(0)
    A = rng.standard_normal((8, 8)) + 1j * rng.standard_normal((8, 8))
    rho = A @ A.conj().T; rho /= np.trace(rho).real / 4
    s = PHFluct(L, bonds, 2.0, np.diag([1, 0, 1, 0, 0, 1, 0, 1.0]), kernel="gwx", nsamp=5)
    V = s.V + 0.1 * (rng.standard_normal(s.V.shape) + 1j * rng.standard_normal(s.V.shape))
    V = V + np.conj(np.transpose(V, (0, 2, 1)))
    _, dV = s.rhs(0.0, rho, V)
    U = 2.0; hs = s.hop(0).toarray(); p = s.partner; a = np.arange(8)
    H = np.zeros((8, 8), complex); H[a, a] = U * rho[p, p]; H[a, p] = -U * rho[a, p]
    dH = np.zeros_like(V); dH[:, a, a] = U * V[:, p, p]; dH[:, a, p] = -U * V[:, a, p]
    ref = -1j * ((hs + H)[None] @ V - V @ (hs + H)[None] + dH @ rho[None] - rho[None] @ dH)
    assert np.max(np.abs(ref - dV)) < 1e-12


def test_mgb2_gaps():
    sys.path.insert(0, os.path.join(ROOT, "materials", "mgb2"))
    import twoband as tb
    wc = tb.fit_wc()
    D = tb.gaps(5.0, wc)
    assert abs(2 * D[0] / tb.THZ - 1.07) < 0.03 and abs(2 * D[1] / tb.THZ - 3.40) < 0.05


def test_wse2_equilibrium_is_stationary():
    sys.path.insert(0, os.path.join(ROOT, "materials", "wse2"))
    import dynamics as wd
    P = dict(wd.PARAMS); P["dE_SK"] = 0.0068; P["g_EEA"] = 0.0
    m = wd.build(P)
    eg = m["eg"]; kT = wd.KB * P["T"]
    N = np.array([m["g"][s] * m["M"][s] * np.exp(-(eg + m["E0"][s]) / kT) for s in m["sp"]]) * 1e10
    d = wd.rhs_factory(m, lambda t: 0.0)(0.0, N.reshape(-1)).reshape(N.shape)
    # interior of the energy grid (edges lose detailed-balance partners at the cutoff);
    # the residual (< 1 % of the outflow) comes from rounding phonon energies to the grid
    out = sum(w[2] for w in m["W"]) * N
    assert np.max(np.abs(d[:, :200]) / out[:, :200]) < 0.02
