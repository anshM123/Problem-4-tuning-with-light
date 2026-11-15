"""Checks of the chain optical-conductivity machinery (run: python -m pytest tests -q)."""
import os, sys
import numpy as np
ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "benchmarks", "hubbard"))
from chain_sigma_ed import Ring, sigma_from_poles, bethe_gap
from lightmat.fluct_ph import DoublonHolonPairs


def test_lanczos_sigma_matches_full_diagonalization():
    R = Ring(6, 4.0, "pbc")
    E0, g = R.ground_state()
    H = np.array([R.H(np.eye(R.n)[:, k]) for k in range(R.n)]).T
    e, v = np.linalg.eigh(H)
    assert abs(E0 - e[0]) < 1e-9
    phi = R.J(g)
    el, al = R.lanczos(phi, 200)
    w = np.linspace(0.1, 12, 300)
    s_full = sigma_from_poles(w, e - E0, (v.T @ phi) ** 2, 6, 0.1)
    s_lanc = sigma_from_poles(w, el - E0, al, 6, 0.1)
    assert np.abs(s_full - s_lanc).max() < 1e-8 * s_full.max()


def test_doublon_holon_exact_at_large_U():
    """At U >> J the doublon-holon closure (Neel background, exact f-sum rule) reproduces the
    exact optical spectrum of the open chain."""
    L, U = 8, 200.0
    R = Ring(L, U, "obc")
    E0, g = R.ground_state()
    e, a = R.lanczos(R.J(g), 300)
    w = np.linspace(U - 6, U + 6, 600)
    se = sigma_from_poles(w, e - E0, a, L, 0.2)
    kb = np.full(L - 1, -8 * 0.5 / U)                 # Neel singlet weight 1/2
    eb = np.full(L - 1, -2.0 / U)
    m = DoublonHolonPairs(L, U, kb, eb)
    ee, vv = np.linalg.eigh(m.H.toarray())
    sm = sigma_from_poles(w, ee, (vv.T @ m.phi) ** 2, L, 0.2)
    dist = np.abs(sm / sm.sum() - se / se.sum()).sum()
    assert dist < 0.25


def test_time_domain_response_equals_pole_sum():
    L, U = 10, 12.0
    kb = np.full(L - 1, -0.33); eb = np.full(L - 1, -0.165)
    m = DoublonHolonPairs(L, U, kb, eb)
    t, j = m.current_response(T=60.0, dt=0.01)
    w = np.linspace(U - 6, U + 6, 200)
    s_t = np.array([-np.sum(j * np.exp((1j * x - 0.1) * t)) * 0.01 for x in w]).real / L
    ee, vv = np.linalg.eigh(m.H.toarray())
    s_p = sigma_from_poles(w, ee, (vv.T @ m.phi) ** 2, L, 0.1)
    assert np.abs(s_t - s_p).max() < 0.02 * s_p.max()
    # exact f-sum rule: sum_n a_n / w_n = -<K> / 2
    an = (vv.T @ m.phi) ** 2
    assert abs(2 * np.sum(an / ee) + kb.sum()) < 1e-10


def test_doublon_holon_band_edge_is_strong_coupling_gap():
    """With the exact bond energy and singlet weight of the Heisenberg chain the second-order
    closure has the exact strong-coupling Mott gap U - 4J + 8 ln2 J^2/U."""
    from scipy.sparse.linalg import eigsh
    L, U = 120, 20.0
    kb = np.full(L, -8 * np.log(2) / U); eb = np.full(L, -4 * np.log(2) / U)
    m = DoublonHolonPairs(L, U, kb, eb, periodic=True)
    emin = eigsh(m.H, k=1, which="SA")[0][0]
    assert abs(emin - (U - 4 + 8 * np.log(2) / U)) < 0.01
    assert abs(bethe_gap(U) - (U - 4 + 8 * np.log(2) / U)) < 0.01
