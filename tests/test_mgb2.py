"""Fast checks of the MgB2 real-time model (run: python -m pytest tests/test_mgb2.py -q)."""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "materials", "mgb2"))
import realtime as rt  # noqa: E402


def test_eliashberg_table_couplings():
    """digitized two-band Eliashberg functions reproduce the published lambda_ab."""
    w = rt.A2F_W
    lam = [[2 * np.trapezoid(rt.A2F[a][b] / w, w) for b in range(2)] for a in range(2)]
    assert np.allclose(lam, [[0.448, 0.155], [0.213, 1.017]], rtol=0.03)
    # lambda_ab / N_b symmetric (same matrix element for a->b and b->a)
    assert abs(lam[0][1] / 0.300 - lam[1][0] / 0.410) < 0.03


def test_random_matrix_coupling_gives_drude():
    """Kubo conductivity of the random velocity operator equals sigma_n L(w) (normal state, 7 K)."""
    b = rt.Band(0, 300, 15.0, 3)
    e = b.eps; f = rt.fermi(e, 7.0)
    de = e[None, :] - e[:, None]; df = f[:, None] - f[None, :]
    for fr in (0.5, 1.0):
        w = fr * rt.THZ; sg = 0.3
        dl = np.exp(-(de - w) ** 2 / (2 * sg ** 2)) / np.sqrt(2 * np.pi) / sg
        s1 = 2 * np.pi * b.cJ ** 2 / (b.Vm * fr * 2 * np.pi * 1e12) * np.sum(b.K ** 2 * df * dl) / rt.MEV
        ref = rt.SIGN[0] * rt.GAM[0] ** 2 / (rt.GAM[0] ** 2 + w ** 2)
        assert abs(s1 / ref - 1) < 0.1


def test_thermal_state_stationary():
    s = rt.Sim(N=(120, 100), W=(12.0, 16.0), collisions=False)
    D0 = s.gaps(s.R)
    for _ in range(30):
        s.step(0.01, lambda t: 0.0)
    assert np.max(abs(s.gaps(s.R) - D0)) < 1e-4


def test_collisions_detailed_balance_and_energy():
    s = rt.Sim(N=(120, 100), W=(12.0, 16.0))
    D = s.gaps(s.R).real
    fs = [rt.level_f(b, Rb, D[b.a])[0] for b, Rb in zip(s.bands, s.R)]
    df, gl, dn = rt.kinetics(s.bands, fs, D, s.ph)
    assert max(np.max(abs(x)) for x in df) < 1e-4        # thermal f and n at the same T: no net flux
    # hot quasiparticles: energy lost by electrons equals energy gained by phonons (no escape)
    fs = [rt.fermi(np.sqrt(b.eps ** 2 + D[b.a] ** 2), 20.0) for b in s.bands]
    tau = rt.TAU_ESC
    rt.TAU_ESC = np.inf
    try:
        df, gl, dn = rt.kinetics(s.bands, fs, D, s.ph)
    finally:
        rt.TAU_ESC = tau
    dUe = sum(2 * rt.NUD[b.a] * b.de * np.sum(np.sqrt(b.eps ** 2 + D[b.a] ** 2) * df[b.a]) for b in s.bands)
    dUp = np.sum(dn * s.ph.D * s.ph.W) * s.ph.dW
    assert abs(dUe + dUp) < 1e-4 * abs(dUp)        # residual: transfers below the first phonon bin
