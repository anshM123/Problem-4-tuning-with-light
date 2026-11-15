"""Fast checks of the LiNbO3 nonlinear-phononics model (materials/phononics)."""
import os, sys
import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PH = os.path.join(HERE, "..", "materials", "phononics")
sys.path.insert(0, PH)
need = pytest.mark.skipif(not all(os.path.exists(os.path.join(PH, f)) for f in
                                  ("anharmonic.npz", "phph.npz", "gamma_modes.npz")),
                          reason="phononics intermediates not generated")


def test_fluence_to_field():
    from params import fluence_to_field
    assert abs(fluence_to_field(95.0) - 21.2) < 0.2


@need
def test_forces_are_gradient_of_potential():
    import drive as d
    rng = np.random.default_rng(0)
    Q = rng.uniform(-1.5, 1.5, (4, 9))
    F = d.forces(Q)
    h = 1e-5
    for j in range(4):
        dq = np.zeros_like(Q); dq[j] = h
        Fn = -(d.potential(Q + dq) - d.potential(Q - dq)) / (2 * h)
        assert np.allclose(F[j], Fn, atol=1e-7)


@need
def test_bath_step_conserves_energy():
    import drive as d
    b = d.Bath(3)
    E0 = b.energies()
    P4 = np.array([1e-6, 0.0, 3e-7])          # J/s per cell
    heat = np.array([0.0, 2e-22, 1e-22])      # J per cell
    dt = 3e-15
    for _ in range(20):
        b.step(dt, P4, heat)
    assert np.allclose(b.energies() - E0, 20 * (P4 * dt + heat), rtol=1e-9, atol=1e-30)
    assert np.all(b.T >= 300.0)


@need
def test_linear_drive_energy_budget():
    import drive as d
    r = d.run(1.0, 17.5, t_end_ps=0.8, record=False, quiet=True)
    bud = r["budget"]
    tot = bud["Uref"] + bud["Uabs"] + bud["Ufield"] + bud["Umech"] + bud["Udiss"]
    assert abs(tot / bud["Uin"] - 1) < 5e-3


@need
def test_literature_surface_swap_is_exact_and_reversible():
    import drive as d
    import diagnostics as g
    rng = np.random.default_rng(1)
    Q = rng.uniform(-1.0, 1.0, (4, 6))
    F_mace = d.forces(Q).copy()
    Vlit, P0, I0 = g.make_surface("L2")
    d.set_surface("literature", Vlit, Q0=P0)
    Qs = Q.copy(); Qs[1:3] = 0.0               # TO2/TO3 at rest: V = Vlit(q1, q4) exactly
    assert np.allclose(d.potential(Qs), Vlit(Qs[0], Qs[3]), atol=1e-8)
    d.set_surface("mace")
    assert np.allclose(d.forces(Q), F_mace)
