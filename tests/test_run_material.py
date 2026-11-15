"""Fast tests of the general pipeline run_material.py (no MACE calls; run: python -m pytest tests -q)."""
import os, sys
import numpy as np
import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, ROOT)
import run_material as rm


def _graphene_tb():
    at, _, _ = rm.load_structure("graphite")
    pd = rm.periodic_dims(at)
    return at, pd, rm.UniversalTB(at.cell, at.get_scaled_positions(), at.get_chemical_symbols(), at.get_masses(), pd)


def test_dimensionality_detection():
    at, _, _ = rm.load_structure("graphite")
    assert rm.periodic_dims(at).tolist() == [True, True, False]
    bp, _, _ = rm.load_structure("black_phosphorus")
    assert rm.periodic_dims(bp).tolist() == [True, True, True]
    assert len(bp) == 8


def test_harrison_graphene_dirac_cone():
    """Universal sp3 model: Dirac point at K at the p level, hbar v_F within 3 % of 6.4 eV A."""
    _, _, tb = _graphene_tb()
    K = np.array([1 / 3, 1 / 3, 0]) @ tb.rec
    e0, _ = tb.eig(K[None])
    assert abs(e0[0, 3] - e0[0, 4]) < 1e-8 and abs(e0[0, 4] - tb.onsite[1]) < 1e-8
    e1, _ = tb.eig((K + [0.005, 0, 0])[None])
    vF = (e1[0, 4] - e0[0, 4]) / 0.005
    assert abs(vF - 6.4) / 6.4 < 0.03
    # bond-stretching derivative of the pi hopping: 2|V_pp_pi|/d = 4.3 eV/A (graphite DFT-calibrated model: 4.68)
    d = 2.464 / np.sqrt(3)
    beta = 2 * abs(rm.harrison_V("ppp", d)) / d
    assert 4.1 < beta < 4.5


def test_eph_vertex_matches_frozen_phonon():
    """g(k, q=0) from the analytic vertex equals the finite-difference change of the Bloch Hamiltonian
    under a frozen optical displacement."""
    at, pd, tb = _graphene_tb()
    k = np.array([[0.31, 0.17, 0.0]]) @ tb.rec
    E, C = tb.eig(k)
    # optical in-plane mode at Gamma: atoms move oppositely along x (unit mass-weighted vector)
    e = np.zeros((2, 3)); e[0, 0] = 1 / np.sqrt(2); e[1, 0] = -1 / np.sqrt(2)
    freq = np.array([48.0])                       # THz, value only sets the zero-point length
    W = None
    # analytic: rebuild the pieces used in eph_W_point
    w_rad = 2 * np.pi * freq * 1e12
    l = np.sqrt(rm.HBAR_Js / (2 * tb.masses[:, None] * rm.AMU * w_rad[None, :])) * 1e10
    U = (e.reshape(2, 3, 1) * l[:, None, :]).reshape(6, 1)
    M = (np.exp(1j * k @ tb.R.T) @ tb.DP.reshape(len(tb.R), -1) - np.exp(1j * k @ tb.R.T) @ tb.DM.reshape(len(tb.R), -1))
    M = M.reshape(6, tb.norb, tb.norb)
    g_an = np.einsum("am,xab,bn,x->mn", C[0].conj(), M, C[0], U[:, 0])
    # finite difference: displace atoms by +-h * u
    h = 1e-4
    Hs = []
    for sgn in (1, -1):
        pos = at.positions.copy() + sgn * h * (U[:, 0].reshape(2, 3)).real
        frac = np.linalg.solve(np.array(at.cell).T, pos.T).T
        tb2 = rm.UniversalTB(at.cell, frac, at.get_chemical_symbols(), at.get_masses(), pd)
        Hs.append(tb2.hk(k)[0])
    g_fd = C[0].conj().T @ ((Hs[0] - Hs[1]) / (2 * h)) @ C[0]
    assert np.allclose(g_an, g_fd, atol=1e-5 * np.abs(g_fd).max())


def test_film_absorption_graphite_repo_value():
    A, R, T = rm.film_absorption(2.72 + 1.47j, 70e-9, 800e-9)
    assert abs(A - 0.4825) < 1e-3 and abs(A + R + T - 1) < 1e-12


def test_fit_recovers_biexponential():
    t = np.linspace(-0.5, 60, 3000)
    y = rm.MODELS["biexp"][0](t, 2.0, 0.4, 1.0, 8.0)
    tu = np.arange(t[0], t[-1], 0.01)
    yc = np.interp(t, tu, rm.gauss_conv(tu, np.interp(tu, t, y), 0.15))
    p = rm.fit_trace(t, yc, "biexp", 0.15, 0.0)
    assert abs(p["tau1_ps"] - 0.4) < 0.03 and abs(p["tau2_ps"] - 8.0) < 0.3


def test_effective_temperature_inverts_bose():
    hw = np.linspace(0.01, 0.2, 50)
    T = np.array([150.0, 800.0])
    E = np.array([np.sum(hw / np.expm1(hw / (rm.KB_eV * x))) for x in T])
    assert np.allclose(rm.eff_temperature(hw, E), T, rtol=1e-6)


def test_input_validation():
    with pytest.raises(SystemExit):
        rm.validate_args(rm.parse_args(["graphite", "--pump-eV", "-1", "--fluence-mJcm2", "5", "--fwhm-fs", "50"]))
    with pytest.raises(SystemExit):
        rm.load_structure("no_such_material_xyz")
    from ase import Atoms
    bad = Atoms("C2", positions=[[0, 0, 0], [0.1, 0, 0]], cell=[5, 5, 5], pbc=True)
    with pytest.raises(SystemExit):
        rm.validate_structure(bad)


def test_semiconductor_carrier_tables_conserve_density():
    """Quasi-Fermi levels hold the photo-carrier density; the carrier energy grows with T_e."""
    bp, _, _ = rm.load_structure("black_phosphorus")
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    u = PhonopyAtoms(symbols=bp.get_chemical_symbols(), cell=np.array(bp.cell), scaled_positions=bp.get_scaled_positions())
    pr = Phonopy(u, supercell_matrix=[1, 1, 1], primitive_matrix="auto").primitive
    pdims = np.array([True] * 3)
    tb = rm.UniversalTB(pr.cell, pr.scaled_positions, pr.symbols, pr.masses, pdims)

    class A:
        nk_dos = [16, 16, 12]
    el = rm.electronic_setup(tb, A, pdims, 0.33)
    assert not el["metal"] and abs(el["ecbm"] - el["evbm"] - 0.33) < 1e-9
    el["nval"] = tb.nval
    tabs = rm.electron_tables(el, 0.1, 100.0)
    assert np.all(np.diff(tabs["E"]) > 0)
    nv = tb.nval // 2
    ec, g = rm.dos_hist(el["E"][:, nv:], el["ecbm"] - 0.01, el["ecbm"] + 6)
    T = tabs["T"][200]
    n = np.sum(g * rm.fermi((ec - tabs["muc"][200]) / (rm.KB_eV * T)))
    assert abs(n - 0.1) < 1e-4


def test_diffuse_structure_factor_independent_of_q_representative():
    """|F(Q)|^2 must not depend on which equivalent q (q or q + G') the eigenvectors are computed at,
    provided G = Q - q_used (phonopy phase convention)."""
    import h5py
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    a = 2.464
    cell = [[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, 15.0]]
    ph = Phonopy(PhonopyAtoms(symbols=["C", "C"], cell=cell, scaled_positions=[[1/3, 2/3, .5], [2/3, 1/3, .5]]),
                 supercell_matrix=[8, 8, 1], primitive_matrix="P")
    ph.force_constants = h5py.File(os.path.join(ROOT, "materials", "graphite", "fc2.hdf5"))["force_constants"][:]
    rec = 2 * np.pi * np.linalg.inv(np.array(cell)).T
    tau = np.array([[1/3, 2/3, .5], [2/3, 1/3, .5]]) @ np.array(cell)
    q1 = np.array([1/3, 1/3, 0.0]); q2 = q1 + np.array([-1, 0, 0])
    Q = (q1 + np.array([2, 1, 0])) @ rec
    out = []
    for qf in (q1, q2):
        f, ev = rm.phonons_at(ph, [qf])
        G = Q - qf @ rec
        I = []
        for j in range(3, 6):
            e = ev[0][:, j].reshape(2, 3)
            I.append(abs(np.sum((e @ Q) * np.exp(-1j * tau @ G))) ** 2)
        out.append(sorted(I))
    # the doubly degenerate pair at K has an arbitrary basis: compare the non-degenerate mode and the sum
    assert abs(out[0][0] - out[1][0]) < 1e-6 * out[0][0]
    assert abs(sum(out[0]) - sum(out[1])) < 1e-6 * sum(out[0])
