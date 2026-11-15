"""Consistency tests of the WSe2 exciton kinetics (run: python -m pytest tests/test_wse2.py -q)."""
import os, sys
import numpy as np
HERE = os.path.join(os.path.dirname(__file__), "..", "materials", "wse2")
sys.path.insert(0, HERE)
from excitons import exciton_table
from dynamics import PARAMS, build, KB
import bloch as B


def _params(**kw):
    tab = exciton_table(); P = dict(PARAMS)
    P["dE_SK"] = P["dCB_SK"] - (tab["Eb_SK"] - tab["Eb_KK"])
    P["M_KK"], P["M_SK"] = tab["M_KK"], tab["M_SK"]
    P.update(kw)
    return P


def test_generator_conserves_particles():
    L = B.generator(build(_params())).toarray()
    assert np.max(abs(L.sum(0))) < 1e-6 * np.max(abs(L))


def test_generator_relaxes_to_boltzmann_valley_ratio():
    P = _params(); m = build(P)
    md = B.relaxation_modes(P)
    kT = KB * P["T"]
    r_boltz = (4 * m["M"]["KK"]) / (12 * m["M"]["SK"] * np.exp(-P["dE_SK"] / kT))
    assert abs(md["r_eq"] / r_boltz - 1) < 0.05


def test_measurement_operator_on_data():
    me = B.as_measured(B.EXP_T, B.EXP_K, B.EXP_S)
    assert me["t_Kmax"] == 0.09 and me["t_Smax"] == 0.4


def test_two_level_fit_recovers_time():
    tp = B.EXP_T
    k, s = B.two_level(tp, 10.0, 0.2, 0.5, 2.5)
    f = B.fit_two_level(tp, k, s)
    assert abs(f["tau_f"] - 0.2) < 0.01


def test_coherent_stage_absorbs_stated_density():
    P = _params(g_EEA=0.0)
    r = B.run_bloch(P, tmax=1e-12, nt=101)
    tot = r["K"][-1] + r["S"][-1]
    assert abs(tot / P["n_exc"] - 1) < 0.02
