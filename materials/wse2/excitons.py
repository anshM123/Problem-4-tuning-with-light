"""Exciton binding energies in monolayer WSe2 on hBN from the Rytova-Keldysh interaction.

V(r) = -(e^2/4 pi eps0) (pi / 2 r0) [H0(eps r / r0) - Y0(eps r / r0)],
r0 = 2 pi chi_2D (DFT 2D polarizability, Berkelbach et al. PRB 88, 045318 (2013)),
eps = (eps_top + eps_bottom)/2 (vacuum / hBN).  The 2D radial Schroedinger equation is
solved on a grid (s states, m = 0).  Effective masses: DFT k.p values
(Kormanyos et al., 2D Mater. 2, 022001 (2015)).
"""
import numpy as np
from scipy.special import struve, y0
from scipy.linalg import eigh_tridiagonal

HB2M = 3.80998   # hbar^2 / (2 m0) in eV A^2
E2 = 14.39965    # e^2 / (4 pi eps0) in eV A

PARAMS = dict(
    r0=45.1,                  # A, WSe2 (Berkelbach 2013)
    eps_top=1.0,
    eps_bot=4.9,              # hBN static (sqrt(eps_par * eps_perp) ~ sqrt(6.9*3.5))
    me_K=0.29, mh_K=0.36,     # K-valley electron / hole (Kormanyos 2015)
    me_S=0.57,                # Sigma (Q) valley electron, geometric mean of 0.5 / 0.65
)


def binding(mu, r0, eps, n_states=1, R=400.0, N=6000):
    r = np.linspace(R / N, R, N)
    h = r[1] - r[0]
    x = eps * r / r0
    V = -E2 * (np.pi / (2 * r0)) * (struve(0, x) - y0(x))
    kin = HB2M / mu
    d = kin * (2 / h ** 2 - 1 / (4 * r ** 2)) + V
    off = -kin / h ** 2 * np.ones(N - 1)
    w = eigh_tridiagonal(d, off, select="i", select_range=(0, n_states - 1))[0]
    return -w


def exciton_table(p=PARAMS):
    eps = 0.5 * (p["eps_top"] + p["eps_bot"])
    muK = p["me_K"] * p["mh_K"] / (p["me_K"] + p["mh_K"])
    muS = p["me_S"] * p["mh_K"] / (p["me_S"] + p["mh_K"])
    return dict(Eb_KK=binding(muK, p["r0"], eps)[0], Eb_SK=binding(muS, p["r0"], eps)[0],
                M_KK=p["me_K"] + p["mh_K"], M_SK=p["me_S"] + p["mh_K"], eps=eps)


if __name__ == "__main__":
    t = exciton_table()
    print({k: round(v, 4) for k, v in t.items()})
    for eps_b in (1.0, 3.0, 4.9, 6.9):
        p = dict(PARAMS); p["eps_bot"] = eps_b
        tt = exciton_table(p)
        print(f"eps_bot={eps_b}: Eb(KK)={tt['Eb_KK']*1e3:.0f} meV  Eb(SK)={tt['Eb_SK']*1e3:.0f} meV  diff={1e3*(tt['Eb_SK']-tt['Eb_KK']):.0f}")
