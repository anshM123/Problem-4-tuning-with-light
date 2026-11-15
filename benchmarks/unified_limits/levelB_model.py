"""Model for the Level B (kinetic) limit: a band of electron levels coupled to one optical
phonon (random e-ph matrix elements, no diagonal coupling), which decays anharmonically into
pairs of a band of acoustic modes -- graphite's e -> optical -> acoustic cascade in miniature.

Also contains the closed-form Boltzmann (golden-rule) collision integrals used as a check,
in the conventions of materials/graphite/dynamics.py.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from lightmat.unified import Model, UnifiedSolver


def build(L=40, M=16, W=4.0, w0=1.0, wac=(0.2, 0.8), g=0.05, k0=0.1, seed=1, jitter=0.0):
    rng = np.random.default_rng(seed)
    eps = np.linspace(-W / 2, W / 2, L)
    if jitter:      # incommensurate levels (avoids degenerate pair energies and early recurrences)
        d = W / (L - 1)
        eps = eps + jitter * d * rng.uniform(-0.5, 0.5, L)
        wj = np.linspace(wac[0], wac[1], M)
        wj = wj + jitter * (wac[1] - wac[0]) / max(M - 1, 1) * rng.uniform(-0.5, 0.5, M)
    Q = 1 + M
    omega = np.concatenate([[w0], wj if jitter else np.linspace(wac[0], wac[1], M)])
    gam = np.zeros((Q, L, L), complex)
    r = rng.standard_normal((L, L)) + 1j * rng.standard_normal((L, L))
    r = (r + r.conj().T) / 2
    np.fill_diagonal(r, 0.0)
    gam[0] = g * r
    kap = np.zeros((Q, Q, Q))
    rk = rng.standard_normal((M, M)); rk = (rk + rk.T) / np.sqrt(2)
    for b in range(M):
        for c in range(M):
            v = k0 * rk[b, c]
            for p in ((0, 1 + b, 1 + c), (1 + b, 0, 1 + c), (1 + b, 1 + c, 0)):
                kap[p] = v
    m = Model(np.diag(eps), U=0.0, omega=omega, gamma=gam, kappa=kap)
    return m, eps


def gauss(x, s):
    return np.exp(-0.5 * (x / s) ** 2) / (np.sqrt(2 * np.pi) * s)


def boltzmann(m, eps, f, N, s):
    """Golden-rule collision integrals (spin-degenerate electrons, f per spin)."""
    Q = m.Q
    df = np.zeros_like(f); dN = np.zeros(Q)
    for a in range(Q):
        g2 = np.abs(m.gamma[a]) ** 2
        if not g2.any():
            continue
        de = eps[:, None] - eps[None, :]                     # e_i - e_j
        em = g2 * gauss(de - m.omega[a], s) * (f[:, None] * (1 - f[None, :]) * (N[a] + 1)
                                               - f[None, :] * (1 - f[:, None]) * N[a])   # i -> j + phonon
        dN[a] += np.pi * 2 * em.sum()
        df += -np.pi * em.sum(1) + np.pi * em.sum(0)
    k2 = m.kappa ** 2
    w = m.omega
    for a in range(Q):
        dec = gauss(w[a] - w[:, None] - w[None, :], s) * (N[a] * (1 + N[:, None]) * (1 + N[None, :])
                                                         - (1 + N[a]) * N[:, None] * N[None, :])
        coa = gauss(w[a] + w[:, None] - w[None, :], s) * (N[a] * N[:, None] * (1 + N[None, :])
                                                         - (1 + N[a]) * (1 + N[:, None]) * N[None, :])
        dN[a] -= np.pi / 8 * np.sum(k2[a] * (dec + 2 * coa))
    return df, dN


def nambu_occ(f):
    """Nambu occupations (up particles, down holes) for spin-degenerate f."""
    return np.concatenate([f, 1 - f])
