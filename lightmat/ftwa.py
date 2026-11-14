"""Fermionic phase-space trajectories (Gaussian fluctuation sampling of the one-body
density matrix) for Hubbard-type models.

Every trajectory rho^(s) evolves under its own time-dependent Hartree-Fock Hamiltonian
    i d rho/dt = [h_eff(rho), rho],   (h_eff)_ab = dH_cl/d rho_ba,
and initial conditions are drawn from a Gaussian whose second moments reproduce the
symmetrized two-point fluctuations of the initial state:
    E[d rho_ab d rho_cd] = 1/2 <{dE_ab, dE_cd}>,   E_ab = c^+_b c_a.
Symmetric products of one-body operators are estimated by products of trajectory
variables, e.g. <n_iu n_id> = E[rho_{iu,iu} rho_{id,id}].

Two representations of the spinful Hubbard model are provided:
  basis='spin'   : orbitals (i up, i dn); Fock term = spin-flip exchange.
  basis='nambu'  : orbitals (c_i up, c^+_i dn); the down spin is particle-hole
                   transformed so the Fock term becomes the s-wave pairing field and
                   pair (T-matrix-like) fluctuations feed back on the densities.
The classical Hamiltonian is H = Tr(h rho) + U sum_i (rho_ii rho_jj - x rho_ij rho_ji),
j = partner orbital of site i, x = exchange weight.
"""
from __future__ import annotations
import numpy as np


class HubbardTWA:
    def __init__(self, h, U, n_up, n_dn, basis="nambu", exchange=1.0, nsamp=2000,
                 seed=0, hfun=None):
        self.L = L = h.shape[0]
        self.U, self.x, self.basis = U, exchange, basis
        self.hfun = hfun if hfun is not None else (lambda t: h)
        self.M = 2 * L
        rng = np.random.default_rng(seed)
        n_up = np.asarray(n_up, float); n_dn = np.asarray(n_dn, float)
        if basis == "nambu":
            occ = np.concatenate([n_up, 1.0 - n_dn])
        else:
            occ = np.concatenate([n_up, n_dn])
        rho0 = np.diag(occ).astype(complex)
        # Gaussian fluctuations in the (diagonal) natural-orbital basis
        M = self.M
        var = 0.5 * (occ[:, None] * (1 - occ[None, :]) + occ[None, :] * (1 - occ[:, None]))
        z = (rng.standard_normal((nsamp, M, M)) + 1j * rng.standard_normal((nsamp, M, M))) / np.sqrt(2)
        d = np.sqrt(var)[None] * z
        d = np.triu(d, 1)
        d = d + np.conj(np.transpose(d, (0, 2, 1)))
        diagvar = occ * (1 - occ)
        dd = rng.standard_normal((nsamp, M)) * np.sqrt(diagvar)[None]
        d[:, np.arange(M), np.arange(M)] = dd
        self.rho = rho0[None] + d
        self.idx = np.arange(L)

    def h_one(self, t):
        h = self.hfun(t)
        L = self.L
        H1 = np.zeros((self.M, self.M), complex)
        H1[:L, :L] = h
        if self.basis == "nambu":
            H1[L:, L:] = -h.T
        else:
            H1[L:, L:] = h
        return H1

    def heff(self, t, rho):
        L, U, x = self.L, self.U, self.x
        H = np.broadcast_to(self.h_one(t), rho.shape).copy()
        i = self.idx; j = self.idx + L
        rii = rho[:, i, i]; rjj = rho[:, j, j]
        rij = rho[:, i, j]; rji = rho[:, j, i]
        if self.basis == "nambu":
            # U n_up n_dn = U n_up - U n_up n'_2
            H[:, i, i] += U - U * rjj
            H[:, j, j] += -U * rii
            H[:, i, j] += x * U * rij
            H[:, j, i] += x * U * rji
        else:
            H[:, i, i] += U * rjj
            H[:, j, j] += U * rii
            H[:, i, j] += -x * U * rij
            H[:, j, i] += -x * U * rji
        return H

    def rhs(self, t, rho):
        H = self.heff(t, rho)
        return -1j * (H @ rho - rho @ H)

    def step(self, t, dt):
        r = self.rho
        k1 = self.rhs(t, r)
        k2 = self.rhs(t + dt / 2, r + dt / 2 * k1)
        k3 = self.rhs(t + dt / 2, r + dt / 2 * k2)
        k4 = self.rhs(t + dt, r + dt * k3)
        self.rho = r + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    # observables -----------------------------------------------------------
    def nup_ndn(self):
        L = self.L
        a = np.real(self.rho[:, self.idx, self.idx]).mean(0)
        b = np.real(self.rho[:, self.idx + L, self.idx + L]).mean(0)
        return (a, 1 - b) if self.basis == "nambu" else (a, b)

    def double_occ(self):
        L = self.L
        a = np.real(self.rho[:, self.idx, self.idx])
        b = np.real(self.rho[:, self.idx + L, self.idx + L])
        if self.basis == "nambu":
            b = 1 - b
        return (a * b).mean(0)
