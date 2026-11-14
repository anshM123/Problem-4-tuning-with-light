"""Time-linear particle-hole fluctuation dynamics for local Hubbard interactions.

Spin-orbital basis a = (site, spin), M = 2L, index a = i (up) or i + L (down).
State:
  rho (M x M)      one-particle density matrix, rho_ab = <c^+_b c_a>
  V   (S x M x M)  Hermitian fluctuation samples of the one-body operators E_ab = c^+_b c_a,
                   normalised so that  (1/S) sum_s V_ab V_cd  estimates the symmetrized
                   covariance  1/2 <{dE_ab, dE_cd}>.
Equations (hbar = 1):
  i dV/dt   = [h^HF(rho), V] + [dh(V), rho]                      (linearized TDHF = RPA kernel)
  i drho/dt = [h^HF(rho), rho] + I,
  I_ab      = U ( G(ab; abar abar) - G(ab; bbar bbar) )          (exact collision term)
  G(ab;cd)  = (1/S) sum_s V_ab V_cd - S0(ab;cd)[rho]             (correlated part)
where abar is the opposite-spin partner on the same site and
  S0(ab;cd) = 1/2 [rho_cb rhobar_ad + rhobar_cb rho_ad].
Kernel variants (paper notation): 'rpa' (Hartree in h and dh), 'gw' (HF in h, Hartree in
dh), 'gwx' (HF in both).  Hopping is applied as a sparse matrix so the cost per step is
O(S M^2) for lattice models.
"""
from __future__ import annotations
import numpy as np
import scipy.sparse as sp


class PHFluct:
    def __init__(self, L, bonds, U, rho0, kernel="gwx", nsamp=None, seed=0,
                 Ufun=None, peierls=None):
        """bonds: list of (i, j, J_ij, dx) with hopping -J_ij c^+_i c_j (+h.c.), dx = x_i - x_j
        used for the Peierls phase exp(i A dx).  Ufun(t) overrides U (adiabatic switching).
        peierls(t) -> A(t) vector potential (scalar, along the bond displacement dx)."""
        self.L = L; self.M = M = 2 * L
        self.bonds = bonds; self.U0 = U
        self.Ufun = Ufun if Ufun is not None else (lambda t: U)
        self.Afun = peierls if peierls is not None else (lambda t: 0.0)
        self.kernel = kernel
        self.rho = rho0.astype(complex)
        self.partner = np.concatenate([np.arange(L, 2 * L), np.arange(L)])
        self.V = self._initial_samples(rho0, nsamp, seed)
        self.S = self.V.shape[0]

    # ----------------------------------------------------------------------
    def _initial_samples(self, rho0, nsamp, seed):
        """Hermitian samples with covariance S0[rho0].  Exact decomposition (deterministic
        sigma points over natural-orbital pairs) if nsamp is None, else Gaussian sampling."""
        M = self.M
        n, W = np.linalg.eigh(rho0)
        n = np.clip(n, 0, 1)
        var = 0.5 * (n[:, None] * (1 - n[None, :]) + n[None, :] * (1 - n[:, None]))
        if nsamp is None:
            samples = []
            for a in range(M):
                for c in range(a + 1, M):
                    if var[a, c] < 1e-14:
                        continue
                    w = np.sqrt(var[a, c] / 2)  # two samples: w^2 + w^2 = var
                    E = np.zeros((M, M), complex); E[a, c] = 1.0
                    samples.append(w * (E + E.T))
                    samples.append(1j * w * (E - E.T))
                if n[a] * (1 - n[a]) > 1e-14:
                    E = np.zeros((M, M), complex); E[a, a] = np.sqrt(n[a] * (1 - n[a]))
                    samples.append(E)
            V = np.array(samples) * np.sqrt(len(samples))  # (1/S) sum = plain sum
        else:
            rng = np.random.default_rng(seed)
            z = (rng.standard_normal((nsamp, M, M)) + 1j * rng.standard_normal((nsamp, M, M))) / np.sqrt(2)
            d = np.triu(np.sqrt(var)[None] * z, 1)
            d = d + np.conj(np.transpose(d, (0, 2, 1)))
            d[:, np.arange(M), np.arange(M)] = rng.standard_normal((nsamp, M)) * np.sqrt(n * (1 - n))[None]
            V = d
        # rotate from natural-orbital basis to site basis
        return np.einsum("ab,sbc,dc->sad", W, V, W.conj())

    # ----------------------------------------------------------------------
    def hop(self, t):
        A = self.Afun(t)
        if getattr(self, "_hopA", None) == A:
            return self._hopM
        L = self.L
        rows, cols, vals = [], [], []
        for (i, j, J, dx) in self.bonds:
            ph = np.exp(1j * A * dx)
            for s in (0, L):
                rows += [i + s, j + s]; cols += [j + s, i + s]
                vals += [-J * ph, np.conj(-J * ph)]
        self._hopA = A
        self._hopM = sp.csr_matrix((vals, (rows, cols)), shape=(self.M, self.M))
        self._hopD = self._hopM.toarray()
        return self._hopM

    def hf_terms(self, rho, U, exchange):
        """Interaction part of the HF Hamiltonian for density matrix rho (dense M x M)."""
        L = self.L
        H = np.zeros((self.M, self.M), complex)
        i = np.arange(L); j = i + L
        H[i, i] = U * rho[j, j]
        H[j, j] = U * rho[i, i]
        if exchange:
            H[i, j] = -U * rho[i, j]
            H[j, i] = -U * rho[j, i]
        return H

    def dh_batch(self, V, U, exchange):
        L = self.L
        dH = np.zeros_like(V)
        i = np.arange(L); j = i + L
        dH[:, i, i] = U * V[:, j, j]
        dH[:, j, j] = U * V[:, i, i]
        if exchange:
            dH[:, i, j] = -U * V[:, i, j]
            dH[:, j, i] = -U * V[:, j, i]
        return dH

    @staticmethod
    def _spmul_left(Hs, V):
        S, M, _ = V.shape
        out = Hs @ V.transpose(1, 0, 2).reshape(M, S * M)
        return out.reshape(M, S, M).transpose(1, 0, 2)

    def collision(self, rho, V, U):
        S = V.shape[0]
        M = self.M
        p = self.partner
        D = np.einsum("saa->sa", V)                  # V_qq per sample
        T1 = np.einsum("sab,sa->ab", V, D[:, p]) / S   # sum_s V_ab V_{abar abar}
        T2 = np.einsum("sab,sb->ab", V, D[:, p]) / S   # sum_s V_ab V_{bbar bbar}
        rb = np.eye(M) - rho
        a = np.arange(M)
        # S0(ab; qq) = 1/2 [rho_qb rhobar_aq + rhobar_qb rho_aq]
        q1 = p[:, None] * np.ones((1, M), int)      # q = abar
        q2 = np.ones((M, 1), int) * p[None, :]      # q = bbar
        A_, B_ = np.meshgrid(a, a, indexing="ij")
        S01 = 0.5 * (rho[q1, B_] * rb[A_, q1] + rb[q1, B_] * rho[A_, q1])
        S02 = 0.5 * (rho[q2, B_] * rb[A_, q2] + rb[q2, B_] * rho[A_, q2])
        return U * ((T1 - S01) - (T2 - S02))

    def rhs(self, t, rho, V):
        U = self.Ufun(t)
        hs = self.hop(t)
        exh = self.kernel in ("gw", "gwx")
        exd = self.kernel == "gwx"
        Hint = self.hf_terms(rho, U, exh)
        hd = self._hopD + Hint
        I = self.collision(rho, V, U) if self.collide else 0.0
        drho = -1j * (hd @ rho - rho @ hd + I)
        # samples
        HV = self._spmul_left(hs, V) + Hint[None] @ V
        VH = np.conj(np.transpose(self._spmul_left(hs, np.conj(np.transpose(V, (0, 2, 1)))), (0, 2, 1))) + V @ Hint[None]
        dH = self.dh_batch(V, U, exd)
        dV = -1j * (HV - VH + dH @ rho[None] - rho[None] @ dH)
        return drho, dV

    collide = True

    def step(self, t, dt):
        r, V = self.rho, self.V
        k1 = self.rhs(t, r, V)
        k2 = self.rhs(t + dt / 2, r + dt / 2 * k1[0], V + dt / 2 * k1[1])
        k3 = self.rhs(t + dt / 2, r + dt / 2 * k2[0], V + dt / 2 * k2[1])
        k4 = self.rhs(t + dt, r + dt * k3[0], V + dt * k3[1])
        self.rho = r + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        self.V = V + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])

    # observables ------------------------------------------------------------
    def cov_diag(self):
        """Symmetrized covariance of site-spin densities  <d n_a d n_b>  (M x M)."""
        D = np.real(np.einsum("saa->sa", self.V))
        return D.T @ D / self.S

    def szsz(self):
        """<dS^z_i dS^z_j> (L x L) from the sample covariance."""
        C = self.cov_diag(); L = self.L
        uu, dd, ud, du = C[:L, :L], C[L:, L:], C[:L, L:], C[L:, :L]
        return 0.25 * (uu + dd - ud - du)

    def double_occ(self):
        L = self.L
        n = np.real(np.diag(self.rho))
        C = self.cov_diag()
        i = np.arange(L)
        # <n_up n_dn> = n_up n_dn + cov (symmetrized = normal-ordered for commuting ops)
        return n[i] * n[i + L] + C[i, i + L]

    def current(self, t):
        """Total current j = -dH/dA = -2 sum_{bonds,spin} J dx Im(e^{iA dx} <c^+_a c_b>)."""
        A = self.Afun(t)
        L = self.L
        j = 0.0
        for (a, b, J, dx) in self.bonds:
            for s in (0, L):
                # <c^+_a c_b> = rho_{b a}
                j += -2 * J * dx * np.imag(np.exp(1j * A * dx) * self.rho[b + s, a + s])
        return j
