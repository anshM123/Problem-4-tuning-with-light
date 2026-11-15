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
                 Ufun=None, peierls=None, onsite=None):
        """bonds: list of (i, j, J_ij, dx) with hopping -J_ij c^+_i c_j (+h.c.), dx = x_i - x_j
        used for the Peierls phase exp(i A dx).  Ufun(t) overrides U (adiabatic switching).
        peierls(t) -> A(t) vector potential (scalar, along the bond displacement dx)."""
        self.L = L; self.M = M = 2 * L
        self.bonds = bonds; self.U0 = U
        self.Ufun = Ufun if Ufun is not None else (lambda t: U)
        self.Afun = peierls if peierls is not None else (lambda t: 0.0)
        self.kernel = kernel
        self.onsite = onsite
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
        p = self.partner
        a = np.arange(self.M)
        # local (on-site) HF potential: diagonal d and spin-partner element o
        d = U * np.real(rho[p, p]).astype(complex)
        o = (-U * rho[a, p]) if exh else np.zeros(self.M, complex)
        if self.onsite is not None:
            d = d + self.onsite(t)
        Hint = np.zeros((self.M, self.M), complex)
        Hint[a, a] = d; Hint[a, p] = o
        hd = self._hopD + Hint
        I = self.collision(rho, V, U) if self.collide else 0.0
        drho = -1j * (hd @ rho - rho @ hd + I)
        # samples: sparse hopping + local terms, O(S M^2)
        HV = self._spmul_left(hs, V) + d[None, :, None] * V + o[None, :, None] * V[:, p, :]
        VH = (np.conj(np.transpose(self._spmul_left(hs, np.conj(np.transpose(V, (0, 2, 1)))), (0, 2, 1)))
              + V * d[None, None, :] + V[:, :, p] * o[p][None, None, :])
        dd = U * V[:, p, p]
        od = (-U * V[:, a, p]) if exd else np.zeros_like(dd)
        dHr = dd[:, :, None] * rho[None, :, :] + od[:, :, None] * rho[None, p, :]
        rdH = rho[None, :, :] * dd[:, None, :] + rho[None][:, :, p] * od[:, p][:, None, :]
        dV = -1j * (HV - VH + dHr - rdH)
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


# ==========================================================================================
# Doublon-holon pair closure of the particle-hole channel (strong coupling, 1D chains)
# ==========================================================================================
def hf_antiferromagnet(h, U, nel=None, tol=1e-12, maxit=20000):
    """Self-consistent collinear Hartree-Fock antiferromagnet of a bipartite chain with hopping
    matrix h at half filling.  Returns rho_up, rho_dn (L x L) and the HF energy."""
    L = h.shape[0]
    N = L // 2 if nel is None else nel
    stag = np.array([(-1) ** i for i in range(L)], float)
    nu = 0.5 + 0.4 * stag; nd = 1 - nu
    for _ in range(maxit):
        _, vu = np.linalg.eigh(h + U * np.diag(nd)); _, vd = np.linalg.eigh(h + U * np.diag(nu))
        ru = vu[:, :N] @ vu[:, :N].T; rd = vd[:, :N] @ vd[:, :N].T
        nu2, nd2 = np.diag(ru).copy(), np.diag(rd).copy()
        if max(abs(nu2 - nu).max(), abs(nd2 - nd).max()) < tol:
            break
        nu, nd = 0.5 * nu + 0.5 * nu2, 0.5 * nd + 0.5 * nd2
    E = float(np.trace(h @ (ru + rd)) + U * np.sum(np.diag(ru) * np.diag(rd)))
    return ru, rd, E


class DoublonHolonPairs:
    """Optical response of the half-filled open Hubbard chain from doublon-holon pair amplitudes.

    At half filling the current operator acting on a Mott ground state creates a doublon (D)
    and a holon (H) on neighbouring sites.  The particle-hole amplitude of such an
    excitation is propagated *exactly as a two-body problem* in the self-interaction-free
    field of the other electrons, the counterpart in the particle-hole channel of the
    self-interaction-free pp ladder (fluct_hubbard.PPLadder, mode 'sic'):

    * the other electrons form a singly occupied background; in one dimension a doublon or
      holon moving through it only shifts the order of the background spins (spin-charge
      separation of the open chain, Ogata and Shiba), so both move with the full hopping J.
      With the operator ordering site by site (up before down) the doublon hop has amplitude
      +J and the holon hop -J; D and H cannot occupy the same site;
    * second order in J/U (Schrieffer-Wolff, antiferromagnetically ordered background):
      - a background bond that touches the doublon or the holon can no longer gain its
        exchange energy, so the pair energy is  E(i, j) = U - sum_{b touching i or j} e_b,
        with e_b the ground-state energy of bond b (its kinetic energy K_b plus the share
        U <n_up n_dn> / z of each of its two sites, z = number of bonds at the site);
      - three-site processes move D or H by two sites across a singly occupied site with
        amplitude t2 = -2 P_s J^2/U = K_b / 4, P_s being the singlet weight of the bond b
        that is crossed (K_b = -8 J^2 P_s / U is its kinetic energy);
      for a uniform ground state the band edge is U - 4J - 4 e_b - |K_b|, and with the exact
      bond energy and singlet weight (e_b = -4 ln2 J^2/U, P_s = ln 2) this is the exact
      strong-coupling Mott gap U - 4J + 8 ln2 J^2/U;
    * the current creates the pair from bond (k, k+1) with amplitude proportional to the
      square root of the bond singlet weight, i.e. of -K_k; for the antiferromagnetic order of
      the (squeezed) background spin chain the amplitudes alternate in sign,
      j_+|s> ~ sum_k (-1)^k sqrt(-K_k) (|H_k D_k+1> - |D_k H_k+1>).  At the same (first) order
      in J/U as the weights, the virtual pairs of the ground state |0> = |s> - K_+|s>/U add the
      vertex correction -(1/U) j_11 K_+|s>, in which the current moves the D or H of a virtual
      pair (it creates pairs at distance 2 that interfere with the direct ones).  The overall
      scale is fixed by the exact f-sum rule  2 <j|(H - E_0)^-1|j> = -<K>.

    All inputs (bond kinetic energies K_b and bond energies e_b) come from the method's own
    ground state, by default the Hartree-Fock antiferromagnet.  The linear response is obtained
    in real time: after a vector-potential step A0 theta(t) the induced current is
        delta j(t) / A0 = <K> + int_0^t chi(s) ds = -2 Re <xi| exp(-i (H - E_0) t) |j>,
    xi = (H - E_0)^-1 j|0>, and Re sigma(w) is its Fourier transform with exp(-eta t), as for
    the other fast methods.  The cost is O(L^2) per time step.
    periodic=True closes the chain into a ring (kbond, ebond of length L; the cyclic shift of the
    background spins when a particle crosses the boundary is neglected, an O(1/L) effect), used
    for long-chain spectra without edge states.
    order=1 keeps only the leading order for comparison: pair energy U - eps_i - eps_j with the
    site energies eps (two sites of correlation energy removed), no three-site hops and no
    vertex correction."""

    def __init__(self, L, U, kbond, ebond, J=1.0, order=2, eps=None, periodic=False):
        import scipy.sparse as sp
        self.L, self.U = L, U
        nb = L if periodic else L - 1                      # bonds b = (b, b+1 mod L)
        site = (lambda x: x % L) if periodic else (lambda x: x if 0 <= x < L else None)
        kbond = np.asarray(kbond, float); ebond = np.asarray(ebond, float)
        idx = -np.ones((L, L), int)
        pairs = [(i, j) for i in range(L) for j in range(L) if i != j]
        for k, (i, j) in enumerate(pairs):
            idx[i, j] = k
        n = len(pairs)
        touch = lambda i: {b % L for b in (i - 1, i) if periodic or 0 <= b < L - 1}
        rows, cols, vals = [], [], []
        for k, (i, j) in enumerate(pairs):
            if order == 1:
                e = U - eps[i] - eps[j]
            else:
                e = U - sum(ebond[b] for b in touch(i) | touch(j))
            rows.append(k); cols.append(k); vals.append(e)
            for d in (-1, 1):
                i1, j1 = site(i + d), site(j + d)
                if i1 is not None and i1 != j:                  # doublon hop
                    rows.append(idx[i1, j]); cols.append(k); vals.append(+J)
                if j1 is not None and j1 != i:                  # holon hop
                    rows.append(idx[i, j1]); cols.append(k); vals.append(-J)
                if order == 1:
                    continue
                # three-site hops of D or H across the singly occupied site a+d to c = a+2d
                for a, other, isD in ((i, j, True), (j, i, False)):
                    m1, c = site(a + d), site(a + 2 * d)
                    if c is not None and m1 != other and c != other:
                        t2 = kbond[m1 if d == 1 else c] / 4.0  # bond between m1 and c
                        rows.append(idx[c, other] if isD else idx[other, c]); cols.append(k)
                        vals.append(t2)
        self.H = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
        # current vertex: j|0> (common factor i J sqrt2 dropped) and, for order 2, its
        # first-order correction from the virtual pairs of the ground state,
        #   j|0> = j_+ |s> - (1/U) j_11 K_+ |s>,   |0> = |s> - K_+|s>/U,
        # where K_+ creates a pair on a bond (even combination), and j_11 = -i J sum_d d T_d moves
        # D or H by d within the one-pair sector.
        phi = np.zeros(n); psiK = np.zeros(n)
        for k in range(nb):
            a = np.sqrt(max(-kbond[k], 0.0)) * (-1) ** k
            phi[idx[(k + 1) % L, k]] += a; phi[idx[k, (k + 1) % L]] -= a
            psiK[idx[(k + 1) % L, k]] -= a; psiK[idx[k, (k + 1) % L]] -= a
        if order >= 2:
            corr = np.zeros(n)
            for k, (i, j) in enumerate(pairs):
                if psiK[k] == 0.0:
                    continue
                for d in (-1, 1):
                    i1, j1 = site(i + d), site(j + d)
                    if i1 is not None and i1 != j:
                        corr[idx[i1, j]] -= d * psiK[k]
                    if j1 is not None and j1 != i:
                        corr[idx[i, j1]] -= d * psiK[k]
            phi = phi - (J / U) * corr
        from scipy.sparse.linalg import spsolve
        xi = spsolve(self.H.tocsc(), phi)
        Ktot = float(np.sum(kbond))
        scale = np.sqrt(-Ktot / (2 * phi @ xi))           # exact f-sum rule
        self.phi, self.xi = phi * scale, xi * scale
        self.n = n

    def current_response(self, T=60.0, dt=0.02):
        """delta j(t)/A0 on the grid t = 0, dt, ..., by RK4 propagation of the pair amplitudes."""
        H = self.H
        f = lambda v: -1j * (H @ v)
        psi = self.phi.astype(complex)
        nt = int(round(T / dt))
        out = np.empty(nt)
        for k in range(nt):
            out[k] = -2.0 * np.real(np.vdot(self.xi, psi))
            k1 = f(psi); k2 = f(psi + dt / 2 * k1); k3 = f(psi + dt / 2 * k2); k4 = f(psi + dt * k3)
            psi = psi + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return np.arange(nt) * dt, out


def local_energies(h, U, ru, rd, d=None):
    """Energy decomposition of a chain ground state (rho_up, rho_dn, double occupation d; the
    Hartree-Fock value n_up n_dn if d is None): bond kinetic energies K_b (bond b = (b, b+1)),
    bond energies e_b = K_b + sum_{sites s of b} U d_s / z_s (z_s = number of bonds at s) and
    site energies eps_s = sum_{b at s} K_b / 2 + U d_s;  sum e_b = sum eps_s = E_0."""
    L = h.shape[0]
    Kb = np.array([2 * h[k, k + 1] * np.real(ru[k, k + 1] + rd[k, k + 1]) for k in range(L - 1)])
    if d is None:
        d = np.real(np.diag(ru) * np.diag(rd))
    d = np.asarray(d, float)
    z = np.full(L, 2.0); z[0] = z[-1] = 1.0
    eb = Kb + U * (d[:-1] / z[:-1] + d[1:] / z[1:])
    eps = U * d.copy()
    eps[:-1] += Kb / 2; eps[1:] += Kb / 2
    return Kb, eb, eps
