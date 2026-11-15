"""One time-linear solver for electrons, phonons and light: the G1-G2 hierarchy for
coupled electron-phonon matter, and the three limits used in this repository.

Hamiltonian (hbar = 1; spin-orbitals i sigma, phonon modes a = 1..Q, x = (b + b^+)/sqrt2,
p = (b - b^+)/(i sqrt2)):

  H = sum_{ij,s} h_ij(t) c^+_is c_js                                   one body, light (Peierls/dipole)
    + U sum_i n_iu n_id + sum_ij V_ij c^+_iu c^+_id c_jd c_ju            e-e (Hubbard; reduced pairing V)
    + sum_a w_a (x_a^2 + p_a^2)/2 + sum_a f_a(t) x_a                    phonons, dipole force
    + sum_a x_a sum_{ij,s} gamma^a_ij c^+_is c_js                        e-ph, linear (gamma = sqrt2 g)
    + (1/6) sum_abc kappa_abc x_a x_b x_c                                ph-ph, cubic anharmonicity

Electrons are described in Nambu form, Psi = (c_1u..c_Lu, c^+_1d..c^+_Ld), with the
generalized density matrix Gam_IJ = <Psi^+_J Psi_I> = [[rho_u, kappa], [kappa^+, 1 - rho_d^T]]
(kappa_ij = <c_jd c_iu>), so normal and anomalous (pairing) averages are one object.

Propagated objects
  Gam                       electron one-body (Nambu) density matrix
  x, y                      factors of the up-down pair matrices X = A + G2, Y = B + G2 (Level A;
                            G2 = connected two-particle correlation, kept positive, lightmat/fluct_hubbard)
  Phi = (<x_a>, <p_a>)      coherent phonon displacements
  S_mn = <{dxi_m, dxi_n}>/2 phonon one-body density matrices (<b^+b>, <bb> as a covariance; vacuum S = 1/2)
  T_mnl                     connected three-phonon correlation (symmetrized third cumulant)
  C^m_IJ = <Psi^+_J Psi_I dxi_m>  connected electron-phonon correlation (2Q Hermitian matrices)
  W                         work done by the time-dependent fields (energy bookkeeping)

Equations of motion.  They follow from the Heisenberg equations (equivalently the Wigner-Moyal
equation for the phonon coordinates, which for linear e-ph coupling and a cubic potential has
only one quantum correction, -(1/24) kappa d^3W/dp^3) by dropping connected correlations with
more than three field operators, except the e-e pair correlation G2, which is kept in the
positive-factorized pp-ladder form of Level A.  With H_N the Hartree-Fock-Bogoliubov Nambu
Hamiltonian (incl. sum_a X_a G^a, G^a = diag(gamma^a, -gamma^aT)) and A_eff the linear phonon
matrix (dx_a = w_a p_a, dp_a = -sum_c (w_a d_ac + kappa_abc X_b) x_c):

 (1) i dGam/dt = [H_N, Gam] + sum_a [G^a, C^{x_a}] + I_ee[G2]                     (exact BBGKY 1)
 (2) dX_a = w_a P_a,  dP_a = -w_a X_a - f_a - Tr(G^a Gam) - Tr gamma^a
                         - 1/2 kappa_abc (X_b X_c + S_bc)                             (exact)
 (3) dS = A S + S A^T + R + R^T,  R_{p_a,n} = -Tr(G^a C^n) - 1/2 kappa_abc T_{b c n}  (exact)
 (4) dC^m = -i([H_N, C^m] + [dH_N[C^m], Gam] + sum_b S_{x_b m} [G^b, Gam]) + sum_n A_mn C^n
            - d_{m,p_a} { 1/2 [(1-Gam) G^a Gam + Gam G^a (1-Gam)] + gamma^a.G2 }     (e-ph source)
 (5) dT_mnl = F_mnl + F_nml + F_lnm + d_{m,p_a} d_{n,p_b} d_{l,p_c} kappa_abc/4,
            F_mnl = sum_k A_mk T_knl - d_{m,p_a} kappa_abc S_{b n} S_{c l}               (ph-ph source)
 (6) G2: i dx/dt = H2 x,  i dy/dt = H2 y,  H2 = h2[H_N] + U (A - B) P   (Level A, h2 incl. the
            phonon field; 'factorized' = T-matrix ladder, 'sic' = self-interaction-free ladder)
 (6') G2 += G2ph:  i dG2ph/dt = H2 G2ph - G2ph H2^+ + sum_e ([gamma^e, rho_u] (x) C_d^e + C_u^e (x) [gamma^e, rho_d])
            the phonon-exchange part of the up-down pair correlation (the retarded phonon-mediated
            interaction acting in the pair channel), kept as a dense L^2 x L^2 matrix; it vanishes
            identically without phonons.

The source in (4) is the Fan-Migdal source of the GD approximation (Karlsson, van Leeuwen,
Perfetto, Stefanucci, PRL 127, 036402 (2021)); with S = (N + 1/2) it contains the
spontaneous and stimulated factors (1-f) f (1+N) - f (1-f) N.  The constant kappa/4 in (5) is
the quantum (hbar^2) term: it produces spontaneous three-phonon decay and is what a
classical-trajectory treatment of the phonons misses (zero-point leakage).  The gamma.G2 term
lets the e-e correlation dress the e-ph vertex.  Energy E = <H> is conserved up to the work W
of the external fields (checked in tests/test_unified.py).

Dropped: connected correlations with four or more field operators except G2 (in particular
<c+ c+ c c dxi>_c and <c+ c dxi dxi>_c and the fourth phonon cumulant).  Consequences
(benchmarks/unified_limits): without e-e interaction the scheme is exact to O(g^2) and all
observables of the driven Holstein dimer have O(g^4) errors; with strong e-e correlation C is
propagated with the TDHF/RPA electron dynamics, so the incoherent phonon emission is only
approximate; the second-order closures are not positivity preserving and become unstable when
relaxation rates approach the phonon frequencies (bipolaron regime, strong three-phonon mixing).

Limits (the three levels of REPORT.md section 2):
 * Level A: gamma = kappa = 0, Q = 0.  (2)-(5) are absent, (1) and (6) are exactly the
   equations of lightmat/fluct_hubbard.PPLadder ('factorized' or 'sic').  UnifiedSolver
   calls the same pair-ladder kernels; the two codes agree to rounding (benchmarks/unified_limits).
 * Level B: Markov + random-phase limit.  (4) and (5) are linear in C, T with sources fixed by
   Gam, S; solving them adiabatically, C = int_0^inf exp(L tau) m(tau) dtau * source (memory
   kernel m = exp(-eta tau) or exp(-sigma^2 tau^2/2)), keeping only populations (Gam diagonal
   in the mean-field eigenbasis, S diagonal in normal modes) and inserting into (1), (3) gives
   the Boltzmann equations: electron-phonon golden rule with
   f_i (1-f_j)(1+N) - f_j (1-f_i) N, and three-phonon decay/coalescence with
   (1+N_a) N_b N_c - N_a (1+N_b)(1+N_c) at rate (pi/8)|kappa_abc|^2 delta(w_a - w_b - w_c) per
   ordered pair.  This is what materials/graphite/dynamics.py integrates (method
   kinetic_rhs; delta width = the memory-kernel width).
 * Level C: V < 0 (phonon-mediated pairing, Eliashberg-derived), no G2.  The stationary point of
   (1) is the BCS/two-band gap equation Delta_a = -sum_b V_ab sum_{k in b} kappa_k with
   kappa_k = -Delta/(2E) tanh(E/2T) (method equilibrium); the linearization of (1) around it
   (time-dependent BdG/RPA) gives the amplitude propagator
   D(w) = [V^-1 - diag(N_a P_a(w))]^-1 of materials/mgb2/twoband.py, with the Higgs mode at 2 Delta.

Material scale (end of this file; used by materials/graphite_unified): coherent Bloch density
matrices over a full Brillouin zone driven by the light field, with the e-e, e-ph and ph-ph
correlations in their Markov / random-phase limit (energy-resolved collision operators that
conserve particle number and energy; tested against kinetic_rhs in tests/test_unified.py).
"""
from __future__ import annotations
import numpy as np
from scipy.special import dawsn

from .fluct_hubbard import PPLadder


# --------------------------------------------------------------------------------------
class Model:
    """Parameters of H.  h: (L,L) array or callable t -> (L,L) (spin independent, Hermitian);
    U: Hubbard U; V: (L,L) Hermitian reduced pairing interaction (pairing channel only) or None;
    pair_field: callable t -> (L,L) external pair field F_ij(t) (adds F c^+_iu c^+_jd + h.c.);
    omega: (Q,) phonon frequencies; gamma: (Q,L,L) e-ph couplings in x units (Hermitian);
    kappa: (Q,Q,Q) symmetric cubic constants; force: callable t -> (Q,) dipole force on x."""

    def __init__(self, h, U=0.0, V=None, pair_field=None, omega=None, gamma=None, kappa=None,
                 force=None):
        self.hfun = h if callable(h) else (lambda t, h=np.asarray(h, complex): h)
        self.L = L = self.hfun(0.0).shape[0]
        self.U = float(U)
        self.V = None if V is None else np.asarray(V, complex)
        self.pair_field = pair_field
        self.omega = np.zeros(0) if omega is None else np.asarray(omega, float)
        self.Q = Q = len(self.omega)
        self.gamma = np.zeros((Q, L, L), complex) if gamma is None else np.asarray(gamma, complex)
        self.kappa = np.zeros((Q, Q, Q)) if kappa is None else np.asarray(kappa, float)
        self.force = force
        if Q:
            k = self.kappa
            assert np.allclose(k, np.transpose(k, (1, 0, 2))) and np.allclose(k, np.transpose(k, (0, 2, 1))), \
                "kappa must be fully symmetric"

    def h(self, t):
        return self.hfun(t)

    def f(self, t):
        return np.zeros(self.Q) if self.force is None else np.asarray(self.force(t), float)


def peierls_chain(L, J, A=None, periodic=False):
    """h(t)_ij for a chain with hopping -J exp(i A(t) (x_i - x_j)) (A in units of 1/lattice)."""
    bonds = [(i, i + 1) for i in range(L - 1)] + ([(L - 1, 0)] if periodic else [])

    def h(t):
        ph = 1.0 if A is None else np.exp(1j * A(t))
        m = np.zeros((L, L), complex)
        for i, j in bonds:
            m[i, j] += -J * ph
            m[j, i] += -J * np.conj(ph)
        return m
    return h


# --------------------------------------------------------------------------------------
class UnifiedSolver:
    """ee: None | 'hf' (mean field only) | 'factorized' (T-matrix pp ladder) | 'sic'
          (self-interaction-free pp ladder); the pp ladder needs kappa = 0 (normal state).
    eph, phph: 'corr' (propagate C, T) | 'meanfield' (Ehrenfest: C = T = 0).
    positivity: after every step, restore the representability of the phonon and mixed
          correlations: S + i Omega/2 >= 0 (uncertainty relation) and positivity of the joint
          covariance of the e-ph coupling operators O_a = sum gamma^a c+c and the phonon quadratures,
          by the minimal uniform rescaling of C (the analogue of Level A's positive factorization;
          not energy conserving, the drift is reported).
    sum_rules (sic only): 'energy' (pair sum rules + energy, as Level A) | 'pairs' | None |
          'auto' (= 'energy' without phonons, i.e. Level A; None with phonons, where eqs. (1)-(6)
          conserve the energy by themselves and the rescaling would act on the pair amplitudes only)."""

    def __init__(self, model, rho_up, rho_dn, kappa=None, X0=None, P0=None, nph0=None,
                 ee="sic", eph="corr", phph="corr", vertex_ee=True, sum_rules="auto", g2ph=None,
                 positivity=False):
        m = self.m = model
        L, Q = m.L, m.Q
        self.L, self.Q = L, Q
        G = np.zeros((2 * L, 2 * L), complex)
        G[:L, :L] = rho_up
        G[L:, L:] = np.eye(L) - np.asarray(rho_dn).T
        if kappa is not None:
            G[:L, L:] = kappa
            G[L:, :L] = np.conj(kappa).T
        self.G = G
        self.eph_corr = (eph == "corr") and Q > 0
        self.phph_corr = (phph == "corr") and Q > 0
        self.vertex_ee = vertex_ee
        # Nambu e-ph vertices
        self.GN = np.zeros((Q, 2 * L, 2 * L), complex)
        for a in range(Q):
            self.GN[a, :L, :L] = m.gamma[a]
            self.GN[a, L:, L:] = -m.gamma[a].T
        self.trg = np.real(np.einsum("aii->a", m.gamma))
        # phonons
        self.Phi = np.zeros(2 * Q)
        if X0 is not None:
            self.Phi[:Q] = X0
        if P0 is not None:
            self.Phi[Q:] = P0
        n0 = np.zeros(Q) if nph0 is None else np.asarray(nph0, float) * np.ones(Q)
        self.S = np.diag(np.concatenate([n0 + 0.5, n0 + 0.5]))
        self.T = np.zeros((2 * Q,) * 3)
        self.C = np.zeros((2 * Q, 2 * L, 2 * L), complex)
        self.W = 0.0
        # e-e correlations (Level A kernels)
        self.ee = ee
        self.pp = None
        if ee in ("factorized", "sic") and m.U != 0.0:
            if kappa is not None and np.abs(kappa).max() > 0:
                raise ValueError("the pp ladder is formulated for the normal state (kappa = 0)")
            self.pp = PPLadder(m.h(0.0), m.U, np.asarray(rho_up, complex), np.asarray(rho_dn, complex),
                               mode=ee)
            self.pp._Ucur = m.U
            self.pp.hfun = m.hfun
            self.x, self.y = self.pp.x, self.pp.y
        else:
            self.x = np.zeros((L * L, 0), complex)
            self.y = np.zeros((L * L, 0), complex)
        if sum_rules == "auto":
            sum_rules = "energy" if Q == 0 else None
        self.sum_rules = sum_rules if ee == "sic" else None
        # phonon-exchange part of the up-down pair correlation (eq. 6'), full L^2 x L^2 matrix
        if g2ph is None:
            g2ph = self.eph_corr and ee in ("hf", "factorized", "sic")
        self.g2ph = bool(g2ph)
        self.Gph = np.zeros((L * L, L * L) if self.g2ph else (0, 0), complex)
        self.Nup = float(np.real(np.trace(rho_up))); self.Ndn = float(np.real(np.trace(rho_dn)))
        self.positivity = positivity
        self.t = 0.0
        self.E0 = None

    # ---------------------------------------------------------------- helpers
    def rhos(self, G=None):
        G = self.G if G is None else G
        L = self.L
        return G[:L, :L], np.eye(L) - G[L:, L:].T

    def hph(self, t, X):
        h = self.m.h(t)
        if self.Q:
            h = h + np.einsum("a,aij->ij", X, self.m.gamma)
        return h

    def hfb(self, t, G, X):
        """Nambu mean-field Hamiltonian H_N[Gam] incl. the coherent phonon field."""
        m, L = self.m, self.L
        h = self.hph(t, X)
        nu = np.real(np.diag(G[:L, :L])); nd = 1.0 - np.real(np.diag(G[L:, L:]))
        H = np.zeros((2 * L, 2 * L), complex)
        H[:L, :L] = h + m.U * np.diag(nd)
        H[L:, L:] = -(h + m.U * np.diag(nu)).T
        D = self._pairing(G)
        if m.pair_field is not None:
            D = D + m.pair_field(t)
        H[:L, L:] = D
        H[L:, :L] = np.conj(D).T
        return H

    def _pairing(self, M, lower=False):
        """Pairing field Delta_ij = d_ij (U kappa_ii + sum_k V_ik kappa_kk) from a Nambu matrix."""
        m, L = self.m, self.L
        k = np.diag(M[L:, :L]) if lower else np.diag(M[:L, L:])
        if m.V is None and m.U == 0.0:
            return np.zeros((L, L), complex)
        d = m.U * k
        if m.V is not None:
            d = d + ((np.conj(m.V) @ k) if lower else (m.V @ k))
        return np.diag(d)

    def dH(self, M):
        """Linear (homogeneous) part of H_N[Gam] applied to a Nambu matrix M (RPA kernel)."""
        m, L = self.m, self.L
        out = np.zeros_like(M)
        if m.U != 0.0:
            out[:L, :L] = -m.U * np.diag(np.diag(M[L:, L:]))
            out[L:, L:] = -m.U * np.diag(np.diag(M[:L, :L]))
        if m.U != 0.0 or m.V is not None:
            out[:L, L:] = self._pairing(M)
            out[L:, :L] = self._pairing(M, lower=True)
        return out

    def A_eff(self, X):
        Q, m = self.Q, self.m
        A = np.zeros((2 * Q, 2 * Q))
        A[np.arange(Q), Q + np.arange(Q)] = m.omega
        A[Q:, :Q] = -np.diag(m.omega) - np.einsum("abc,b->ac", m.kappa, X)
        return A

    # ---------------------------------------------------------------- e-e (Level A kernels)
    def _ee(self, G, x, y, h):
        """Collision terms I_up, I_dn of the first BBGKY equation and the pair-factor
        derivatives, from the Level A pp ladder with one-body field h (incl. phonons)."""
        pp, U, L = self.pp, self.m.U, self.L
        ru, rd = self.rhos(G)
        nu = np.real(np.diag(ru)); nd = np.real(np.diag(rd))
        hu = h + U * np.diag(nd); hd = h + U * np.diag(nu)
        pp._Ucur = U
        if self.ee == "sic":
            dru, drd, dx, dy = pp._rhs_sic(ru, rd, x, y, h)
            Iu = 1j * dru - (hu @ ru - ru @ hu)
            Id = 1j * drd - (hd @ rd - rd @ hd)
        else:
            g = pp._gsel(x, y, ru, rd)
            Iu = U * (g["iiji"] - g["ijjj"])
            Id = U * (g["iiij"] - g["jijj"])
            dx = -1j * (pp._h2(x, hu, hd) + U * pp._pauli_P(x, ru, rd))
            dy = -1j * (pp._h2(y, hu, hd) + U * pp._pauli_P(y, ru, rd))
        return Iu, Id, dx, dy

    def _vertex_ee(self, G, x, y):
        """sum_kl gamma^a_kl <c^+_j c^+_k c_l c_i>_c (opposite-spin k,l) from the pair factors,
        as Nambu matrices (Q, 2L, 2L): the e-e correlation dressing of the e-ph source."""
        L, Q = self.L, self.Q
        ru, rd = self.rhos(G)
        I = np.eye(L)
        out = np.zeros((Q, 2 * L, 2 * L), complex)
        for smp, (Mu, Md), w in ((x, (I - ru, I - rd), self.pp.gw[0]), (y, (ru, rd), self.pp.gw[1])):
            if w == 0 or smp.shape[1] == 0:
                continue
            s = smp.reshape(L, L, -1)
            for a in range(Q):
                g = self.m.gamma[a]
                up = np.einsum("ilr,kl,jkr->ij", s, g, s.conj()) - Mu * np.trace(g @ Md)
                dn = np.einsum("lir,kl,kjr->ij", s, g, s.conj()) - Md * np.trace(g @ Mu)
                out[a, :L, :L] += w * up
                out[a, L:, L:] += -w * dn.T
        return out

    def _dcorr(self, G, x, y, Gph=None):
        L = self.L
        out = np.zeros(L)
        if self.pp is not None:
            ru, rd = self.rhos(G)
            out = out + self.pp._gsel(x if self.ee != "sic" else None, y, ru, rd)["diag"]
        Gph = self.Gph if Gph is None else Gph
        if Gph.size:
            d = np.arange(L) * (L + 1)
            out = out + np.real(Gph[d, d])
        return out

    # ---------------------------------------------------------------- phonon-exchange pair correlation
    def _coll_dense(self, Gm):
        """Collision terms of eq. (1) from a dense up-down pair correlation (as PPLadder g1g2)."""
        L, U = self.L, self.m.U
        G4 = Gm.reshape(L, L, L, L)
        i = np.arange(L)
        Iu = U * (G4[i[:, None], i[:, None], i[None, :], i[:, None]] - G4[i[:, None], i[None, :], i[None, :], i[None, :]])
        Id = U * (G4[i[:, None], i[:, None], i[:, None], i[None, :]] - G4[i[None, :], i[:, None], i[None, :], i[None, :]])
        return Iu, Id

    def _vertex_dense(self, Gm):
        L, Q = self.L, self.Q
        G4 = Gm.reshape(L, L, L, L)
        out = np.zeros((Q, 2 * L, 2 * L), complex)
        for a in range(Q):
            g = self.m.gamma[a]
            out[a, :L, :L] = np.einsum("iljk,kl->ij", G4, g)
            out[a, L:, L:] = -np.einsum("likj,kl->ij", G4, g).T
        return out

    def _rhs_gph(self, t, G, Gph, C, X):
        """eq. (6'): i dG2ph/dt = H2 G2ph - G2ph H2^+ + sum_e ([g^e, rho_u] (x) C_d^e + C_u^e (x) [g^e, rho_d]).
        H2 is the pair Hamiltonian of the e-e ladder (T-matrix; for 'sic' the Hartree and blocking
        fields exclude the pair's own share rho/N, exact for one up and one down electron)."""
        m, L, Q = self.m, self.L, self.Q
        U = m.U
        ru, rd = self.rhos(G)
        cu = cd = 1.0
        if self.ee == "sic":
            cu = 1.0 - 1.0 / self.Nup if self.Nup > 0 else 0.0
            cd = 1.0 - 1.0 / self.Ndn if self.Ndn > 0 else 0.0
        h = self.hph(t, X)
        ou, od = cu * ru, cd * rd
        hu = h + U * np.diag(np.real(np.diag(od))); hd = h + U * np.diag(np.real(np.diag(ou)))
        V = Gph.reshape(L, L, -1)
        HG = (np.einsum("ab,bjr->ajr", hu, V) + np.einsum("jb,abr->ajr", hd, V)).reshape(L * L, -1)
        if U != 0.0:
            d = np.arange(L) * (L + 1)
            I = np.eye(L)
            Blk = np.eye(L * L) - np.kron(ou, I) - np.kron(I, od)
            HG = HG + U * Blk[:, d] @ Gph[d, :]
        src = np.zeros_like(Gph)
        for e in range(Q):
            g = m.gamma[e]
            Cu = C[e, :L, :L]; Cd = -C[e, L:, L:].T
            src += np.kron(g @ ru - ru @ g, Cd) + np.kron(Cu, g @ rd - rd @ g)
        return -1j * (HG - HG.conj().T + src)

    # ---------------------------------------------------------------- sources (shared with the kinetic limit)
    def c_source(self, G, S, x=None, y=None, Gph=None):
        """Inhomogeneous part of eq. (4) for given Gam, S (array (2Q, 2L, 2L))."""
        Q, L = self.Q, self.L
        GN = self.GN
        GNG = GN @ G
        GGN = G @ GN
        comm = GNG - GGN                                                    # [G^b, Gam]
        src = -1j * np.einsum("bm,bij->mij", S[:Q, :], comm)
        # 1/2 [(1-Gam) G^a Gam + Gam G^a (1-Gam)] = 1/2 (G^a Gam + Gam G^a) - Gam G^a Gam
        sym = 0.5 * (GNG + GGN) - G @ GNG
        src[Q:] -= sym
        if self.vertex_ee and self.pp is not None and x is not None:
            src[Q:] -= self._vertex_ee(G, x, y)
        if self.vertex_ee and Gph is not None and Gph.size:
            src[Q:] -= self._vertex_dense(Gph)
        return src

    def t_source(self, S):
        """Inhomogeneous part of eq. (5) (Gaussian + quantum term), symmetrized."""
        Q = self.Q
        F = np.zeros((2 * Q,) * 3)
        Sx = S[:Q, :]
        F[Q:] = -np.einsum("abc,bn,cl->anl", self.m.kappa, Sx, Sx, optimize=True)
        out = F + np.transpose(F, (1, 0, 2)) + np.transpose(F, (2, 1, 0))
        out[Q:, Q:, Q:] += self.m.kappa / 4.0
        return out

    def feedback(self, C, T):
        """Contributions of C, T to dGam/dt (eq. 1) and dS/dt (eq. 3)."""
        Q = self.Q
        dG = -1j * (np.sum(self.GN @ C[:Q], 0) - np.sum(C[:Q] @ self.GN, 0))
        R = np.zeros((2 * Q, 2 * Q))
        R[Q:, :] = -np.real(np.einsum("aij,nji->an", self.GN, C))
        R[Q:, :] -= 0.5 * np.einsum("abc,bcn->an", self.m.kappa, T[:Q, :Q, :])
        return dG, R + R.T

    # ---------------------------------------------------------------- equations of motion
    def state(self):
        return [self.G, self.x, self.y, self.Gph, self.Phi, self.S, self.T, self.C, self.W]

    def set_state(self, s):
        self.G, self.x, self.y, self.Gph, self.Phi, self.S, self.T, self.C, self.W = s

    def rhs(self, t, st):
        G, x, y, Gph, Phi, S, T, C, W = st
        m, L, Q = self.m, self.L, self.Q
        X = Phi[:Q]
        H = self.hfb(t, G, X)
        dG = -1j * (H @ G - G @ H)
        dx = 0.0 * x; dy = 0.0 * y
        if self.pp is not None:
            Iu, Id, dx, dy = self._ee(G, x, y, self.hph(t, X))
            dG[:L, :L] += -1j * Iu
            dG[L:, L:] += 1j * Id.T
        dGph = 0.0 * Gph
        if self.g2ph:
            if m.U != 0.0:
                Iu, Id = self._coll_dense(Gph)
                dG[:L, :L] += -1j * Iu
                dG[L:, L:] += 1j * Id.T
            dGph = self._rhs_gph(t, G, Gph, C, X)
        dPhi = np.zeros(2 * Q); dS = np.zeros_like(S); dT = np.zeros_like(T); dC = np.zeros_like(C)
        if Q:
            A = self.A_eff(X)
            Fe = np.real(np.einsum("aij,ji->a", self.GN, G)) + self.trg
            dPhi[:Q] = m.omega * Phi[Q:]
            dPhi[Q:] = -m.omega * X - m.f(t) - Fe - 0.5 * np.einsum("abc,bc->a", m.kappa,
                                                                     np.outer(X, X) + S[:Q, :Q])
            dS = A @ S + S @ A.T
            if self.eph_corr or self.phph_corr:
                fG, fS = self.feedback(C, T)
                dG += fG
                dS += fS
            if self.eph_corr:
                dC = -1j * (H @ C - C @ H)
                dHC = np.array([self.dH(c) for c in C]) if (m.U != 0.0 or m.V is not None) else None
                if dHC is not None:
                    dC += -1j * (np.einsum("mij,jk->mik", dHC, G) - np.einsum("ij,mjk->mik", G, dHC))
                dC += np.einsum("mn,nij->mij", A, C)
                dC += self.c_source(G, S, x, y, Gph)
            if self.phph_corr:
                dT = (np.einsum("mk,knl->mnl", A, T) + np.einsum("nk,mkl->mnl", A, T)
                      + np.einsum("lk,mnk->mnl", A, T)) + self.t_source(S)
        # work of the external fields: <dH/dt>
        eps = 1e-5
        hd = (m.h(t + eps) - m.h(t - eps)) / (2 * eps)
        ru, rd = self.rhos(G)
        dW = float(np.real(np.trace(hd @ (ru + rd))))
        if Q and m.force is not None:
            dW += float((m.f(t + eps) - m.f(t - eps)) @ Phi[:Q] / (2 * eps))
        return [dG, dx, dy, dGph, dPhi, dS, dT, dC, dW]

    def step(self, dt):
        t = self.t
        if self.E0 is None:
            self.E0 = self.energy()
        s0 = self.state()
        k1 = self.rhs(t, s0)
        k2 = self.rhs(t + dt / 2, [a + dt / 2 * b for a, b in zip(s0, k1)])
        k3 = self.rhs(t + dt / 2, [a + dt / 2 * b for a, b in zip(s0, k2)])
        k4 = self.rhs(t + dt, [a + dt * b for a, b in zip(s0, k3)])
        self.set_state([a + dt / 6 * (b + 2 * c + 2 * d + e) for a, b, c, d, e in zip(s0, k1, k2, k3, k4)])
        self.t = t + dt
        if self.pp is not None and self.sum_rules is not None:
            self._apply_sum_rules()
        if self.positivity and self.Q and self.eph_corr:
            self._project_positive()

    def _coupling_cov(self):
        """<dO_a dO_b> for the e-ph coupling operators (Wick part + e-e pair correlation)."""
        Q, L = self.Q, self.L
        G = self.G
        Gb = np.eye(2 * L) - G
        M = np.real(np.einsum("aij,jk,bkl,li->ab", self.GN, Gb, self.GN, G, optimize=True))
        vert = None
        if self.pp is not None:
            vert = self._vertex_ee(G, self.x, self.y)
        if self.Gph.size:
            v2 = self._vertex_dense(self.Gph)
            vert = v2 if vert is None else vert + v2
        if vert is not None:
            # opposite-spin pair part: 2 Re sum_ij gamma^b_ji vertex^a_ij (up block), as <dO dO> is symmetric in spin
            M = M + 2 * np.real(np.einsum("bji,aij->ab", self.m.gamma, vert[:, :L, :L]))
        return 0.5 * (M + M.T)

    def _project_positive(self):
        Q = self.Q
        Om = np.block([[np.zeros((Q, Q)), np.eye(Q)], [-np.eye(Q), np.zeros((Q, Q))]])
        Hs = self.S + 0.5j * Om
        w, V = np.linalg.eigh(Hs)
        if w.min() < 0:
            Hs = (V * np.maximum(w, 0)) @ V.conj().T
            self.S = np.real(0.5 * (Hs + Hs.conj().T))
            self.S = 0.5 * (self.S + self.S.T)
            Hs = self.S + 0.5j * Om
        Moo = self._coupling_cov()
        Mox = np.real(np.einsum("aij,nji->an", self.GN, self.C))

        def mineig(lam):
            M = np.block([[Moo, lam * Mox], [lam * Mox.T, Hs]])
            return np.linalg.eigvalsh(M).min()
        if mineig(1.0) >= -1e-12 or mineig(0.0) < -1e-12:
            return
        lo, hi = 0.0, 1.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if mineig(mid) >= -1e-12:
                lo = mid
            else:
                hi = mid
        self.C = self.C * lo

    def _apply_sum_rules(self):
        """Level A constraints (pair sum rules; energy) on the sic pair amplitudes.  The energy
        target is the total energy E0 + W minus all non-Hubbard terms, so that with phonons
        and fields switched off it is exactly Level A's E0."""
        pp, L = self.pp, self.L
        ru, rd = self.rhos()
        pp.ru, pp.rd, pp.y = ru, rd, self.y
        pp._tcur = self.t
        if self.sum_rules == "energy":
            parts = self.energy(parts=True)
            pp.E0 = self.E0 + self.W - (parts["total"] - parts["kin"] - parts["U"])
            pp._Ucur = self.m.U
        else:
            pp.E0 = 0.0
            pp._Ucur = 0.0
        pp._sum_rules()
        pp._Ucur = self.m.U
        self.y = pp.y

    # ---------------------------------------------------------------- observables
    def densities(self):
        ru, rd = self.rhos()
        return np.real(np.diag(ru)), np.real(np.diag(rd))

    def double_occ(self):
        nu, nd = self.densities()
        L = self.L
        kap = np.diag(self.G[:L, L:])
        return nu * nd + np.abs(kap) ** 2 + self._dcorr(self.G, self.x, self.y)

    def phonon_number(self):
        """<b^+_a b_a> = (X^2 + P^2 + S_xx + S_pp - 1)/2 per mode."""
        Q = self.Q
        X, P = self.Phi[:Q], self.Phi[Q:]
        return 0.5 * (X ** 2 + P ** 2 + np.diag(self.S)[:Q] + np.diag(self.S)[Q:] - 1.0)

    def gap(self):
        """Pairing field Delta_i (diagonal of the anomalous block of H_N without external field)."""
        return np.diag(self._pairing(self.G))

    def energy(self, parts=False):
        m, L, Q = self.m, self.L, self.Q
        t = self.t
        G, Phi, S, T, C = self.G, self.Phi, self.S, self.T, self.C
        ru, rd = self.rhos()
        kap = np.diag(G[:L, L:])
        e = {}
        e["kin"] = float(np.real(np.trace(m.h(t) @ (ru + rd))))
        nu, nd = np.real(np.diag(ru)), np.real(np.diag(rd))
        e["U"] = m.U * float(np.sum(nu * nd + np.abs(kap) ** 2 + self._dcorr(G, self.x, self.y)))
        e["V"] = 0.0 if m.V is None else float(np.real(np.conj(kap) @ m.V @ kap))
        if m.pair_field is not None:
            Fp = m.pair_field(t)
            e["V"] += 2 * float(np.real(np.sum(Fp * np.conj(G[:L, L:]))))
        if Q:
            X, P = Phi[:Q], Phi[Q:]
            e["ph"] = float(np.sum(m.omega * (X ** 2 + P ** 2 + np.diag(S)[:Q] + np.diag(S)[Q:] - 1.0)) / 2)
            Fe = np.real(np.einsum("aij,ji->a", self.GN, G)) + self.trg
            e["eph"] = float(X @ Fe + np.sum(np.real(np.einsum("aij,aji->a", self.GN, C[:Q]))))
            k = m.kappa
            e["anh"] = float((np.einsum("abc,a,b,c->", k, X, X, X) + 3 * np.einsum("abc,a,bc->", k, X, S[:Q, :Q])
                              + np.einsum("abc,abc->", k, T[:Q, :Q, :Q])) / 6)
            e["force"] = float(m.f(t) @ X)
        e["total"] = sum(e.values())
        return e if parts else e["total"]

    # ---------------------------------------------------------------- Level C: equilibrium
    def equilibrium(self, T, mu=0.0, iters=2000, tol=1e-12, mix=0.5):
        """Self-consistent finite-temperature HFB state (stationary point of eq. 1 without
        correlations): Gam = f(H_N[Gam]), f the Fermi function at T (Nambu, chemical potential
        mu folded into h).  Returns the pairing field."""
        L = self.L
        X = self.Phi[:self.Q]
        G = self.G
        shift = np.diag(np.concatenate([-mu * np.ones(L), mu * np.ones(L)]))
        for it in range(iters):
            H = self.hfb(self.t, G, X) + shift
            E, V = np.linalg.eigh(H)
            f = 0.5 * (1 - np.tanh(E / (2 * T))) if T > 0 else (E < 0).astype(float)
            Gn = (V * f) @ V.conj().T
            err = np.abs(Gn - G).max()
            G = mix * Gn + (1 - mix) * G
            if err < tol:
                break
        self.G = G
        return self.gap()

    # ---------------------------------------------------------------- Level B: kinetic limit
    def kinetic_rhs(self, f, N, eta, kernel="gauss", t=None):
        """Markov + random-phase limit of eqs. (1)-(5).

        f: occupations of the eigenstates of H_N (2L Nambu), N: phonon occupations (Q).  The
        state is built with Gam diagonal in the mean-field eigenbasis and S = (N + 1/2) in the
        normal modes (random phases), the sources of (4) and (5) are evaluated with the same
        functions as in the full dynamics, the linear equations are solved adiabatically with
        memory kernel m(tau) = exp(-eta tau) ('lorentz') or exp(-eta^2 tau^2 / 2) ('gauss'),
        and C, T are inserted into (1) and (3).  Returns df/dt, dN/dt and the eigenbasis.
        Valid for harmonic phonon normal modes (X = 0)."""
        t = self.t if t is None else t
        L, Q = self.L, self.Q
        H = self.hfb(t, self.G, np.zeros(Q))
        eps, Ue = np.linalg.eigh(H)
        G = (Ue * f) @ Ue.conj().T
        S = np.diag(np.concatenate([N + 0.5, N + 0.5]))
        A = self.A_eff(np.zeros(Q))
        lam, Vp = np.linalg.eig(A)
        Vi = np.linalg.inv(Vp)

        def K(z):          # int_0^inf exp(z tau) m(tau) dtau, z purely imaginary: z = -i Delta
            D = np.imag(-z)
            if kernel == "lorentz":
                return 1.0 / (eta + 1j * D)
            return np.sqrt(np.pi / 2) / eta * np.exp(-D ** 2 / (2 * eta ** 2)) - 1j * np.sqrt(2) / eta * dawsn(D / (np.sqrt(2) * eta))
        C = np.zeros((2 * Q, 2 * L, 2 * L), complex)
        T = np.zeros((2 * Q,) * 3)
        if self.eph_corr:
            src = self.c_source(G, S)
            st = np.einsum("mn,nij->mij", Vi, src)
            st = Ue.conj().T @ st @ Ue
            z = -1j * (eps[None, :, None] - eps[None, None, :]) + lam[:, None, None]
            Ct = K(z) * st
            Ct = Ue @ Ct @ Ue.conj().T
            C = np.einsum("mn,nij->mij", Vp, Ct)
        if self.phph_corr:
            src = self.t_source(S)
            st = np.einsum("am,bn,cl,mnl->abc", Vi, Vi, Vi, src, optimize=True)
            z = lam[:, None, None] + lam[None, :, None] + lam[None, None, :]
            Tt = K(z) * st
            T = np.real(np.einsum("am,bn,cl,mnl->abc", Vp, Vp, Vp, Tt, optimize=True))
        dG, dS = self.feedback(C, T)
        df = np.real(np.einsum("ia,ij,ja->a", Ue.conj(), dG, Ue))
        dN = 0.5 * (np.diag(dS)[:Q] + np.diag(dS)[Q:])
        return df, dN, eps, Ue


# ======================================================================================
# Material scale: Bloch electrons + correlations in the Markov / random-phase limit
# ======================================================================================
# The same eqs. (1)-(6), for a crystal with many k-points and phonon modes, in the limit that
# is tractable over a full Brillouin zone:
#  * electrons: coherent Bloch density matrices rho_k (eq. 1 with the light field), with the
#    collision terms of eqs. (4) and (6) in the Markov limit, evaluated with energy-binned
#    occupations F(e) (isotropic random-phase limit in k at fixed energy) and applied to each k
#    as in/out scattering of populations and the corresponding dephasing of coherences;
#  * e-e: Markov limit of the second-order G2 equation with a local (on-site) interaction U,
#        Gamma_out(e) = 2 pi U^2 int de2 de3 rho(e2) rho(e3) rho(e4) F2 (1-F3)(1-F4),  e4 = e + e2 - e3,
#        Gamma_in(e)  = 2 pi U^2 int ...                    (1-F2) F3 F4         (rho per site, per spin);
#  * e-ph: Markov limit of eq. (4) with an energy-resolved coupling kernel W_q(e) (transitions
#        e -> e + w_q), i.e. dN_q/dt = sum_e W_q(e) [F(e+w)(1-F(e))(N+1) - F(e)(1-F(e+w)) N];
#  * ph-ph: Markov limit of eq. (5), three-phonon decay/coalescence on process lists
#        (a, b, c, w_abc) with w = (pi/8)|kappa_abc|^2 delta(...) (kinetic_rhs normalization).
# All three conserve particle number and energy (up to the energy-bin discretization of the
# phonon energies, handled by splitting each transition between the two neighbouring bins).

def markov_ee_local(F, rho, U, de):
    """In/out scattering rates per state (energy units, hbar = 1) on a uniform energy grid.
    F: occupations per bin, rho: DOS per site and spin (1/energy), U: on-site interaction."""
    a = rho * F
    b = rho * (1 - F)
    n = len(F)
    bb = np.convolve(b, b)
    aa = np.convolve(a, a)                  # index m = i3 + i4
    idx = np.arange(n)[:, None] + np.arange(n)[None, :]
    gout = 2 * np.pi * U ** 2 * de ** 2 * (bb[idx] @ a)
    gin = 2 * np.pi * U ** 2 * de ** 2 * (aa[idx] @ b)
    return gin, gout


class EphKernel:
    """Energy-resolved e-ph Markov collision operator.  W: (Q, nbins) transition kernel per
    mode for e -> e + w (absorption direction) on the electron grid (spacing de); w: (Q,) phonon
    energies; g: states per bin per cell (both spins); nmodes: modes per branch (Nq), which
    normalizes the electron rates per cell.  Fractional energy shifts are split between the two
    neighbouring bins; each part adds (s de / w) phonons of energy w, so that the electron
    number and the total energy are conserved exactly."""

    def __init__(self, W, w, de, g, nmodes):
        self.W, self.w, self.de, self.g, self.Nq = W, w, de, g, nmodes
        s = w / de
        self.s0 = np.floor(s).astype(int)
        self.r = s - self.s0
        self.groups = []
        act = W.sum(1) > 0
        nb = W.shape[1]
        for k, wt in ((0, 1 - self.r), (1, self.r)):
            sh = self.s0 + k
            for s_ in np.unique(sh[act]):
                qs = np.nonzero((sh == s_) & act & (wt > 1e-12))[0]
                if len(qs) and 0 < s_ < nb:
                    self.groups.append((s_, qs, W[qs, :nb - s_] * wt[qs, None], s_ * de / w[qs]))

    def rates(self, F, N):
        """dN/dt (Q,) and electron in/out rates per state (nbins,) (time units of W)."""
        nb = len(F)
        dN = np.zeros(len(N)); rin = np.zeros(nb); rout = np.zeros(nb)
        for s, qs, Wq, fac in self.groups:
            Flo, Fhi = F[:nb - s], F[s:]
            Nq_ = N[qs][:, None]
            em = Wq * (Fhi * (1 - Flo))[None, :] * (Nq_ + 1)        # hi -> lo + phonon
            ab = Wq * (Flo * (1 - Fhi))[None, :] * Nq_              # lo + phonon -> hi
            dN[qs] += fac * (em - ab).sum(1)
            emt = em.sum(0) / self.Nq
            abt = ab.sum(0) / self.Nq
            rout[s:] += emt; rin[:nb - s] += emt
            rout[:nb - s] += abt; rin[s:] += abt
        g = self.g
        Gout = np.where(F * g > 1e-300, rout / np.maximum(F * g, 1e-300), 0.0)
        Gin = np.where((1 - F) * g > 1e-300, rin / np.maximum((1 - F) * g, 1e-300), 0.0)
        return dN, Gin, Gout


def markov_phph(N, procs):
    """Three-phonon Markov collision integral on process lists (kinds 0: a -> b + c,
    1: a + b -> c, 2: a + c -> b) with weights in the kinetic_rhs normalization."""
    out = np.zeros(len(N))
    for kind, (a, b, c, v) in enumerate(procs):
        na, n1, n2 = N[a], N[b], N[c]
        if kind == 0:
            r = (na + 1) * n1 * n2 - na * (n1 + 1) * (n2 + 1)
        elif kind == 1:
            r = (na + 1) * (n1 + 1) * n2 - na * n1 * (n2 + 1)
        else:
            r = (na + 1) * n1 * (n2 + 1) - na * (n1 + 1) * n2
        out += np.bincount(a, weights=v * r, minlength=len(N))
    return out


class BlochCarriers:
    """Coherent two-band Bloch electrons rho_k (2x2, orbital basis, spin degenerate) driven by a
    vector potential through h(k + a(t)) (Peierls), with Markov collision terms applied in the
    instantaneous band basis: populations relax as dF_k = -Gout F_k + Gin (1 - F_k), coherences
    decay at [(Gin + Gout)_c + (Gin + Gout)_v] / 2 (random-phase Markov limit of eqs. 4 and 6).
    hk_bands: (kc + a) -> (energies (...,2) ascending, eigenvectors (...,2,2) as columns)."""

    def __init__(self, kc, hk_bands, T0, egrid):
        self.kc = kc
        self.bands = hk_bands
        e, v = hk_bands(kc)
        f = 1 / (np.exp(np.clip(e / T0, -700, 700)) + 1)
        self.rho = np.einsum("kam,km,kbm->kab", v, f, v.conj())
        self.eg = egrid
        self.de = egrid[1] - egrid[0]

    def step_coherent(self, a, dt):
        """exact propagation with the instantaneous h(k + a) over dt (hbar = 1)."""
        e, v = self.bands(self.kc + a)
        ph = np.exp(-1j * e * dt)
        Uk = np.einsum("kam,km,kbm->kab", v, ph, v.conj())
        self.rho = Uk @ self.rho @ np.conj(np.transpose(Uk, (0, 2, 1)))

    def band_occupations(self, a):
        e, v = self.bands(self.kc + a)
        rb = np.einsum("kam,kab,kbn->kmn", v.conj(), self.rho, v)
        return e, v, rb

    def _bins(self, e):
        return np.clip(np.rint((e - self.eg[0]) / self.de).astype(int), 0, len(self.eg) - 1)

    def binned(self, a):
        """occupation summed per energy bin and number of states per bin (random-phase projection)."""
        e, v, rb = self.band_occupations(a)
        ib = self._bins(e)
        occ = np.real(np.stack([rb[:, 0, 0], rb[:, 1, 1]], -1))
        cnt = np.bincount(ib.ravel(), minlength=len(self.eg)).astype(float)
        s = np.bincount(ib.ravel(), weights=occ.ravel(), minlength=len(self.eg))
        return s, cnt

    def step_collisions(self, a, Gin, Gout, dt):
        """exact exponential relaxation over dt with bin rates Gin, Gout (per state)."""
        e, v, rb = self.band_occupations(a)
        ib = self._bins(e)
        gi, go = Gin[ib], Gout[ib]
        tot = gi + go
        feq = np.where(tot > 0, gi / np.where(tot > 0, tot, 1), 0.0)
        dec = np.exp(-tot * dt)
        for m in range(2):
            rb[:, m, m] = feq[:, m] + (np.real(rb[:, m, m]) - feq[:, m]) * dec[:, m]
        deph = np.exp(-0.5 * (tot[:, 0] + tot[:, 1]) * dt)
        rb[:, 0, 1] *= deph
        rb[:, 1, 0] *= deph
        self.rho = np.einsum("kam,kmn,kbn->kab", v, rb, v.conj())

    def energy(self, a):
        """band energy per cell (both spins)."""
        e, v, rb = self.band_occupations(a)
        return 2 * float(np.mean(np.real(rb[:, 0, 0]) * e[:, 0] + np.real(rb[:, 1, 1]) * e[:, 1]))

    def coherence(self, a):
        _, _, rb = self.band_occupations(a)
        return float(np.mean(np.abs(rb[:, 0, 1])))
