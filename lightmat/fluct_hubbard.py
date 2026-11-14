"""Positive-factorized time-linear NEGF for Hubbard models (pp ladder channel).

Notation (hbar = 1). rho^s_ij = <c^+_js c_is>, rhobar = 1 - rho.
Two-particle pp objects on pairs (i up, j down), stored as L x L matrices / vectors of
length L^2:
    A = rhobar^up (x) rhobar^dn     (uncorrelated hole-hole pair matrix, chi^>_0)
    B = rho^up    (x) rho^dn        (uncorrelated particle-particle pair matrix, chi^<_0)
    G = connected part of <c^+_k up c^+_l dn c_j dn c_i up>   (pair correlation)
The G1-G2 T-matrix (TPP) equation
    i dG/dt = [h2, G] + U (A-B) P G - U G P (A-B) + U (A P B - B P A)
is identical to the pair of *linear* equations
    i dX/dt = (h2 - U Y P) X - X (h2 - U Y P)^+ ,   X = A + G
    i dY/dt = (h2 + U X P) Y - Y (h2 + U X P)^+ ,   Y = B + G
(P = projector on on-site pairs, h2 = h^HF_up (x) 1 + 1 (x) h^HF_dn).
Writing X = x x^+, Y = y y^+ with rectangular sample matrices x, y that evolve as
    i dx/dt = (h2 - U Y P) x,     i dy/dt = (h2 + U X P) y
keeps X, Y positive semidefinite for all times (the instability of plain G1-G2 is a
loss of positivity). For uncorrelated initial states the factorization is exact and of
low rank (rank X = #hole pairs, rank Y = #particle pairs), so no sampling noise.
The correlation used for observables is the Keldysh average G = (X - A + Y - B)/2, and
the one-particle density matrices evolve with the exact collision integral
    i d rho^up_ij/dt = [h^HF_up, rho^up]_ij + U (G_{ii,ji} - G_{ij,jj}),
    i d rho^dn_ij/dt = [h^HF_dn, rho^dn]_ij + U (G_{ii,ij} - G_{ji,jj}).
G_{ab,cd} indices: a = up annihilated, b = down annihilated, c = up created,
d = down created.  Double occupation d_i = n_i^up n_i^dn + G_{ii,ii}.
"""
from __future__ import annotations
import numpy as np


class PPLadder:
    def __init__(self, h, U, rho_up, rho_dn, rank_tol=1e-12, mode="factorized",
                 contract=False, selfconsistent=False, lyap_eps=1e-2,
                 gweights=None):
        self.selfconsistent = selfconsistent
        self.gw = (0.5, 0.5) if gweights is None else gweights
        self.lyap_eps = lyap_eps
        """h: (L,L) one-body hopping matrix (may be replaced by a callable h(t)).
        mode: 'factorized' (positive x,y samples) or 'g1g2' (plain G1-G2 TPP for G)."""
        self.hfun = h if callable(h) else (lambda t, h=h: h)
        self.U = U
        self.L = L = rho_up.shape[0]
        self.mode = mode
        self.contract = contract
        self.Nu = float(np.real(np.trace(rho_up)))
        self.Nd = float(np.real(np.trace(rho_dn)))
        self.ru = rho_up.astype(complex)
        self.rd = rho_dn.astype(complex)
        if mode == "factorized":
            I = np.eye(L)
            fu, fd = self._factor(I - self.ru, rank_tol), self._factor(I - self.rd, rank_tol)
            self.x = np.einsum("ar,bs->abrs", fu, fd).reshape(L * L, -1)
            fu, fd = self._factor(self.ru, rank_tol), self._factor(self.rd, rank_tol)
            self.y = np.einsum("ar,bs->abrs", fu, fd).reshape(L * L, -1)
        else:
            self.G = np.zeros((L * L, L * L), complex)
        self.diag = np.arange(L) * (L + 1)  # flat indices of on-site pairs (i,i)

    # -- helpers ---------------------------------------------------------
    @staticmethod
    def _factor(M, tol):
        w, v = np.linalg.eigh(M)
        keep = w > tol
        return (v[:, keep] * np.sqrt(w[keep])).astype(complex)

    def _AB(self, ru=None, rd=None):
        ru = self.ru if ru is None else ru
        rd = self.rd if rd is None else rd
        L = self.L
        I = np.eye(L)
        return np.kron(I - ru, I - rd), np.kron(ru, rd)

    def _h2(self, v, hu, hd):
        """apply h2 to columns of v (L^2 x r): reshape each column to L x L matrix V_{ij}
        (i up, j down): h2 V = hu V + V hd^T."""
        L = self.L
        V = v.reshape(L, L, -1)
        out = np.einsum("ab,bjr->ajr", hu, V) + np.einsum("jb,abr->ajr", hd, V)
        return out.reshape(L * L, -1)

    def contract_rho(self, x, y):
        """One-particle density matrices from the pair matrices via the exact
        contraction sum rules  Tr_dn Y = N_dn rho_up,  Tr_dn X = (L - N_dn) rhobar_up
        (and the same with up <-> dn).  Positivity of X, Y then implies 0 <= rho <= 1
        whenever the two estimates agree; we use their average."""
        L = self.L
        X = (x @ x.conj().T).reshape(L, L, L, L)
        Y = (y @ y.conj().T).reshape(L, L, L, L)
        I = np.eye(L)
        def est(Nother, trY, trX):
            r1 = trY / Nother if Nother > 0 else None
            r2 = I - trX / (L - Nother) if Nother < L else None
            if r1 is None:
                return r2
            if r2 is None:
                return r1
            return 0.5 * (r1 + r2)
        ru = est(self.Nd, np.einsum("ajcj->ac", Y), np.einsum("ajcj->ac", X))
        rd = est(self.Nu, np.einsum("jajc->ac", Y), np.einsum("jajc->ac", X))
        return ru, rd

    def correlation(self, x=None, y=None, ru=None, rd=None):
        if self.mode != "factorized":
            return self.G
        x = self.x if x is None else x
        y = self.y if y is None else y
        A, B = self._AB(ru, rd)
        return 0.5 * (x @ x.conj().T - A + y @ y.conj().T - B)

    # -- equations of motion --------------------------------------------
    def rhs(self, t, state):
        L, U = self.L, self.U
        h = self.hfun(t)
        if self.mode == "factorized":
            ru, rd, x, y = state
            if self.contract:
                ru, rd = self.contract_rho(x, y)
        else:
            ru, rd, G = state
        nu = np.real(np.diag(ru)); nd = np.real(np.diag(rd))
        hu = h + U * np.diag(nd)
        hd = h + U * np.diag(nu)
        if self.mode == "factorized":
            return self._rhs_fact(ru, rd, x, y, hu, hd)
        if self.mode == "factorized_old":
            A, B = self._AB(ru, rd)
            Gc = 0.5 * (x @ x.conj().T - A + y @ y.conj().T - B)
        else:
            Gc = G
        G4 = Gc.reshape(L, L, L, L)  # [a,b,c,d]
        idx = np.arange(L)
        # collision integrals
        # I^up_ij = U (G_{ii,ji} - G_{ij,jj})
        Giiji = G4[idx[:, None], idx[:, None], idx[None, :], idx[:, None]]  # [i,j] -> G[i,i,j,i]
        Gijjj = G4[idx[:, None], idx[None, :], idx[None, :], idx[None, :]]  # G[i,j,j,j]
        Iu = U * (Giiji - Gijjj)
        # I^dn_ij = U (G_{ii,ij} - G_{ji,jj})
        Giiij = G4[idx[:, None], idx[:, None], idx[:, None], idx[None, :]]  # G[i,i,i,j]
        Gjijj = G4[idx[None, :], idx[:, None], idx[None, :], idx[None, :]]  # G[j,i,j,j]
        Id = U * (Giiij - Gjijj)
        dru = -1j * (hu @ ru - ru @ hu + Iu)
        drd = -1j * (hd @ rd - rd @ hd + Id)
        if self.mode == "factorized":
            X = x @ x.conj().T
            Y = y @ y.conj().T
            # (Y P x): Y[:, diag] @ x[diag, :]
            dx = -1j * (self._h2(x, hu, hd) - U * Y[:, self.diag] @ x[self.diag, :])
            dy = -1j * (self._h2(y, hu, hd) + U * X[:, self.diag] @ y[self.diag, :])
            return [dru, drd, dx, dy]
        A, B = self._AB(ru, rd)
        D = A - B
        h2G = self._h2(Gc, hu, hd)
        Gh2 = self._h2(Gc.conj().T, hu, hd).conj().T
        PG = np.zeros_like(Gc); PG[self.diag] = Gc[self.diag]
        DPG = D[:, self.diag] @ Gc[self.diag, :]
        GPD = Gc[:, self.diag] @ D[self.diag, :]
        Psi = A[:, self.diag] @ B[self.diag, :] - B[:, self.diag] @ A[self.diag, :]
        dG = -1j * (h2G - Gh2 + U * (DPG - GPD) + U * Psi)
        return [dru, drd, dG]

    # -- efficient factorized right-hand side (never forms L^2 x L^2 objects) ----
    def _gsel(self, x, y, ru, rd):
        """Elements of G = (X-A+Y-B)/2 needed for collision integrals and d_corr."""
        L = self.L
        idx = np.arange(L)
        ubar = np.eye(L) - ru; dbar = np.eye(L) - rd
        out = {}
        for key in ("iiji", "ijjj", "iiij", "jijj", "diag"):
            out[key] = 0.0
        for s, (Mu, Md), wgt in ((x, (ubar, dbar), self.gw[0]), (y, (ru, rd), self.gw[1])):
            s4 = s.reshape(L, L, -1)
            sd = s4[idx, idx]          # (L, r): pair (p,p)
            cs4 = s4.conj(); csd = sd.conj()
            du = np.diag(Mu); dd = np.diag(Md)
            out["iiji"] = out["iiji"] + wgt * (np.einsum("ir,jir->ij", sd, cs4) - Mu * dd[:, None])
            out["ijjj"] = out["ijjj"] + wgt * (np.einsum("ijr,jr->ij", s4, csd) - Mu * dd[None, :])
            out["iiij"] = out["iiij"] + wgt * (np.einsum("ir,ijr->ij", sd, cs4) - du[:, None] * Md)
            out["jijj"] = out["jijj"] + wgt * (np.einsum("jir,jr->ij", s4, csd) - du[None, :] * Md)
            out["diag"] = out["diag"] + wgt * (np.real(np.einsum("ir,ir->i", sd, csd)) - np.real(du * dd))
        return out

    def _pauli_P(self, s, ru, rd):
        """((A - B) P s) for samples s (L^2 x r)."""
        L = self.L
        idx = np.arange(L)
        ubar = np.eye(L) - ru; dbar = np.eye(L) - rd
        sd = s.reshape(L, L, -1)[idx, idx]
        out = np.einsum("ap,bp,pr->abr", ubar, dbar, sd) - np.einsum("ap,bp,pr->abr", ru, rd, sd)
        return out.reshape(L * L, -1)

    def _rhs_fact(self, ru, rd, x, y, hu, hd):
        U = self.U
        g = self._gsel(x, y, ru, rd)
        Iu = U * (g["iiji"] - g["ijjj"])
        Id = U * (g["iiij"] - g["jijj"])
        dru = -1j * (hu @ ru - ru @ hu + Iu)
        drd = -1j * (hd @ rd - rd @ hd + Id)
        dx = -1j * (self._h2(x, hu, hd) + U * self._pauli_P(x, ru, rd))
        dy = -1j * (self._h2(y, hu, hd) + U * self._pauli_P(y, ru, rd))
        if self.selfconsistent:
            # Keep X - A and Y - B equal to the same correlation G when rho changes by
            # collisions: add generators K with K A + A K^+ = dA/dt|_coll (Lyapunov
            # equations on the L x L one-particle factors), applied to the samples.
            cu = -1j * Iu; cd = -1j * Id            # collisional d rho / dt (Hermitian)
            L = self.L
            I = np.eye(L)
            Ku_h = self._lyap(I - ru, -cu); Kd_h = self._lyap(I - rd, -cd)
            Ku_p = self._lyap(ru, cu); Kd_p = self._lyap(rd, cd)
            dx = dx + self._h2(x, Ku_h, Kd_h)
            dy = dy + self._h2(y, Ku_p, Kd_p)
        return [dru, drd, dx, dy]

    def _lyap(self, M, C):
        """Hermitian K with K M + M K = C (M Hermitian PSD), regularized."""
        w, V = np.linalg.eigh(M)
        Cp = V.conj().T @ C @ V
        den = w[:, None] + w[None, :]
        Kp = Cp / np.maximum(den, self.lyap_eps)
        return V @ Kp @ V.conj().T

    def state(self):
        if self.mode == "factorized":
            return [self.ru, self.rd, self.x, self.y]
        return [self.ru, self.rd, self.G]

    def set_state(self, s):
        if self.mode == "factorized":
            self.ru, self.rd, self.x, self.y = s
        else:
            self.ru, self.rd, self.G = s

    def step(self, t, dt):
        s0 = self.state()
        k1 = self.rhs(t, s0)
        k2 = self.rhs(t + dt / 2, [a + dt / 2 * b for a, b in zip(s0, k1)])
        k3 = self.rhs(t + dt / 2, [a + dt / 2 * b for a, b in zip(s0, k2)])
        k4 = self.rhs(t + dt, [a + dt * b for a, b in zip(s0, k3)])
        self.set_state([a + dt / 6 * (b + 2 * c + 2 * d + e)
                        for a, b, c, d, e in zip(s0, k1, k2, k3, k4)])

    # -- observables -----------------------------------------------------
    def rhos(self):
        if self.mode == "factorized" and self.contract:
            return self.contract_rho(self.x, self.y)
        return self.ru, self.rd

    def densities(self):
        ru, rd = self.rhos()
        return np.real(np.diag(ru)), np.real(np.diag(rd))

    def dcorr(self):
        ru, rd = self.rhos()
        if self.mode == "factorized":
            return self._gsel(self.x, self.y, ru, rd)["diag"]
        G = self.correlation(ru=ru, rd=rd)
        return np.real(G[self.diag, self.diag])
