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
                 gweights=None, Ufun=None):
        self.Ufun = Ufun
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
        self.sic = (mode == "sic")
        if self.sic:
            mode = self.mode = "factorized"
            self.gw = (0.0, 1.0)
            self.sum_rules = True          # exact pair sum rules + energy conservation
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
        L = self.L
        U = self.U if self.Ufun is None else self.Ufun(t)
        self._Ucur = U
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
            if getattr(self, "sic", False):
                return self._rhs_sic(ru, rd, x, y, h)
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
            if s is None or wgt == 0:
                continue
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

    def _rhs_sic(self, ru, rd, x, y, h):
        """Self-interaction-free pair ladder (mode 'sic').

        Each column psi_r of y is a pair amplitude (up index a, down index b).  The third
        particles seen by a pair are all electrons *except the pair's own two*.  With the
        pair's own one-body matrices o^up_r = Tr_dn psi psi^+ / |psi|^2 (o^dn_r alike) and
        r_r = rho - o_r,
          i dpsi/dt = [h (x) 1 + 1 (x) h + U diag(n_dn - o^dn_r) (x) 1 + 1 (x) U diag(n_up - o^up_r)] psi
                      + U (1 - r^up_r (x) 1 - 1 (x) r^dn_r) P psi.
        For one up and one down electron the 'other' densities vanish and this is the exact
        two-body Schroedinger equation; the T-matrix closure is not exact there, because each
        electron feels the Hartree field of its own partner (which makes breaking a doublon
        resonant instead of costing U).  The correlation is G = Y - B, and the one-body
        matrices follow the exact first BBGKY equation.  After every step the exact sum rules
        and energy conservation are restored (see _sum_rules)."""
        U = self._Ucur
        L = self.L
        r = y.shape[1]
        Y4 = y.reshape(L, L, r)
        nrm = np.maximum(np.real(np.einsum("abr,abr->r", Y4, Y4.conj())), 1e-300)
        ou = np.einsum("abr,cbr->acr", Y4, Y4.conj()) / nrm
        od = np.einsum("abr,adr->bdr", Y4, Y4.conj()) / nrm
        nu = np.real(np.diag(ru)); nd = np.real(np.diag(rd))
        vu = nd[:, None] - np.real(np.einsum("aar->ar", od))      # others' dn density seen by up
        vd = nu[:, None] - np.real(np.einsum("bbr->br", ou))      # others' up density seen by dn
        idx = np.arange(L)
        s = Y4[idx, idx]                                           # (L, r) on-site amplitudes
        ruo = ru[:, :, None] - ou; rdo = rd[:, :, None] - od
        Pterm = np.zeros_like(Y4)
        Pterm[idx, idx] += s
        Pterm -= ruo * s[None, :, :] + s[:, None, :] * np.transpose(rdo, (1, 0, 2))
        dY4 = (np.einsum("ac,cbr->abr", h, Y4) + np.einsum("bd,adr->abr", h, Y4)
               + U * (vu[:, None, :] + vd[None, :, :]) * Y4 + U * Pterm)
        dy = -1j * dY4.reshape(L * L, r)
        hu = h + U * np.diag(nd); hd = h + U * np.diag(nu)
        g = self._gsel(None, y, ru, rd)
        Iu = U * (g["iiji"] - g["ijjj"])
        Id = U * (g["iiij"] - g["jijj"])
        dru = -1j * (hu @ ru - ru @ hu + Iu)
        drd = -1j * (hd @ rd - rd @ hd + Id)
        return [dru, drd, 0.0 * x, dy]

    def _rhs_fact(self, ru, rd, x, y, hu, hd):
        U = self._Ucur
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

    def _sum_rules(self, iters=20):
        """Exact constraints of the up-down pair matrix restored after each step by a
        positive rescaling of the pair amplitudes (keeps Y = y y^+ positive semidefinite):
          sum_b Y_(ab),(ab) = N_dn n^up_a,   sum_a Y_(ab),(ab) = N_up n^dn_b   (pair sum rules),
          U sum_a Y_(aa),(aa) = E_0 - Tr h (rho^up + rho^dn)                  (energy conservation).
        The first two use y -> (D_up (x) D_dn) y with positive diagonal D (Sinkhorn scaling);
        the third rescales the on-site amplitudes; the three are iterated jointly."""
        L = self.L
        Y4 = self.y.reshape(L, L, -1)
        nu = np.clip(np.real(np.diag(self.ru)), 1e-12, None)
        nd = np.clip(np.real(np.diag(self.rd)), 1e-12, None)
        du = np.ones(L); dd = np.ones(L); lam = 1.0
        diag = np.real(np.einsum("abr,abr->ab", Y4, Y4.conj()))
        h = self.hfun(self._tcur)
        Ekin = float(np.real(np.trace(h @ (self.ru + self.rd))))
        Dtarget = max((self.E0 - Ekin) / self._Ucur, 0.0) if self._Ucur > 0 else None
        idx = np.arange(L)
        for _ in range(iters):
            M = diag * (du[:, None] ** 2) * (dd[None, :] ** 2); M[idx, idx] *= lam ** 2
            du = du * np.sqrt(self.Nd * nu / np.maximum(M.sum(1), 1e-300))
            M = diag * (du[:, None] ** 2) * (dd[None, :] ** 2); M[idx, idx] *= lam ** 2
            dd = dd * np.sqrt(self.Nu * nd / np.maximum(M.sum(0), 1e-300))
            if Dtarget is not None:
                M = diag * (du[:, None] ** 2) * (dd[None, :] ** 2)
                lam = np.sqrt(Dtarget / max(np.trace(M), 1e-300))
        Y4 = Y4 * du[:, None, None] * dd[None, :, None]
        Y4[idx, idx] *= lam
        self.y = Y4.reshape(L * L, -1)

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

    def purify(self):
        """Positivity projection for the G1-G2 mode: X = A + G and Y = B + G are the
        hole-hole and particle-particle pair matrices, which must be positive
        semidefinite (N-representability).  Negative eigenvalues are removed and G is
        rebuilt as the average of the two projected estimates."""
        A, B = self._AB()
        G = 0.5 * (self.G + self.G.conj().T)
        out = []
        for M0 in (A + G, B + G):
            w, v = np.linalg.eigh(M0)
            out.append((v * np.clip(w, 0, None)) @ v.conj().T)
        self.G = 0.5 * ((out[0] - A) + (out[1] - B))

    def step(self, t, dt):
        self._tcur = t + dt
        if getattr(self, "sum_rules", False) and not hasattr(self, "E0"):
            U0 = self.U if self.Ufun is None else self.Ufun(t)
            nu0, nd0 = np.real(np.diag(self.ru)), np.real(np.diag(self.rd))
            self.E0 = float(np.real(np.trace(self.hfun(t) @ (self.ru + self.rd)))) + U0 * float(np.sum(nu0 * nd0 + self.dcorr()))
        s0 = self.state()
        k1 = self.rhs(t, s0)
        k2 = self.rhs(t + dt / 2, [a + dt / 2 * b for a, b in zip(s0, k1)])
        k3 = self.rhs(t + dt / 2, [a + dt / 2 * b for a, b in zip(s0, k2)])
        k4 = self.rhs(t + dt, [a + dt * b for a, b in zip(s0, k3)])
        self.set_state([a + dt / 6 * (b + 2 * c + 2 * d + e)
                        for a, b, c, d, e in zip(s0, k1, k2, k3, k4)])
        if self.mode == "g1g2" and getattr(self, "project", False):
            self.purify()
        if getattr(self, "sum_rules", False):
            self._sum_rules()

    # -- observables -----------------------------------------------------
    def pair_corr(self):
        """G_{ij,ij} = connected <n_i up n_j dn> for all i, j (factorized mode)."""
        L = self.L
        ru, rd = self.rhos()
        ubar = np.eye(L) - ru; dbar = np.eye(L) - rd
        out = 0.0
        for smp, (Mu, Md), wgt in ((self.x, (ubar, dbar), self.gw[0]), (self.y, (ru, rd), self.gw[1])):
            if wgt == 0:
                continue
            s4 = smp.reshape(L, L, -1)
            out = out + wgt * (np.real(np.einsum("ijr,ijr->ij", s4, s4.conj())) - np.real(np.outer(np.diag(Mu), np.diag(Md))))
        return out

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
