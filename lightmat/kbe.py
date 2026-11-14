"""Reference two-time Kadanoff-Baym solver (second-Born self-energy) for spin-symmetric
Hubbard models, used to benchmark accuracy and cost against the time-linear methods.

G^<_ij(t,t') = i<c^+_j(t') c_i(t)>,  G^>_ij(t,t') = -i<c_i(t) c^+_j(t')>   (one spin; spin-symmetric)
  i d/dt G^{<,>}(t,t') = h^HF(t) G^{<,>}(t,t') + I^{<,>}(t,t')
  I^<(t,t') = int_0^t ds [S^>-S^<](t,s) G^<(s,t') - int_0^t' ds S^<(t,s) [G^>-G^<](s,t')
  I^>(t,t') = int_0^t ds [S^>-S^<](t,s) G^>(s,t') - int_0^t' ds S^>(t,s) [G^>-G^<](s,t')
  i d/dt G^<(t,t) = [h^HF, G^<(t,t)] + I^<(t,t) + I^<(t,t)^+
Second Born (local U, opposite spin):  S^{<,>}_ij(t,t') = U^2 G^{<,>}_ij(t,t') G^{<,>}_ij(t,t') G^{>,<}_ji(t',t).
Uncorrelated initial state; Heun (predictor-corrector) time stepping, trapezoidal memory
integrals.  Cost O(N_t^3 L^3), memory O(N_t^2 L^2).
"""
import numpy as np


class KBE2B:
    def __init__(self, h, U, rho0, dt, nmax):
        self.h, self.U, self.dt = np.asarray(h, complex), U, dt
        L = h.shape[0]; self.L = L
        self.Gl = np.zeros((nmax + 1, nmax + 1, L, L), complex)
        self.Gg = np.zeros_like(self.Gl)
        self.Gl[0, 0] = 1j * rho0
        self.Gg[0, 0] = -1j * (np.eye(L) - rho0)
        self.Sl = np.zeros_like(self.Gl); self.Sg = np.zeros_like(self.Gl)
        self.n = 0
        self._sigma(0, 0)

    def _sigma(self, a, b):
        U2 = self.U ** 2
        for (x, y) in ((a, b), (b, a)):
            gl, gg = self.Gl[x, y], self.Gg[x, y]
            glr, ggr = self.Gl[y, x], self.Gg[y, x]
            self.Sl[x, y] = U2 * gl * gl * ggr.T
            self.Sg[x, y] = U2 * gg * gg * glr.T

    def _hHF(self, n):
        rho = np.real(np.diag(-1j * self.Gl[n, n]))
        return self.h + self.U * np.diag(rho)

    @staticmethod
    def _trap(w_n):
        w = np.ones(w_n); w[0] = w[-1] = 0.5
        if w_n == 1:
            w[0] = 0.0
        return w

    def _coll(self, n):
        """I^<(n, m), I^>(n, m) for m = 0..n (arrays (n+1, L, L))."""
        dt = self.dt
        SR = self.Sg[n, :n + 1] - self.Sl[n, :n + 1]             # (s, L, L)
        w = self._trap(n + 1) * dt
        A = np.einsum("s,sab,smbc->mac", w, SR, self.Gl[:n + 1, :n + 1])
        Bg = np.einsum("s,sab,smbc->mac", w, SR, self.Gg[:n + 1, :n + 1])
        GR = self.Gg[:n + 1, :n + 1] - self.Gl[:n + 1, :n + 1]     # (s, m)
        Il = np.empty((n + 1, self.L, self.L), complex); Ig = np.empty_like(Il)
        for m in range(n + 1):
            wm = self._trap(m + 1) * dt
            Il[m] = A[m] - np.einsum("s,sab,sbc->ac", wm, self.Sl[n, :m + 1], GR[:m + 1, m])
            Ig[m] = Bg[m] - np.einsum("s,sab,sbc->ac", wm, self.Sg[n, :m + 1], GR[:m + 1, m])
        return Il, Ig

    def _rhs(self, n):
        h = self._hHF(n)
        Il, Ig = self._coll(n)
        dl = -1j * (np.einsum("ab,mbc->mac", h, self.Gl[n, :n + 1]) + Il)
        dg = -1j * (np.einsum("ab,mbc->mac", h, self.Gg[n, :n + 1]) + Ig)
        Gd = self.Gl[n, n]
        dd = -1j * (h @ Gd - Gd @ h + Il[n] + Il[n].conj().T)
        return dl, dg, dd

    def step(self):
        n, dt = self.n, self.dt
        dl, dg, dd = self._rhs(n)
        # predictor (Euler)
        self.Gl[n + 1, :n + 1] = self.Gl[n, :n + 1] + dt * dl
        self.Gg[n + 1, :n + 1] = self.Gg[n, :n + 1] + dt * dg
        self.Gl[n + 1, n + 1] = self.Gl[n, n] + dt * dd
        self.Gg[n + 1, n + 1] = self.Gl[n + 1, n + 1] - 1j * np.eye(self.L)
        self._fill(n + 1)
        # corrector (Heun) on the new row
        dl2, dg2, dd2 = self._rhs(n + 1)
        self.Gl[n + 1, :n + 1] = self.Gl[n, :n + 1] + 0.5 * dt * (dl + dl2[:n + 1])
        self.Gg[n + 1, :n + 1] = self.Gg[n, :n + 1] + 0.5 * dt * (dg + dg2[:n + 1])
        self.Gl[n + 1, n + 1] = self.Gl[n, n] + 0.5 * dt * (dd + dd2)
        self.Gg[n + 1, n + 1] = self.Gl[n + 1, n + 1] - 1j * np.eye(self.L)
        self._fill(n + 1)
        self.n = n + 1

    def _fill(self, k):
        for m in range(k):
            self.Gl[m, k] = -self.Gl[k, m].conj().T
            self.Gg[m, k] = -self.Gg[k, m].conj().T
        for m in range(k + 1):
            self._sigma(k, m)

    def density(self, n=None):
        n = self.n if n is None else n
        return np.real(np.diag(-1j * self.Gl[n, n]))
