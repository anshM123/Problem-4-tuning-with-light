"""Exact time evolution of a small Hubbard-Holstein cluster with anharmonic phonons
(reference for lightmat/unified.py).

H(t) = sum_<ij>,s h_ij(t) c^+_is c_js + U sum_i n_iu n_id + sum_a w_a b^+_a b_a
     + sum_a x_a sum_i gamma^a_ii n_i + sum_a f_a(t) x_a + (1/6) sum_abc kappa_abc x_a x_b x_c,
x = (b + b^+)/sqrt2, phonon Fock spaces truncated at ncut[a] quanta (convergence checked by
the caller).  Electrons: lightmat.ed.Hubbard basis (fixed N_up, N_dn).
"""
import os, sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from lightmat.ed import Hubbard


class EPhED:
    def __init__(self, L, Nup, Ndn, bonds, J, U, omega, gamma, kappa, ncut, force=None, Afun=None):
        self.el = el = Hubbard(L, Nup, Ndn, bonds, U, J={b: J for b in bonds})
        self.bonds = bonds
        self.Afun = Afun
        self.force = force
        Q = len(omega)
        self.Q = Q
        de = el.du * el.dd
        ops = []
        for a in range(Q):
            n = ncut[a]
            b = sp.diags(np.sqrt(np.arange(1, n)), 1, format="csr")
            ops.append(b)
        def embed(op, a):
            mats = [sp.identity(ncut[k], format="csr") for k in range(Q)]
            mats[a] = op
            out = mats[0]
            for M in mats[1:]:
                out = sp.kron(out, M, format="csr")
            return out
        self.dph = int(np.prod(ncut))
        self.x = [embed((o + o.T) / np.sqrt(2), a) for a, o in enumerate(ops)]
        self.nb = [embed(o.T @ o, a) for a, o in enumerate(ops)]
        Iph = sp.identity(self.dph, format="csr")
        Ie = sp.identity(de, format="csr")
        Hph = sum(omega[a] * self.nb[a] for a in range(Q))
        for a in range(Q):
            for b in range(Q):
                for c in range(Q):
                    if kappa[a, b, c] != 0:
                        Hph = Hph + (kappa[a, b, c] / 6.0) * (self.x[a] @ self.x[b] @ self.x[c])
        nsite = el.nu[:, None, :] + el.nd[None, :, :]
        nsite = nsite.reshape(de, L)
        Heph = sp.csr_matrix((de * self.dph, de * self.dph))
        for a in range(Q):
            g = np.real(np.diag(gamma[a]))
            assert np.allclose(gamma[a], np.diag(np.diag(gamma[a]))), "site-diagonal coupling only"
            Heph = Heph + sp.kron(sp.diags(nsite @ g), self.x[a], format="csr")
        self.Hstat = (sp.kron(el.Hint, Iph) + sp.kron(Ie, Hph) + Heph).tocsr()
        self.Ie, self.Iph = Ie, Iph
        self.de = de
        self.Xop = [sp.kron(Ie, x, format="csr") for x in self.x]
        self.Nop = [sp.kron(Ie, n, format="csr") for n in self.nb]
        dbl = (el.nu[:, None, :] * el.nd[None, :, :]).reshape(de, L)
        self.Dop = [sp.kron(sp.diags(dbl[:, i]), Iph, format="csr") for i in range(L)]
        self.nop = [sp.kron(sp.diags(nsite[:, i]), Iph, format="csr") for i in range(L)]

    def H(self, t):
        ph = None
        if self.Afun is not None:
            A = self.Afun(t)
            ph = {b: np.exp(1j * A) for b in self.bonds}
        H = self.Hstat + sp.kron(self.el.kinetic(ph), self.Iph, format="csr")
        if self.force is not None:
            f = self.force(t)
            H = H + sum(f[a] * self.Xop[a] for a in range(self.Q))
        return H

    def initial(self, up, dn):
        pe = self.el.product_state(up, dn)
        v = np.zeros(self.dph); v[0] = 1.0
        return np.kron(pe, v).astype(complex)

    def run(self, psi, times, nsub=4):
        out = []
        t = times[0]
        for k, T in enumerate(times):
            if k > 0:
                dt = (T - t) / nsub
                for s in range(nsub):
                    psi = expm_multiply(-1j * self.H(t + (s + 0.5) * dt) * dt, psi)
                t = T
            ev = lambda O: float(np.real(np.vdot(psi, O @ psi)))
            out.append(dict(d=[ev(O) for O in self.Dop], n=[ev(O) for O in self.nop],
                            N=[ev(O) for O in self.Nop], X=[ev(O) for O in self.Xop],
                            E=ev(self.H(T))))
        return out
