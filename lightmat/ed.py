"""Exact diagonalization / exact time evolution for small Hubbard clusters.

H(t) = -sum_{<ij>,s} J_ij(t) c^+_is c_js + U sum_i n_iu n_id
States: fixed (N_up, N_dn); basis = (up configuration) x (down configuration).
Peierls phases enter through complex J_ij(t) = J_ij exp(i A(t) * (x_i - x_j)).
"""
from __future__ import annotations
import itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply, eigsh


def configs(L, N):
    out = []
    for occ in itertools.combinations(range(L), N):
        s = 0
        for o in occ:
            s |= 1 << o
        out.append(s)
    return np.array(out, dtype=np.int64)


def hop_matrix(L, N, bonds):
    """For each bond (i,j) return sparse matrix of c^+_i c_j on the N-particle space
    of a single spin species (fermionic sign from sites strictly between)."""
    cf = configs(L, N)
    index = {c: k for k, c in enumerate(cf)}
    mats = {}
    for (i, j) in bonds:
        for (a, b) in ((i, j), (j, i)):
            rows, cols, vals = [], [], []
            for k, c in enumerate(cf):
                if (c >> b) & 1 and not (c >> a) & 1:
                    c2 = c ^ (1 << b) ^ (1 << a)
                    lo, hi = min(a, b), max(a, b)
                    between = bin(c & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))).count("1")
                    sign = -1.0 if between % 2 else 1.0
                    rows.append(index[c2]); cols.append(k); vals.append(sign)
            mats[(a, b)] = sp.csr_matrix((vals, (rows, cols)), shape=(len(cf), len(cf)))
    return cf, mats


class Hubbard:
    def __init__(self, L, Nup, Ndn, bonds, U, J=None):
        self.L, self.Nup, self.Ndn, self.U = L, Nup, Ndn, U
        self.bonds = list(bonds)
        self.J = {b: 1.0 for b in self.bonds} if J is None else dict(J)
        self.cu, self.Mu = hop_matrix(L, Nup, self.bonds)
        self.cd, self.Md = hop_matrix(L, Ndn, self.bonds)
        self.du, self.dd = len(self.cu), len(self.cd)
        self.Iu = sp.identity(self.du, format="csr")
        self.Id = sp.identity(self.dd, format="csr")
        # occupations per site for each spin species
        self.nu = np.array([[(c >> i) & 1 for i in range(L)] for c in self.cu], float)
        self.nd = np.array([[(c >> i) & 1 for i in range(L)] for c in self.cd], float)
        dbl = (self.nu[:, None, :] * self.nd[None, :, :]).sum(-1).reshape(-1)
        self.Hint = sp.diags(U * dbl)

    def kinetic(self, phase=None):
        """phase: dict bond->complex factor multiplying J for c^+_i c_j (i,j)=bond."""
        Tu = sp.csr_matrix((self.du, self.du), dtype=complex)
        Td = sp.csr_matrix((self.dd, self.dd), dtype=complex)
        for (i, j) in self.bonds:
            f = self.J[(i, j)] * (1.0 if phase is None else phase[(i, j)])
            Tu = Tu - f * self.Mu[(i, j)] - np.conj(f) * self.Mu[(j, i)]
            Td = Td - f * self.Md[(i, j)] - np.conj(f) * self.Md[(j, i)]
        return sp.kron(Tu, self.Id, format="csr") + sp.kron(self.Iu, Td, format="csr")

    def H(self, phase=None):
        return self.kinetic(phase) + self.Hint

    def product_state(self, up_sites, dn_sites):
        su = sum(1 << i for i in up_sites); sd = sum(1 << i for i in dn_sites)
        iu = int(np.where(self.cu == su)[0][0]); idn = int(np.where(self.cd == sd)[0][0])
        psi = np.zeros(self.du * self.dd, complex)
        psi[iu * self.dd + idn] = 1.0
        return psi

    def ground_state(self, phase=None):
        w, v = eigsh(self.H(phase).real if phase is None else self.H(phase), k=1, which="SA")
        return w[0], v[:, 0].astype(complex)

    # observables -------------------------------------------------------
    def densities(self, psi):
        p = np.abs(psi.reshape(self.du, self.dd)) ** 2
        nup = p.sum(1) @ self.nu
        ndn = p.sum(0) @ self.nd
        return nup, ndn

    def double_occ(self, psi):
        p = np.abs(psi.reshape(self.du, self.dd)) ** 2
        return np.einsum("ab,ai,bi->i", p, self.nu, self.nd)

    def szsz(self, psi):
        """<S^z_i S^z_j> matrix."""
        p = np.abs(psi.reshape(self.du, self.dd)) ** 2
        sz = 0.5 * (self.nu[:, None, :] - self.nd[None, :, :])  # (du,dd,L)
        return np.einsum("ab,abi,abj->ij", p, sz, sz)

    def sz(self, psi):
        nu, nd = self.densities(psi)
        return 0.5 * (nu - nd)


def evolve(model, psi0, times, hfun=None, observe=None, substeps=4):
    """Propagate psi0 over the time grid. hfun(t)->sparse H (time dependent) or None
    for static H. observe(psi, t) -> dict of arrays recorded at each grid time."""
    out = []
    psi = psi0.copy()
    t_prev = times[0]
    H = model.H() if hfun is None else None
    out.append(observe(psi, t_prev))
    for t in times[1:]:
        if hfun is None:
            psi = expm_multiply(-1j * H * (t - t_prev), psi)
        else:
            dt = (t - t_prev) / substeps
            for k in range(substeps):
                tm = t_prev + (k + 0.5) * dt  # midpoint (exponential midpoint rule)
                psi = expm_multiply(-1j * hfun(tm) * dt, psi)
        t_prev = t
        out.append(observe(psi, t))
    return out


def chain_bonds(L, periodic=False):
    b = [(i, i + 1) for i in range(L - 1)]
    if periodic:
        b.append((L - 1, 0))
    return b


def ladder_bonds(Lx, periodic=True):
    """Sites (x, r) -> x + Lx*r, r in {0,1}. Returns leg bonds, rung bonds."""
    legs = []
    for r in (0, 1):
        for x in range(Lx if periodic else Lx - 1):
            legs.append((x + Lx * r, (x + 1) % Lx + Lx * r))
    rungs = [(x, x + Lx) for x in range(Lx)]
    return legs, rungs
