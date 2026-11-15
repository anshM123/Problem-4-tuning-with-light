"""Exact optical conductivity of the half-filled Hubbard ring (reference for Fig. 10 of
Schroedter et al.): Lanczos ground state plus a continued-fraction (Gauss-quadrature)
representation of the current spectral function A(w) = sum_n |<n|j|0>|^2 delta(w - E_n + E_0).

Ring (bc = pbc) or open chain (bc = obc) of L sites, N_up = N_dn = L/2, matrix-free H on
Psi[up config, down config]:
    H Psi = T Psi + Psi T^T + U * docc o Psi,
with true fermionic boundary conditions: periodic for odd N_sigma, antiperiodic for even
N_sigma (closed shell, non-degenerate ground state) for the ring.  The current is j = -dH/dA for a uniform
Peierls phase.  The real part of sigma uses exactly the same definition and Lorentzian
broadening as the time-domain fast method (kick A(t) = A0 theta(t), sigma = -j(w)/(A0 L)
with exp(-eta t)):
    Re sigma(w) = (1/L) sum_n (|j_n0|^2 / w_n) [eta/(eta^2+(w-w_n)^2) + eta/(eta^2+(w+w_n)^2)]
                  + (D/L) eta/(eta^2+w^2),
with the Drude part D = -<K> - 2 sum_n |j_n0|^2 / w_n (f-sum rule).
Also the exact (L = infinity) Bethe-ansatz Mott gap (Ovchinnikov / Lieb-Wu):
    Delta(U) = U - 4J + 8J int_0^inf J1(w) / (w (1 + exp(w U / 2J))) dw.
Usage: python chain_sigma_ed.py pbc|obc L U [U ...]   (writes data/chain_sigma_ed_{bc}_L{L}_U{U}.npz)
"""
import sys, os, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, LinearOperator
from scipy.special import j1
from scipy.integrate import quad
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lightmat.ed import hop_matrix

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def bethe_gap(U):
    f = lambda w: j1(w) / (w * (1.0 + np.exp(min(w * U / 2.0, 700.0))))
    return U - 4.0 + 8.0 * quad(f, 0, 200, limit=400)[0]


def sigma_from_poles(w, wn, an, L, eta, drude=0.0):
    """Lorentzian-broadened Re sigma(w) from current poles wn (> 0) and weights an = |j_n0|^2."""
    m = wn > 1e-8
    wn, an = wn[m], an[m]
    lor = lambda x: eta / (eta ** 2 + x ** 2)
    s = (an / wn)[None, :] * (lor(w[:, None] - wn[None, :]) + lor(w[:, None] + wn[None, :]))
    return (s.sum(1) + drude * lor(w)) / L


class Ring:
    def __init__(self, L, U, bc="pbc"):
        self.L, self.U = L, U
        N = L // 2
        phase = 1.0 if N % 2 == 1 else -1.0           # closed shell
        bonds = [(i, i + 1) for i in range(L - 1)] + ([(L - 1, 0)] if bc == "pbc" else [])
        cf, mats = hop_matrix(L, N, bonds)
        self.dim1 = len(cf)
        T = sp.csr_matrix((self.dim1, self.dim1))
        Jr = sp.csr_matrix((self.dim1, self.dim1))
        for (i, j) in bonds:
            f = phase if (i, j) == (L - 1, 0) else 1.0
            # hopping -J f c+_i c_j + h.c.; bond displacement x_i - x_j = -1 (i -> j = i+1)
            T = T - f * mats[(i, j)] - f * mats[(j, i)]
            # j = -dH/dA with phase exp(i A dx): j = i J dx f (c+_i c_j - c+_j c_i), dx = -1
            # store the real antisymmetric part Jr (j = i Jr up to a global sign)
            Jr = Jr - f * (mats[(i, j)] - mats[(j, i)])
        self.T, self.Jr = T.tocsr(), Jr.tocsr()
        cfu = cf.astype(np.int64)
        dbl = np.zeros((self.dim1, self.dim1))
        for i in range(L):
            b = ((cfu >> i) & 1).astype(float)
            dbl += np.outer(b, b)
        self.Ud = U * dbl
        self.n = self.dim1 ** 2

    def H(self, v):
        P = v.reshape(self.dim1, self.dim1)
        out = self.T @ P + (self.T @ P.T).T + self.Ud * P
        return out.reshape(-1)

    def K(self, v):
        P = v.reshape(self.dim1, self.dim1)
        return (self.T @ P + (self.T @ P.T).T).reshape(-1)

    def J(self, v):
        P = v.reshape(self.dim1, self.dim1)
        return (self.Jr @ P + (self.Jr @ P.T).T).reshape(-1)

    def ground_state(self):
        op = LinearOperator((self.n, self.n), matvec=self.H, dtype=float)
        rng = np.random.default_rng(1)
        E, V = eigsh(op, k=1, which="SA", tol=1e-11, ncv=16, v0=rng.standard_normal(self.n))
        return E[0], V[:, 0]

    def lanczos(self, phi, M):
        """Continued fraction of <phi| (z - H)^-1 |phi>: returns poles and weights."""
        nrm = np.linalg.norm(phi)
        q = phi / nrm; qold = np.zeros_like(q); b = 0.0
        al, be = [], []
        for k in range(M):
            w = self.H(q)
            a = q @ w
            w -= a * q + b * qold
            al.append(a)
            b = np.linalg.norm(w)
            if b < 1e-10:
                break
            be.append(b)
            qold, q = q, w / b
        n = len(al)
        Tm = np.diag(al) + np.diag(be[:n - 1], 1) + np.diag(be[:n - 1], -1)
        e, v = np.linalg.eigh(Tm)
        return e, nrm ** 2 * v[0] ** 2


def run(L, U, bc="pbc", M=400, eta=0.1):
    t0 = time.time()
    R = Ring(L, U, bc)
    E0, g = R.ground_state()
    phi = R.J(g)
    e, a = R.lanczos(phi, M)
    wn = e - E0
    Kexp = g @ R.K(g)
    docc = (g.reshape(R.dim1, R.dim1) ** 2 * R.Ud).sum() / U / L if U > 0 else np.nan
    drude = -Kexp - 2 * np.sum(a[wn > 1e-8] / wn[wn > 1e-8])
    w = np.linspace(0.02, 1.5 * U + 6, 700)
    sig = sigma_from_poles(w, wn, a, L, eta, drude=0.0)          # regular part only
    gap = bethe_gap(U)
    out = os.path.join(D, f"chain_sigma_ed_{bc}_L{L}_U{U:g}.npz")
    np.savez(out, w=w, sigma=sig, wn=wn, an=a, E0=E0, Ekin=Kexp, docc=docc, drude=drude,
             gap=gap, L=L, U=U, eta=eta, bc=bc)
    re = sig
    print(f"ED {bc} L={L} U={U:g}: E0/L={E0/L:.6f} <K>/L={Kexp/L:.5f} d={docc:.5f} Drude/L={drude/L:.2e} "
          f"lowest pole w/ weight>1e-3 of max: {wn[(wn>1e-6)&(a>1e-3*a.max())].min():.3f} "
          f"peak {w[np.argmax(re)]:.2f} Bethe gap {gap:.3f}  {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    bc = sys.argv[1]
    L = int(sys.argv[2])
    for U in sys.argv[3:]:
        run(L, float(U), bc)
