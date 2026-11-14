"""Electron-phonon coupling of graphene pi electrons to all phonon branches.

pi bands: nearest-neighbour tight binding, H_AB(k) = t sum_j exp(i k.d_j).
Coupling: bond-stretching modulation of the hopping, dt = beta * d_hat.(u_B - u_A),
with phonon eigenvectors and frequencies from the MACE force constants (phonopy
convention u_s = Q e_s exp(i q.r_s)/sqrt(M_s)).  beta is fixed by the DFT value of the
A1' coupling at K, <g^2_K>_F = 0.0994 eV^2 (Piscanec et al., PRL 93, 185503 (2004));
<g^2_Gamma>_F is then a prediction of the model.

Phonon generation rates for thermal electrons at temperature Te (mu = 0):
  dn_{nu q}/dt = A_{nu q}(Te) [ (n+1) exp(-hw/kTe) - n ]
  A_{nu q}(Te) = (2 pi/hbar) (2/N_k) sum_{k,m,n} |g^nu_mn(k,q)|^2
                 delta(e_m(k+q) - e_n(k) - hw) f_n(k)(1 - f_m(k+q)).
"""
from __future__ import annotations
import numpy as np

HBAR_eVs = 6.582119569e-16
HBAR_Js = 1.054571817e-34
AMU = 1.66053906660e-27
M_C = 12.011 * AMU
KB_eV = 8.617333262e-5
THz2eV = 4.135667696e-3  # h * 1 THz in eV (phonopy frequencies are nu in THz)

a = 2.464
A1 = np.array([a, 0.0]); A2 = np.array([-a / 2, a * np.sqrt(3) / 2])
B = 2 * np.pi * np.linalg.inv(np.array([A1, A2])).T          # rows b1, b2
TAU = {"A": (1 / 3) * A1 + (2 / 3) * A2, "B": (2 / 3) * A1 + (1 / 3) * A2}
# A -> B nearest-neighbour vectors
_cands = []
for n1 in (-1, 0, 1):
    for n2 in (-1, 0, 1):
        d = TAU["B"] + n1 * A1 + n2 * A2 - TAU["A"]
        if abs(np.linalg.norm(d) - a / np.sqrt(3)) < 1e-6:
            _cands.append(d)
D = np.array(_cands)                       # (3,2)
DHAT = D / np.linalg.norm(D, axis=1)[:, None]
T_HOP = -3.0                               # eV  (hbar v_F = 3|t| a_cc / 2 = 6.4 eV A)


def frac_to_cart(f):
    return np.asarray(f) @ B


def bands(kc):
    """kc: (...,2) Cartesian k. Returns energies (...,2) ascending and eigvecs (...,2,2)."""
    f = T_HOP * np.exp(1j * kc @ D.T).sum(-1)
    af = np.abs(f)
    ph = np.where(af > 1e-14, f / np.where(af > 1e-14, af, 1), 1.0)
    e = np.stack([-af, af], -1)
    # eigenvectors (A,B) components: lower (1, -conj(ph))/sqrt2, upper (1, conj(ph))/sqrt2
    v = np.empty(kc.shape[:-1] + (2, 2), complex)
    v[..., 0, 0] = 1 / np.sqrt(2); v[..., 1, 0] = -np.conj(ph) / np.sqrt(2)
    v[..., 0, 1] = 1 / np.sqrt(2); v[..., 1, 1] = np.conj(ph) / np.sqrt(2)
    return e, v


def coupling(kc, qc, evec, freq_THz, beta):
    """g_mn(k,q) for one phonon mode.
    kc: (Nk,2); qc: (2,); evec: (2 atoms, 3) complex (phonopy, atoms ordered A,B);
    freq_THz: phonon frequency nu (THz).  Returns g (Nk, 2 (m at k+q), 2 (n at k)) in eV."""
    w = 2 * np.pi * freq_THz * 1e12
    l = np.sqrt(HBAR_Js / (2 * M_C * w)) * 1e10            # Angstrom
    eA = evec[0, :2]; eB = evec[1, :2]
    pk = np.exp(1j * kc @ D.T)                             # (Nk,3)
    pkq = np.exp(1j * (kc + qc) @ D.T)
    dA = DHAT @ eA; dB = DHAT @ eB                          # (3,)
    M_AB = beta * (pkq * dB - pk * dA).sum(-1)             # <A,k+q|dH|B,k>
    M_BA = beta * (np.conj(pk) * dB - np.conj(pkq) * dA).sum(-1)  # <B,k+q|dH|A,k>
    _, vk = bands(kc)
    _, vkq = bands(kc + qc)
    Hm = np.zeros((len(kc), 2, 2), complex)
    Hm[:, 0, 1] = M_AB; Hm[:, 1, 0] = M_BA
    g = np.einsum("kam,kab,kbn->kmn", np.conj(vkq), Hm, vk)
    return l * g


def electron_energy(Te, nk=600):
    """Electronic energy per unit cell (eV, both spins) relative to T=0, mu=0."""
    f = (np.arange(nk) + 0.5) / nk
    kk = np.stack(np.meshgrid(f, f, indexing="ij"), -1).reshape(-1, 2) @ B
    e, _ = bands(kk)
    eu = e[:, 1]
    x = eu / (KB_eV * Te)
    occ = 1 / (np.exp(np.clip(x, -700, 700)) + 1)
    # electrons above + holes below (particle-hole symmetric): 2 * sum e f, times spin 2
    return 2 * 2 * np.mean(eu * occ)
