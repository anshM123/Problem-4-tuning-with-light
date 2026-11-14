"""Electron-phonon generation coefficients A_{q nu}(Te) for AB graphite (4 atoms, 12 branches).

Electrons: pi bands of each layer (nearest-neighbour tight binding, as in ../graphite);
interlayer hopping is neglected for the coupling.  Each layer couples to the in-plane
displacements of its own two atoms through the bond-length modulation beta d_hat.(u_B - u_A),
with the layer's own nearest-neighbour vectors and the phonopy eigenvector components of
that layer.  A_{q nu} is the sum of the two layer contributions.  Couplings are calibrated
exactly as for the monolayer (DFT <g^2_K>, then GW scaling in dynamics).  Grid points within
two spacings of Gamma, K, K' are cell-averaged (NSUB x NSUB in-plane sub-points).
Output: eph_rates.npz (qfrac, freq, Te, A[Nq, 12, nT] in 1/ps).
"""
import os, sys, time
import numpy as np, h5py
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "graphite"))
from eph_tb import bands, B as B2, HBAR_eVs, HBAR_Js, M_C, KB_eV, THz2eV
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

HERE = os.path.dirname(os.path.abspath(__file__))
NSUB = 4; Nk = int(sys.argv[1]) if len(sys.argv) > 1 else 450
SIG = 0.03; BETA2 = 0.0994 / 0.00454
a, c = 2.464, 6.711
A1 = np.array([a, 0.0]); A2 = np.array([-a / 2, a * np.sqrt(3) / 2])
cell = [[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]]
pos = np.array([[0, 0], [1/3, 2/3], [0, 0], [2/3, 1/3]])          # fractional in-plane
layers = [(0, 1), (2, 3)]
unit = PhonopyAtoms(symbols=["C"] * 4, cell=cell,
                    scaled_positions=[[0, 0, 0.25], [1/3, 2/3, 0.25], [0, 0, 0.75], [2/3, 1/3, 0.75]])
ph = Phonopy(unit, supercell_matrix=[6, 6, 2])
ph.force_constants = h5py.File(os.path.join(HERE, "fc2.hdf5"))["force_constants"][:]
P = np.load(os.path.join(HERE, "phph.npz"))
qf3, freq = P["qfrac"], P["freq"]
Te = np.geomspace(300, 40000, 32); kT = KB_eV * Te


def nn_vectors(iA, iB):
    rA = pos[iA] @ np.array([A1, A2]); out = []
    for n1 in (-1, 0, 1):
        for n2 in (-1, 0, 1):
            d = pos[iB] @ np.array([A1, A2]) + n1 * A1 + n2 * A2 - rA
            if abs(np.linalg.norm(d) - a / np.sqrt(3)) < 1e-6:
                out.append(d)
    return np.array(out)


DL = [nn_vectors(*l) for l in layers]


def bands_layer(kc, D):
    f = -3.0 * np.exp(1j * kc @ D.T).sum(-1)
    af = np.abs(f); ph_ = np.where(af > 1e-14, f / np.where(af > 1e-14, af, 1), 1.0)
    e = np.stack([-af, af], -1)
    v = np.empty(kc.shape[:-1] + (2, 2), complex)
    v[..., 0, 0] = 1 / np.sqrt(2); v[..., 1, 0] = -np.conj(ph_) / np.sqrt(2)
    v[..., 0, 1] = 1 / np.sqrt(2); v[..., 1, 1] = np.conj(ph_) / np.sqrt(2)
    return e, v


def g_layer(kc, qc, eA, eB, D, fTHz):
    w = 2 * np.pi * fTHz * 1e12
    l = np.sqrt(HBAR_Js / (2 * M_C * w)) * 1e10
    Dh = D / np.linalg.norm(D, axis=1)[:, None]
    pk = np.exp(1j * kc @ D.T); pkq = np.exp(1j * (kc + qc) @ D.T)
    dA = Dh @ eA; dB = Dh @ eB
    MAB = (pkq * dB - pk * dA).sum(-1)
    MBA = (np.conj(pk) * dB - np.conj(pkq) * dA).sum(-1)
    _, vk = bands_layer(kc, D); _, vkq = bands_layer(kc + qc, D)
    Hm = np.zeros((len(kc), 2, 2), complex); Hm[:, 0, 1] = MAB; Hm[:, 1, 0] = MBA
    return l * np.einsum("kam,kab,kbn->kmn", np.conj(vkq), Hm, vk)


f = (np.arange(Nk) + 0.5) / Nk
kc = np.stack(np.meshgrid(f, f, indexing="ij"), -1).reshape(-1, 2) @ B2
EK = [bands_layer(kc, D)[0] for D in DL]


def A_at(qfrac3_list):
    """Average of A over the given 3D q-points (sub-sampling); returns (nb, nT)."""
    ph.run_qpoints(qfrac3_list, with_eigenvectors=True)
    fs, evs = ph.qpoints.frequencies, ph.qpoints.eigenvectors
    nb = fs.shape[1]
    acc = np.zeros((nb, len(Te)))
    for isub, q3 in enumerate(qfrac3_list):
        qc = np.array(q3[:2]) @ B2
        for il, (iA, iB) in enumerate(layers):
            ek = EK[il]; ekq = bands_layer(kc + qc, DL[il])[0]
            for nu in range(nb):
                hw = fs[isub, nu] * THz2eV
                if hw < 1e-3:
                    continue
                e = evs[isub][:, nu].reshape(4, 3)
                g = g_layer(kc, qc, e[iA, :2], e[iB, :2], DL[il], fs[isub, nu])
                g2 = BETA2 * np.abs(g) ** 2
                for m in range(2):
                    for n in range(2):
                        de = ekq[:, m] - ek[:, n] - hw
                        s_ = np.abs(de) < 4 * SIG
                        if not s_.any():
                            continue
                        wgt = g2[s_, m, n] * np.exp(-0.5 * (de[s_] / SIG) ** 2) / (np.sqrt(2 * np.pi) * SIG)
                        en = ek[s_, n][:, None] / kT[None]; em = ekq[s_, m][:, None] / kT[None]
                        occ = 1 / (np.exp(np.clip(en, -700, 700)) + 1) * (1 - 1 / (np.exp(np.clip(em, -700, 700)) + 1))
                        acc[nu] += (wgt[:, None] * occ).sum(0)
    return acc * (2 * np.pi / HBAR_eVs) * 2 / len(kc) * 1e-12 / len(qfrac3_list)


Nq1 = int(round(np.sqrt(len(qf3) / 2)))
M2 = np.array([[1, 0], [-0.5, np.sqrt(3) / 2]])
def near(q, q0, r):
    return min(np.linalg.norm((q - np.array(q0) - [s1, s2]) @ M2) for s1 in (-1, 0, 1) for s2 in (-1, 0, 1)) <= r + 1e-9
offs = (np.arange(NSUB) + 0.5) / NSUB - 0.5
def task(iq):
    q = qf3[iq]
    if any(near(q[:2], c0, 2.0 / Nq1) for c0 in ([0, 0], [1/3, 1/3], [2/3, 2/3])):
        subs = [[q[0] + o1 / Nq1, q[1] + o2 / Nq1, q[2]] for o1 in offs for o2 in offs]
    else:
        subs = [list(q)]
    return A_at(subs)


if __name__ == "__main__":
    from multiprocessing import Pool
    t0 = time.time()
    with Pool(int(os.environ.get("NPROC", "20"))) as pool:
        res = pool.map(task, range(len(qf3)), chunksize=1)
    A = np.array(res)
    np.savez(os.path.join(HERE, "eph_rates.npz"), qfrac=qf3, freq=freq, Te=Te, A=A)
    print("done", time.time() - t0)
