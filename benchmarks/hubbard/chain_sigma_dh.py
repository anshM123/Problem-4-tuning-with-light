"""Optical conductivity of the half-filled open Hubbard chain from the doublon-holon pair
closure of the particle-hole channel (lightmat.fluct_ph.DoublonHolonPairs).

Ground state: the method's own Hartree-Fock antiferromagnet of the same chain (no exact input).
From it: bond kinetic energies K_b (current vertex, three-site hops, f-sum rule) and bond
energies e_b (the exchange energy a doublon-holon pair removes from the background).  Linear
response to a vector-potential step, Re sigma(w) with the same exp(-eta t) Lorentzian
broadening (eta = 0.1 J) and frequency grid as chain_sigma.py.
Usage: python chain_sigma_dh.py [o1] [sw] [ring] L U [U ...]
  writes data/chain_sigma_dh_L{L}_U{U}.npz (open chain, second order in J/U: the method used);
  o1: leading order only (file tag _o1); ring: periodic ring without edge states, for the
  long-chain threshold (file tag _ring); sw: Hartree-Fock bond energies corrected by the
  zero-point spin fluctuations (RPA / linear spin waves) of the antiferromagnet (file tag _sw).
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lightmat.fluct_ph import DoublonHolonPairs, hf_antiferromagnet, local_energies

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


# Bond singlet weight 1/4 - <S_i.S_i+1> of the antiferromagnet: 1/2 for the Hartree-Fock (Neel)
# state; 1/4 + 1/4 + (1/2)(1 - 2/pi) with the zero-point transverse spin fluctuations of RPA
# (linear spin waves, Anderson 1952) on top of it.
SW_FACTOR = (0.25 + 0.25 + 0.5 * (1 - 2 / np.pi)) / 0.5


def run(L, U, order=2, ring=False, sw=False, eta=0.1, T=60.0, dt=0.02):
    t0 = time.time()
    h = np.zeros((L, L))
    for i in range(L - 1 + int(ring)):
        h[i, (i + 1) % L] = h[(i + 1) % L, i] = -1.0
    ru, rd, E = hf_antiferromagnet(h, U)
    if ring:
        b = np.arange(L); c = (b + 1) % L
        Kb = 2 * h[b, c] * (ru[b, c] + rd[b, c])
        d = np.diag(ru) * np.diag(rd)
        eb = Kb + U * (d[b] + d[c]) / 2
        eps = U * d + (Kb + np.roll(Kb, 1)) / 2
    else:
        Kb, eb, eps = local_energies(h, U, ru, rd)
    if sw:   # bond singlet weights, hence K_b and e_b, enhanced by the spin-wave zero-point factor
        Kb, eb, eps = Kb * SW_FACTOR, eb * SW_FACTOR, eps * SW_FACTOR
    m = DoublonHolonPairs(L, U, Kb, eb, order=order, eps=eps, periodic=ring)
    tt, j = m.current_response(T, dt)
    w = np.linspace(0.02, 1.5 * U + 6, 700)
    sig = np.array([-np.sum(j * np.exp((1j * x - eta) * tt)) * dt for x in w]) / L
    re = sig.real
    b = L // 2
    edge = U - 4 - 4 * eb[b] + Kb[b] if order == 2 else U - 4 - 2 * eps[b]
    tag = ("" if order == 2 else "_o1") + ("_sw" if sw else "") + ("_ring" if ring else "")
    np.savez(os.path.join(D, f"chain_sigma_dh{tag}_L{L}_U{U:g}.npz"), w=w, sigma=sig, t=tt, j=j,
             E0=E, Kb=Kb, eb=eb, eps=eps, edge=edge, L=L, U=U, eta=eta, order=order)
    print(f"dh{tag} L={L} U={U:g}: E_HF/L={E/L:.4f}  bulk band edge {edge:.3f}  "
          f"onset(10%) {w[np.argmax(re > 0.1 * re.max())]:.2f}  peak {w[np.argmax(re)]:.2f}  "
          f"{time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    order, ring, sw = 2, False, False
    while args[0] in ("o1", "ring", "sw"):
        order = 1 if args[0] == "o1" else order
        ring = ring or args[0] == "ring"
        sw = sw or args[0] == "sw"
        args = args[1:]
    L = int(args[0])
    for U in args[1:]:
        run(L, float(U), order, ring, sw)
