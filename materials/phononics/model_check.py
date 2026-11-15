"""Cross-check of the anharmonicity along Q_IR with two other MACE foundation models.

Same structure and eigenvector (medium-omat-0); frozen-phonon energies for |Q| <= 1.5 fitted
with the polynomial form of von Hoegen et al.  Writes model_check.log (about 1 min).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
import torch
torch.set_num_threads(4)
from ase import Atoms
from mace.calculators import mace_mp

HERE = os.path.dirname(os.path.abspath(__file__))
an = np.load(os.path.join(HERE, "anharmonic.npz")); g = np.load(os.path.join(HERE, "gamma_modes.npz"))
e4 = an["evecA1"][3]; m = g["masses"]
at0 = Atoms(numbers=g["numbers"], cell=g["cell"], scaled_positions=g["scaled_positions"], pbc=True)
qs = np.arange(-1.5, 1.5001, 0.1)
out = open(os.path.join(HERE, "model_check.log"), "w")
for model in ("medium-omat-0", "medium-mpa-0", "medium"):
    calc = mace_mp(model=model, default_dtype="float64", device="cpu")
    a = at0.copy(); a.calc = calc; E0 = a.get_potential_energy(); p0 = a.get_positions()
    V = []
    for q in qs:
        b = at0.copy(); b.calc = calc
        b.set_positions(p0 + q * e4 / np.sqrt(m)[:, None])
        V.append(b.get_potential_energy() - E0)
    A = np.vstack([qs ** 2 / 2, qs ** 3 / 3, qs ** 4 / 4, qs ** 5 / 5]).T
    w2, a3, a4, a5 = np.linalg.lstsq(A, np.array(V), rcond=None)[0]
    nu = np.sqrt(w2 * 1.602176634e-19 / (1e-20 * 1.66053906660e-27)) / (2 * np.pi) / 1e12
    s = f"{model:15s} nu {nu:.2f} THz  a3 {a3*1e3:.0f}  a4 {a4*1e3:.0f}  a5 {a5*1e3:.0f} meV"
    print(s, flush=True); out.write(s + "\n")
out.write("DFT (von Hoegen et al. 2018): nu 19 THz (measured value used there), a3 1568, a4 901, a5 7.1 meV\n")
out.close()
