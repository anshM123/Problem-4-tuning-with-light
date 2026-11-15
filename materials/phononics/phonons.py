"""LiNbO3 (R3c, ferroelectric) from the MACE foundation potential medium-omat-0.

1. Experimental structure (Abrahams, Reddy & Bernstein, J. Phys. Chem. Solids 27, 997 (1966)):
   room-temperature lattice kept, internal coordinates relaxed with MACE under R3c symmetry.
   (A full MACE relaxation expands the cell by 3.4 %, the usual PBE-type overestimate, and
   softens the phonons; the measured lattice is used as for graphite.)
2. Harmonic (3x3x3 rhombohedral supercell, 270 atoms) and cubic (2x2x2, 80 atoms, pair
   cutoff 4.0 A) force constants with phono3py finite displacements (0.01 A).
3. Zone-centre modes: frequencies, irreducible representations and eigenvectors.

Writes (this directory): lno_relaxed.json, fc2.hdf5, fc3.hdf5, phono3py_params.yaml,
gamma_modes.npz.  Run time 20-60 min on 4 threads (1558 fc3 supercells of 80 atoms).
"""
import os, time, json
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
import torch
torch.set_num_threads(4)
from ase import Atoms
from ase.spacegroup import crystal
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from phonopy.structure.atoms import PhonopyAtoms
from phono3py import Phono3py
from phonopy import Phonopy
from mace.calculators import mace_mp

HERE = os.path.dirname(os.path.abspath(__file__))
calc = mace_mp(model="medium-omat-0", default_dtype="float64", device="cpu")

# ---- experimental structure, hexagonal setting of R3c (No. 161)
a_h, c_h = 5.14829, 13.8631
hexcell = crystal(["Nb", "Li", "O"],
                  basis=[(0, 0, 0), (0, 0, 0.2829), (0.0492, 0.3446, 0.0647)],
                  spacegroup=161, cellpar=[a_h, a_h, c_h, 90, 90, 120])
# rhombohedral primitive cell (10 atoms), polar axis along z
import spglib
lat, pos, num = spglib.standardize_cell((hexcell.cell[:], hexcell.get_scaled_positions(), hexcell.numbers),
                                        to_primitive=True)
atoms = Atoms(numbers=num, cell=lat, scaled_positions=pos, pbc=True)
print("primitive cell:", atoms.get_chemical_formula(), "V =", round(atoms.get_volume(), 3), "A^3")

t0 = time.time()
atoms.calc = calc
e0 = atoms.get_potential_energy()
atoms.set_constraint(FixSymmetry(atoms))
opt = BFGS(atoms, logfile=None)   # room-temperature experimental lattice, internal coordinates relaxed
opt.run(fmax=1e-4, steps=500)
atoms.set_constraint()
cp = atoms.cell.cellpar()
print(f"relaxed: a_r = {cp[0]:.4f} A, alpha = {cp[3]:.3f} deg, V = {atoms.get_volume():.3f} A^3 "
      f"(exp 106.08), dE = {atoms.get_potential_energy()-e0:.4f} eV, max F = "
      f"{np.abs(atoms.get_forces()).max():.1e} ({time.time()-t0:.0f}s)", flush=True)
# hexagonal lattice constants of the relaxed cell
a_r = cp[0]; al = np.radians(cp[3])
a_hex = 2 * a_r * np.sin(al / 2); c_hex = 3 * a_r * np.sqrt(1 - 4 / 3 * np.sin(al / 2) ** 2)
print(f"hexagonal: a = {a_hex:.4f} (exp 5.1483), c = {c_hex:.4f} (exp 13.8631)")
json.dump({"symbols": atoms.get_chemical_symbols(), "cell": atoms.cell[:].tolist(),
           "scaled_positions": atoms.get_scaled_positions().tolist(),
           "a_hex": a_hex, "c_hex": c_hex, "volume": atoms.get_volume()},
          open(os.path.join(HERE, "lno_relaxed.json"), "w"), indent=1)

unit = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.cell[:],
                    scaled_positions=atoms.get_scaled_positions())
ph3 = Phono3py(unit, supercell_matrix=[2, 2, 2], phonon_supercell_matrix=[3, 3, 3],
               primitive_matrix="P")
ph3.generate_displacements(distance=0.01, cutoff_pair_distance=4.0)
ph3.generate_fc2_displacements(distance=0.01)


def forces_of(sc):
    at = Atoms(sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True)
    at.calc = calc
    return at.get_forces()


sc2 = ph3.phonon_supercells_with_displacements
ph3.phonon_forces = np.array([forces_of(s) for s in sc2])
print(f"fc2 displacements: {len(sc2)} ({time.time()-t0:.0f}s)", flush=True)
allsc = ph3.supercells_with_displacements
n3 = sum(s is not None for s in allsc)
print(f"fc3 displacements: {n3}", flush=True)
F = []
for i, s in enumerate(allsc):
    F.append(np.zeros((len(ph3.supercell), 3)) if s is None else forces_of(s))
    if i % 500 == 0:
        print(f"  {i}/{len(allsc)} ({time.time()-t0:.0f}s)", flush=True)
ph3.forces = np.array(F)
ph3.produce_fc2(symmetrize_fc2=True, is_compact_fc=True)
ph3.produce_fc3(symmetrize_fc3r=True, is_compact_fc=True)
ph3.save(os.path.join(HERE, "phono3py_params.yaml"),
         settings={"force_sets": False, "displacements": False})   # fc2/fc3 are stored separately
from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5
write_fc2_to_hdf5(ph3.fc2, filename=os.path.join(HERE, "fc2.hdf5"))
write_fc3_to_hdf5(ph3.fc3, filename=os.path.join(HERE, "fc3.hdf5"))

# ---- zone-centre modes
phn = Phonopy(unit, supercell_matrix=[3, 3, 3], primitive_matrix="P")
phn.force_constants = ph3.fc2
phn.set_irreps([0, 0, 0], degeneracy_tolerance=1e-3)
irr = phn.irreps
labels = [""] * 30
try:
    for deg, lab in zip(irr._degenerate_sets, irr._ir_labels):
        for b in deg:
            labels[b] = lab or ""
except Exception as exc:  # labels are informative only
    print("irrep labels unavailable:", exc)
phn.run_qpoints([[0, 0, 0]], with_eigenvectors=True)
f = phn.qpoints.frequencies[0]
ev = phn.qpoints.eigenvectors[0]
masses = np.array(unit.masses)
print("Gamma modes (THz / cm-1):")
for b in range(30):
    print(f"  {b:2d} {labels[b]:>4s} {f[b]:7.3f} THz {f[b]*33.35641:7.1f} cm-1")
np.savez(os.path.join(HERE, "gamma_modes.npz"), freq=f, eigvec=ev.real, labels=np.array(labels),
         masses=masses, cell=atoms.cell[:], scaled_positions=atoms.get_scaled_positions(),
         numbers=atoms.get_atomic_numbers())
print(f"total {time.time()-t0:.0f}s")
