"""AB-stacked graphite (4 atoms): harmonic and cubic force constants from MACE (medium-omat-0)
at the experimental lattice (a = 2.464 A, c = 6.711 A)."""
import os, time
import numpy as np
from ase import Atoms
from phonopy.structure.atoms import PhonopyAtoms
from phono3py import Phono3py
from mace.calculators import mace_mp
HERE = os.path.dirname(os.path.abspath(__file__))
a, c = 2.464, 6.711
cell = [[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]]
# layer A at z=1/4: atoms (0,0) and (1/3,2/3); layer B at z=3/4: atoms (0,0) and (2/3,1/3)
unit = PhonopyAtoms(symbols=["C"] * 4, cell=cell,
                    scaled_positions=[[0, 0, 0.25], [1/3, 2/3, 0.25], [0, 0, 0.75], [2/3, 1/3, 0.75]])
calc = mace_mp(model="medium-omat-0", default_dtype="float64", device="cpu")
ph3 = Phono3py(unit, supercell_matrix=[4, 4, 2], phonon_supercell_matrix=[6, 6, 2], primitive_matrix="auto")
ph3.generate_displacements(distance=0.01, cutoff_pair_distance=3.45)
ph3.generate_fc2_displacements(distance=0.01)
def F(sc):
    at = Atoms(sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True); at.calc = calc
    return at.get_forces()
t0 = time.time()
ph3.phonon_forces = np.array([F(s) for s in ph3.phonon_supercells_with_displacements])
sc3 = ph3.supercells_with_displacements
print("fc3 displacements:", sum(s is not None for s in sc3), flush=True)
ph3.forces = np.array([F(s) if s is not None else np.zeros((len(ph3.supercell), 3)) for s in sc3])
print(f"forces {time.time()-t0:.0f}s", flush=True)
ph3.produce_fc2(symmetrize_fc2=True); ph3.produce_fc3(symmetrize_fc3r=True)
ph3.save(os.path.join(HERE, "phono3py_params.yaml"))
from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5
write_fc2_to_hdf5(ph3.fc2, filename=os.path.join(HERE, "fc2.hdf5"))
write_fc3_to_hdf5(ph3.fc3, filename=os.path.join(HERE, "fc3.hdf5"))
print("done", time.time() - t0)
