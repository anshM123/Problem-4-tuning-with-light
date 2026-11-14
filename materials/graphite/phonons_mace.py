"""Harmonic and cubic force constants of graphene from the MACE-MP universal potential.
Writes fc2.hdf5 / fc3.hdf5 / phono3py_disp.yaml into this directory.
In-plane phonons of graphite and graphene agree closely; the weak interlayer coupling
mainly affects the lowest out-of-plane modes near Gamma.
"""
import os, time
import numpy as np
from ase import Atoms
from phonopy.structure.atoms import PhonopyAtoms
from phono3py import Phono3py
from mace.calculators import mace_mp

HERE = os.path.dirname(os.path.abspath(__file__))
a = 2.464  # graphite in-plane lattice constant (Angstrom)
c = 15.0
cell = [[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]]
unit = PhonopyAtoms(symbols=["C", "C"], cell=cell,
                    scaled_positions=[[1 / 3, 2 / 3, 0.5], [2 / 3, 1 / 3, 0.5]])
calc = mace_mp(model="medium-omat-0", default_dtype="float64", device="cpu")

ph3 = Phono3py(unit, supercell_matrix=[5, 5, 1], phonon_supercell_matrix=[8, 8, 1],
               primitive_matrix="auto")
ph3.generate_displacements(distance=0.01, cutoff_pair_distance=4.0)
ph3.generate_fc2_displacements(distance=0.01)


def forces(scells):
    out = []
    for sc in scells:
        at = Atoms(sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True)
        at.calc = calc
        out.append(at.get_forces())
    return np.array(out)


t0 = time.time()
sc2 = ph3.phonon_supercells_with_displacements
ph3.phonon_forces = forces(sc2)
print(f"fc2 displacements: {len(sc2)} ({time.time()-t0:.0f}s)", flush=True)
sc3 = [s for s in ph3.supercells_with_displacements if s is not None]
print(f"fc3 displacements: {len(sc3)}", flush=True)
# phono3py expects forces for every displacement entry (None -> zeros)
allsc = ph3.supercells_with_displacements
F = []
for s in allsc:
    if s is None:
        F.append(np.zeros((len(ph3.supercell), 3)))
    else:
        at = Atoms(s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True)
        at.calc = calc
        F.append(at.get_forces())
ph3.forces = np.array(F)
print(f"forces done ({time.time()-t0:.0f}s)", flush=True)
ph3.produce_fc2(symmetrize_fc2=True)
ph3.produce_fc3(symmetrize_fc3r=True)
ph3.save(os.path.join(HERE, "phono3py_params.yaml"))
from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5
write_fc2_to_hdf5(ph3.fc2, filename=os.path.join(HERE, "fc2.hdf5"))
write_fc3_to_hdf5(ph3.fc3, filename=os.path.join(HERE, "fc3.hdf5"))
# quick check: frequencies at Gamma, M, K (THz)
ph3.mesh_numbers = [12, 12, 1]
from phonopy import Phonopy
phn = Phonopy(unit, supercell_matrix=[8, 8, 1], primitive_matrix="auto")
phn.force_constants = ph3.fc2
for name, q in (("G", [0, 0, 0]), ("M", [0.5, 0, 0]), ("K", [1 / 3, 1 / 3, 0])):
    f = phn.get_frequencies(q)
    print(name, np.round(f * 4.135667, 1), "meV")
print(f"total {time.time()-t0:.0f}s")
