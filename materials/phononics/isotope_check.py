"""Isotope (mass-disorder) scattering width of the driven A1(TO4) mode, natural abundances.

phono3py's Tamura model with the MACE harmonic force constants; adds to the three-phonon width.
Writes isotope_check.log (a few seconds).
"""
import os, io, contextlib
import numpy as np
import phono3py
from phono3py.api_isotope import Phono3pyIsotope

HERE = os.path.dirname(os.path.abspath(__file__))
ph3 = phono3py.load(os.path.join(HERE, "phono3py_params.yaml"), fc2_filename=os.path.join(HERE, "fc2.hdf5"),
                    fc3_filename=os.path.join(HERE, "fc3.hdf5"), produce_fc=False, log_level=0)
iIR = int(np.load(os.path.join(HERE, "anharmonic.npz"))["ia1"][-1])
lines = []
for mesh in (8, 12):
    iso = Phono3pyIsotope([mesh] * 3, ph3.phonon_primitive, sigmas=[None])
    iso.init_dynamical_matrix(ph3.fc2, ph3.phonon_supercell, ph3.phonon_primitive)
    with contextlib.redirect_stdout(io.StringIO()):
        iso.run([0])
    g = float(np.array(iso.gamma).reshape(-1)[iIR])
    lines.append(f"mesh {mesh}^3: isotope Gamma_TO4 = {g:.5f} THz half width -> FWHM {2*g*33.356:.3f} cm-1")
open(os.path.join(HERE, "isotope_check.log"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
