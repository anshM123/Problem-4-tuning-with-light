# Black phosphorus: experimental conditions used for the blind prediction

Only the experimental conditions and the analysis definitions of the experiments below were used
as inputs. Measured results (time constants, amplitudes, curves) were not consulted before the
predictions in `preregistered_predictions.json` were written.

## Experiments

1. H. Seiler, D. Zahn, M. Zacharias, P.-N. Hildebrandt, T. Vasileiadis, Y. W. Windsor, Y. Qi,
   C. Carbogno, C. Draxl, R. Ernstorfer, F. Caruso, "Accessing the anisotropic nonthermal phonon
   populations in black phosphorus", Nano Lett. 21, 6171 (2021), doi:10.1021/acs.nanolett.1c01786
   (diffuse scattering; primary target).
2. D. Zahn, P.-N. Hildebrandt, T. Vasileiadis, Y. W. Windsor, Y. Qi, H. Seiler, R. Ernstorfer,
   "Anisotropic nonequilibrium lattice dynamics of black phosphorus", Nano Lett. 20, 3728 (2020),
   doi:10.1021/acs.nanolett.0c00734 (Bragg-peak Debye-Waller analysis of the same setup and sample
   type; defines the mean-square-displacement analysis).

## Sample
- Bulk black phosphorus (HQ Graphene), mechanically exfoliated, free-standing film on a Cu TEM grid.
- Thickness 39 ± 5 nm (about 80 layers).
- Electron beam along the out-of-plane axis b; the probed reciprocal plane is (h0l), i.e. X–Γ–A.
- Axes (Cmce): a = 3.31 Å zigzag, b = 10.46 Å out of plane, c = 4.37 Å armchair.
- Base temperature 100 K.

## Pump
- 1.61 eV (770 nm), 50 fs, polarized along armchair (c).
- Incident fluence 9.8 ± 1.4 mJ/cm².
- Absorbed fluence from transfer-matrix optics with n(armchair, 800 nm) = 3.19 + 0.29i (the
  optical constant stated in ref. 1); the stated absorbed energy density is 380 ± 70 J/cm³.
- Repetition rate and angle of incidence not stated (normal incidence assumed).

## Probe and analysis
- 70 keV electrons; fits use a Gaussian time response of 150 fs FWHM (ref. 2); electron pulse
  ≈ 200 fs (ref. 1).
- Bragg peaks (400) (zigzag) and (004) (armchair): change of mean-square displacement along the
  peak direction from I(t)/I0 = exp(-ΔU·G²), fitted with a biexponential convolved with the 150 fs
  Gaussian (ref. 2).
- Diffuse intensity ΔI(Q,t) = I(Q,t) − I(Q,t<0) in the Brillouin-zone rectangle around (400), at
  the A point (biexponential: rise and slower relaxation) and the X point (single exponential)
  (ref. 1). Integration boxes not stated.

## Inputs to the pipeline
```
python run_material.py black_phosphorus --pump-eV 1.61 --fluence-mJcm2 9.8 --fwhm-fs 50 --T 100 \
    --thickness-nm 39 --nk 3.19+0.29j --pol 0,0,1 --tmax-ps 100 --beam "0 1 0" \
    --bragg "4 0 0;0 0 4;2 0 0;0 0 2" --diffuse-zone "4 0 0" --irf-fs 150 --qmax 10
```
Material inputs that are not experiment-specific: the crystal structure (Brown & Rundqvist,
Acta Cryst. 19, 684 (1965)) and the bulk band gap 0.33 eV (Keyes, Phys. Rev. 92, 580 (1953)),
which sets the scissor of the tight-binding bands and the split between carrier heat
(hν − E_g per pair) and energy stored across the gap (E_g per pair).

Special-point labels follow Setyawan and Curtarolo (base-centred orthorhombic): in the (h0l) plane
X = (0.55 a*, 0, 0), Z = (0, 0, c*/2), A = (0.55 a*, 0, c*/2) relative to a Bragg peak, with
a* = 2π/a and c* = 2π/c; Y and T appear in this plane at (a*, 0, 0) and (a*, 0, c*/2).
