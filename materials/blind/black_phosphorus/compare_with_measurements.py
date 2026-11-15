"""Comparison of the preregistered black-phosphorus predictions with the published measurements
(reference_measurements.md, compiled after the predictions were committed in 9389d30).

    python materials/blind/black_phosphorus/compare_with_measurements.py

Inputs: preregistered_predictions.json (unchanged), results.json (the preregistered run, used only
for quantities that need the full traces, e.g. the normalized-MSD equalization time), and the
measured values transcribed below from reference_measurements.md. Writes comparison.json.

Point labels. The pipeline labels special points with the Setyawan-Curtarolo convention (ASE), in
which X lies along zigzag (a*) and A is the zone corner. Seiler et al. 2021 use another convention:
"no local minima arise in the conduction band along Gamma-X and Z-Q (armchair direction)" and
"emission of phonons with momenta q along the Gamma-A directions" = "along the zigzag direction".
So the paper's X is along armchair (l) and its A along zigzag (h). The preregistered file defines
every point by its coordinates, so the paper's points are matched by coordinates:
  paper A (zigzag, around (400)) <-> preregistered "X" = (4 +- 0.5, 0, 0)
  paper X (armchair, around (400)) <-> preregistered "Z" = (4, 0, +-0.5)
Both the literal-label and the coordinate-matched comparisons are reported.
"""
import json
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(os.path.join(HERE, "preregistered_predictions.json")))
R = json.load(open(os.path.join(HERE, "results.json")))
pr = P["predictions"]

# measured (Zahn et al. 2020 Table I, Fig. 2; Seiler et al. 2021 text and Fig. 2), see reference_measurements.md
M = {
    "zz_tau1": (0.58, 0.04), "zz_tau2": (20.0, 2.0), "zz_A1": (0.68e-3, 0.02e-3), "zz_A2": (0.85e-3, 0.03e-3),
    "ac_tau1": (0.48, 0.05), "ac_tau2": (22.0, 1.0), "ac_A1": (0.66e-3, 0.03e-3), "ac_A2": (1.84e-3, 0.03e-3),
    "zz_dU_inf": (1.53e-3, 0.04e-3), "ac_dU_inf": (2.50e-3, 0.04e-3), "aniso_inf": (1.63, 0.05),
    "zz_fast_frac": (0.44, None), "ac_fast_frac": (0.26, None),
    "I400_plateau_Z20": (0.911, 0.003), "I004_plateau_Z20": (0.910, 0.003),
    "I200_plateau_Z20": (0.982, 0.003), "I002_plateau_Z20": (0.976, 0.003),
    "I400_2ps_S21": (0.961, 0.003), "I004_2ps_S21": (0.977, 0.003),
    "I400_10ps_S21": (0.939, 0.003), "I004_10ps_S21": (0.953, 0.003),
    "I400_100ps_S21": (0.921, 0.003), "I004_100ps_S21": (0.924, 0.003),
    "diffA_rise": (1.7, 0.1), "diffA_decay": (30.0, 2.0), "diffX_tau": (14.3, 0.1),
    "diffA_peak_rel": (0.017, 0.002), "diffX_plateau_rel": (0.039, 0.002),
    "nonthermal_lifetime_Z20": (60.0, None), "thermal_S21": ("50-100", None),
    "final_T_Z20": (370.0, 50.0),
}

zz = pr["Bragg_400_zigzag"]; ac = pr["Bragg_004_armchair"]
fz = zz["biexp_fit_of_delta_MSD"]; fa = ac["biexp_fit_of_delta_MSD"]
lnI = lambda b, t: b["minus_ln_I_over_I0_at_times"][t]
d = pr["diffuse"]
# Z20 definition of the non-thermal lifetime: MSD normalized to its 100 ps value is the same in both
# directions (evaluated on the preregistered traces, threshold 2 % and 5 %)
t = np.array(R["traces"]["t_ps"])
br = {tuple(int(round(x)) for x in b["hkl"]): b for b in R["bragg"]}
uz = np.array(br[(4, 0, 0)]["trace_dU_along_G_A2"]); ua = np.array(br[(0, 0, 4)]["trace_dU_along_G_A2"])
k100 = int(np.argmin(abs(t - 100)))
diff = np.abs(uz / uz[k100] - ua / ua[k100])
def t_equal(th):
    bad = np.where((diff > th) & (t > 0.2))[0]
    return float(t[bad[-1] + 1]) if len(bad) and bad[-1] + 1 < len(t) else None

rows = []
def row(name, pred, meas, unit="", note="", posthoc=None, form=""):
    m, e = meas
    ratio = (pred / m) if isinstance(pred, (int, float)) and isinstance(m, (int, float)) and m else None
    rp = (posthoc / m) if isinstance(posthoc, (int, float)) and isinstance(m, (int, float)) and m else None
    rows.append(dict(observable=name, unit=unit, fit_form=form, preregistered=pred,
                     post_hoc_relabelled=posthoc, measured=m, measured_err=e,
                     ratio_pred_over_meas=ratio, ratio_posthoc_over_meas=rp, note=note))

BI = "biexp A1(1-e^-t/t1)+A2(1-e^-t/t2) (x) 150 fs Gaussian, as Z20"
row("dMSD zigzag tau1", fz["tau1_ps"], M["zz_tau1"], "ps", form=BI)
row("dMSD zigzag tau2", fz["tau2_ps"], M["zz_tau2"], "ps", form=BI)
row("dMSD armchair tau1", fa["tau1_ps"], M["ac_tau1"], "ps", form=BI)
row("dMSD armchair tau2", fa["tau2_ps"], M["ac_tau2"], "ps", form=BI)
row("dMSD zigzag A1", fz["A1"], M["zz_A1"], "A^2", form=BI)
row("dMSD zigzag A2", fz["A2"], M["zz_A2"], "A^2", form=BI)
row("dMSD armchair A1", fa["A1"], M["ac_A1"], "A^2", form=BI)
row("dMSD armchair A2", fa["A2"], M["ac_A2"], "A^2", form=BI)
row("dMSD zigzag at 100 ps", zz["delta_MSD_along_G_A2_at_times"]["100"], M["zz_dU_inf"], "A^2")
row("dMSD armchair at 100 ps", ac["delta_MSD_along_G_A2_at_times"]["100"], M["ac_dU_inf"], "A^2")
row("late anisotropy dU_ac/dU_zz (100 ps)", ac["delta_MSD_along_G_A2_at_times"]["100"] / zz["delta_MSD_along_G_A2_at_times"]["100"],
    M["aniso_inf"])
row("fast fraction A1/(A1+A2), zigzag", fz["A1"] / (fz["A1"] + fz["A2"]), M["zz_fast_frac"])
row("fast fraction A1/(A1+A2), armchair", fa["A1"] / (fa["A1"] + fa["A2"]), M["ac_fast_frac"])
for hkl, key, lab in (((4, 0, 0), "I400_plateau_Z20", "Bragg_400_zigzag"), ((0, 0, 4), "I004_plateau_Z20", "Bragg_004_armchair"),
                      ((2, 0, 0), "I200_plateau_Z20", "Bragg_200_zigzag"), ((0, 0, 2), "I002_plateau_Z20", "Bragg_002_armchair")):
    row(f"I/I0 {hkl} at 100 ps (Z20 plateau)", float(np.exp(-lnI(pr[lab], "100"))), M[key],
        note="ratio of -ln(I/I0): %.2f" % (lnI(pr[lab], "100") / -np.log(M[key][0])))
for tt in ("2", "10", "100"):
    for hkl, lab in (("400", "Bragg_400_zigzag"), ("004", "Bragg_004_armchair")):
        row(f"I/I0 ({hkl}) at {tt} ps (S21)", float(np.exp(-lnI(pr[lab], tt))), M[f"I{hkl}_{tt}ps_S21"],
            note="ratio of -ln(I/I0): %.2f" % (lnI(pr[lab], tt) / -np.log(M[f"I{hkl}_{tt}ps_S21"][0])))
# diffuse: literal labels vs coordinate-matched labels
Xp, Ap, Zp = d["X"], d["A"], d["Z"]
row("diffuse A (S21: zigzag point, (400) zone) rise", Ap["fit_biexponential"]["tau1_ps"], M["diffA_rise"], "ps",
    posthoc=Xp["fit_biexponential"]["tau1_ps"], form="biexp, rise then relaxation",
    note="preregistered column: point labelled A (SC corner (4.5,0,0.5)); post hoc: same point type by coordinates, (4.5,0,0)")
row("diffuse A (S21) relaxation", Ap["fit_biexponential"]["tau2_ps"], M["diffA_decay"], "ps",
    posthoc=Xp["fit_biexponential"]["tau2_ps"], form="biexp, rise then relaxation",
    note="preregistered A biexp has a negative first amplitude (no rise-then-decay shape)")
row("diffuse X (S21: armchair point, (400) zone) single-exp tau", Xp["fit_single_exponential"]["tau_ps"], M["diffX_tau"], "ps",
    posthoc=Zp["fit_single_exponential"]["tau_ps"], form="single exponential A(1-e^-t/tau)",
    note="preregistered column: point labelled X (zigzag, (4.5,0,0)); post hoc: armchair point (4,0,0.5), labelled Z in the file")
row("diffuse A (S21) peak dI/I0", max(Ap["dI_over_I0_at_times"].values()), M["diffA_peak_rel"],
    posthoc=max(Xp["dI_over_I0_at_times"].values()),
    note="not like-for-like: simulated one-phonon intensity only; measured I includes background (Bragg tails, "
         "multi-phonon, inelastic), which dilutes relative changes")
row("diffuse X (S21) plateau dI/I0", Xp["dI_over_I0_at_times"]["100"], M["diffX_plateau_rel"], posthoc=Zp["dI_over_I0_at_times"]["100"],
    note="not like-for-like (see above)")
row("non-thermal lifetime (normalized MSD equal in both directions)", t_equal(0.02), M["nonthermal_lifetime_Z20"], "ps",
    note="derived from the preregistered traces with Z20's definition, 2 %% threshold (5 %%: %s ps); preregistered "
         "lattice_equilibration_ps (branch-temperature spread < 10 %%) = %.1f ps"
         % (t_equal(0.05), pr["characteristic_times"]["lattice_equilibration_ps"]))
row("lattice equilibration (preregistered definition)", pr["characteristic_times"]["lattice_equilibration_ps"],
    (60.0, None), "ps", note="measured: ~60 ps (Z20), 50-100 ps (S21)")
row("final lattice temperature", P["predictions"]["branch_temperatures"]["T_acoustic_K"][-1], M["final_T_Z20"], "K",
    note="20 % of the absorbed energy stays across the gap (no recombination); S21 quotes 300 K")
row("absorbed energy density", P["simulation"]["absorbed"]["energy_density_J_cm3"], (380.0, 70.0), "J/cm^3",
    note="the paper's optical constant was an input")

audit = dict(
    axes="pipeline Cartesian x = a = 3.31 A = zigzag, y = b = out of plane, z = c = 4.37 A = armchair; "
         "(h00) -> U_xx = U11, (00l) -> U_zz = U33, identical to Z20",
    labels="pipeline: Setyawan-Curtarolo (ASE) labels, X = (0.55 a*,0,0) zigzag, Z = (0,0,c*/2) armchair, A = corner. "
           "S21 text: Gamma-X is armchair, Gamma-A is zigzag. The two conventions differ; matched by coordinates. "
           "(reference_measurements.md reads the S21 Fig. 2b inset as the opposite assignment; the S21 text, quoted "
           "in the module docstring, is used here.)",
    msd_definition="per-atom <u_a u_b> = hbar/(2 M N) sum (2n+1)/omega Re(e_a e_b*), projected on G-hat; "
                   "I/I0 = exp(-G.dU.G) as Z20 Eqs. 2-3. Cross-check with phonopy thermal displacements: "
                   "U_xx/U_yy/U_zz(100 K) = 3.19/4.55/4.74 e-3 A^2 (phonopy, same grid) vs 3.36/4.33/4.85 e-3 "
                   "(pipeline); converged 24x24x18 grid: 3.24/4.71/4.92 e-3; dU(100->304 K) = 3.83/6.50/6.89 e-3 "
                   "(phonopy, thermal) vs pipeline at 100 ps 3.87/6.29/7.38 e-3. No bookkeeping error.",
    equilibrium="harmonic MACE: dU(100->370 K) = 5.2e-3 (zigzag), 9.4e-3 (armchair) A^2 from phonopy thermal "
                "displacements; the measured asymptotic dMSD (1.53e-3, 2.50e-3) is 3.4-3.8x smaller than a fully "
                "thermal lattice at the experiment's own dT = 270 K. In the classical limit dU/dT = 3 hbar^2/(M k "
                "Theta^2) for a Debye solid; the measured slopes imply Theta_eff ~ 900 K (zigzag) and ~710 K "
                "(armchair), the MACE phonons 490 K and 370 K.",
    verdict="no convention or bookkeeping bug in run_material.py; the diffuse-point labels differ between conventions "
            "and are matched by coordinates; the amplitude excess is a physics/measurement question, not fixed here",
)
out = dict(rows=rows, audit=audit, preregistered_timestamp=P["timestamp_utc"], preregistered_sha1=P["pipeline"]["sha1"])
json.dump(out, open(os.path.join(HERE, "comparison.json"), "w"), indent=1)
for r in rows:
    f = lambda v: ("%.3g" % v) if isinstance(v, (int, float)) else str(v)
    print(f"{r['observable'][:62]:62s} pre {f(r['preregistered']):>9s}  post {f(r['post_hoc_relabelled']):>9s}  "
          f"meas {f(r['measured']):>7s}  ratio {f(r['ratio_pred_over_meas']):>6s} / {f(r['ratio_posthoc_over_meas'])}")
