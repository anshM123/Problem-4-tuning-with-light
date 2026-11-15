"""Collects the run_material.py output for black phosphorus into preregistered_predictions.json.

    python materials/blind/black_phosphorus/write_predictions.py

Reads results.json (central prediction) and, when present, results_eph_x0.5.json /
results_eph_x2.json (electron-phonon coupling scaled by 0.5 and 2, the stated uncertainty of the
universal tight-binding coupling) and results_recomb.json (carrier recombination within 10 ps).
"""
import datetime
import hashlib
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..")


def load(name):
    p = os.path.join(HERE, name)
    return json.load(open(p)) if os.path.exists(p) else None


def pick(res):
    if res is None:
        return None
    br = {tuple(int(round(x)) for x in b["hkl"]): b for b in res["bragg"]}
    out = {}
    for hkl, lab in (((4, 0, 0), "Bragg_400_zigzag"), ((0, 0, 4), "Bragg_004_armchair"),
                     ((2, 0, 0), "Bragg_200_zigzag"), ((0, 0, 2), "Bragg_002_armchair")):
        b = br[hkl]
        out[lab] = dict(G_invA=b["G_invA"],
                        delta_MSD_along_G_A2_at_times=dict(zip(map(str, res["debye_waller"]["times_ps"]), b["dU_along_G_A2"])),
                        minus_ln_I_over_I0_at_times=dict(zip(map(str, res["debye_waller"]["times_ps"]), b["neg_ln_I_over_I0"])),
                        biexp_fit_of_delta_MSD=b["fit_dU_biexp"], biexp_fit_of_minus_ln_I=b["fit_neg_ln_I_biexp"])
    d = res["diffuse"]
    out["diffuse"] = {}
    for lab in d:
        out["diffuse"][lab] = dict(Q_example_invA=d[lab]["example_Q_invA"], n_equivalent_Q=d[lab]["n_Q_points"],
                                   grid_offset_invA=d[lab]["grid_offset_invA"],
                                   dI_over_I0_at_times=dict(zip(map(str, res["debye_waller"]["times_ps"]), d[lab]["dI_over_I0_at_times"])),
                                   fit_single_exponential=d[lab]["fit_single"], fit_biexponential=d[lab]["fit_biexp"],
                                   fit_rise_times_decay=d[lab]["fit_rise_decay"],
                                   rise_63_ps=d[lab]["rise_63_ps"], peak_ps=d[lab]["peak_ps"],
                                   one_phonon_weight_by_branch=d[lab]["branch_weights"])
    ct = res["characteristic_times"]
    out["characteristic_times"] = {k: v for k, v in ct.items() if k != "hot_phonon_fit_rise_decay"}
    out["characteristic_times"]["hot_phonon_fit_rise_decay"] = ct["hot_phonon_fit_rise_decay"]
    out["electron_temperature_peak_K"] = res["electrons"]["Te_peak_K"]
    out["branch_temperatures"] = {k: res["branches"][k] for k in ("times_ps", "Te_K", "T_optical_K", "T_acoustic_K")}
    return out


R = load("results.json")
pipe = os.path.join(ROOT, "run_material.py")
sha = hashlib.sha1(open(pipe, "rb").read()).hexdigest()
pred = dict(
    material="black phosphorus, bulk (39 nm free-standing exfoliated film)",
    experiment_citation=[
        "H. Seiler et al., Accessing the anisotropic nonthermal phonon populations in black phosphorus, "
        "Nano Lett. 21, 6171 (2021), doi:10.1021/acs.nanolett.1c01786",
        "D. Zahn et al., Anisotropic nonequilibrium lattice dynamics of black phosphorus, Nano Lett. 20, 3728 (2020), "
        "doi:10.1021/acs.nanolett.0c00734"],
    conditions=dict(thickness_nm=39, free_standing=True, base_temperature_K=100, pump_photon_energy_eV=1.61,
                    pump_fwhm_fs=50, pump_polarization="armchair (c)", incident_fluence_mJcm2=9.8,
                    refractive_index_armchair_800nm="3.19+0.29i (stated in Seiler et al.)",
                    electron_beam="along b (out of plane); probed plane (h0l)", instrument_response_fwhm_fs=150,
                    see="conditions.md"),
    timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    pipeline=dict(script="run_material.py", sha1=sha,
                  command=("python run_material.py black_phosphorus --pump-eV 1.61 --fluence-mJcm2 9.8 --fwhm-fs 50 "
                           "--T 100 --thickness-nm 39 --nk 3.19+0.29j --pol 0,0,1 --tmax-ps 100 --beam \"0 1 0\" "
                           "--bragg \"4 0 0;0 0 4;2 0 0;0 0 2\" --diffuse-zone \"4 0 0\" --irf-fs 150 --qmax 10"),
                  phonons="MACE medium-omat-0 force constants (fc2, fc3), full-grid three-phonon kernel",
                  electron_phonon="universal sp3 Slater-Koster tight binding (Harrison 1980), scissor to 0.33 eV",
                  no_fitting_to_this_experiment=True),
    simulation=dict(absorbed=R["optics"], electrons=R["electrons"], phonons=R["phonons"],
                    energy_check=R["energy_check"], warnings=R["warnings"]),
    definitions=dict(
        delta_MSD="change of the atom-averaged mean-square displacement projected on the Bragg vector G, "
                  "Delta U_G(t) = G^.Delta U(t).G^ (A^2), equal to -ln(I/I0)/|G|^2 for black phosphorus "
                  "(one Wyckoff site); zigzag = (h00), armchair = (00l)",
        biexp="Delta(t) = A1 (1 - exp(-t/tau1)) + A2 (1 - exp(-t/tau2)) for t > 0, convolved with a 150 fs FWHM "
              "Gaussian; tau1 < tau2; the simulated curve (50 fs pump) is first convolved with the 141 fs probe so that "
              "the total response is 150 fs",
        diffuse="one-phonon intensity at Q = q + G in the (h0l) plane, in the zone around (400), at the "
                "symmetry-equivalent points closest to (400); dI/I0 = I(t)/I(t<0) - 1 (time constants are the same "
                "for Delta I = I(t) - I(t<0)); single exponential: A (1 - exp(-t/tau)); biexponential as above "
                "(a rise and a slower relaxation appear as one positive and one negative amplitude); "
                "rise x decay: A (1 - exp(-t/tau_r)) exp(-t/tau_d) + C (1 - exp(-t/tau_d)); all convolved with 150 fs",
        points="Setyawan-Curtarolo labels of the base-centred orthorhombic zone, relative to (400): "
               "X = (4 +- 0.55, 0, 0), Z = (4, 0, +-0.5), A = (4 +- 0.55, 0, +-0.5), Y = (4 +- 1, 0, 0), "
               "T = (4 +- 1, 0, +-0.5) in conventional reciprocal-lattice units; X uses the nearest q-grid point "
               "(0.5 a* from the Bragg peak); the labels X1 and A1 are equivalent to X and A (identical traces); Y and T sit "
               "on glide-forbidden positions (h odd) of the (h0l) plane"),
    predictions=pick(R),
    sensitivity=dict(
        eph_coupling_x0p5=pick(load("results_eph_x0.5.json")),
        eph_coupling_x2=pick(load("results_eph_x2.json")),
        note="the universal tight-binding coupling is uncertain by about a factor of 2 in |g|^2; these runs bracket "
             "the central prediction"),
    expected_accuracy=(
        "Phonon frequencies and anharmonicity from the MACE foundation model: frequencies within ~5 %, three-phonon "
        "rates within a factor ~1.5-2 (graphite: lifetimes 2x shorter than DFT). Electron-phonon coupling from "
        "Harrison's universal hoppings: within a factor ~2 in |g|^2; the tight-binding bands of black phosphorus have "
        "a 3.8 eV gap before the scissor, so the band curvature near the edges, and hence the carrier phase space, is "
        "the least reliable input. Absorbed energy: from the stated optical constant (+-15 %). No recombination: the "
        "gap energy (0.33/1.61 = 20 % of the absorbed energy) stays in the carriers, so late-time Debye-Waller "
        "amplitudes are a lower bound if carriers recombine within 100 ps. Time constants are expected within a "
        "factor ~2, amplitudes within ~30 %."),
)
json.dump(pred, open(os.path.join(HERE, "preregistered_predictions.json"), "w"), indent=1)
print("written", os.path.join(HERE, "preregistered_predictions.json"))
