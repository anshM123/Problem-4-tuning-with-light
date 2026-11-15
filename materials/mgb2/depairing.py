"""Static depairing of the MgB2 film by a uniform vector potential A (supercurrent state, 7 K).

Two independent calculations with the same inputs (sigma_n,a, 1/tau_a, dressed DOS nu_a, pairing
couplings of twoband.py):

 (a) the random-matrix Bogoliubov-de Gennes model of realtime.py, in thermal equilibrium with a
     static A: H_a = eps tau3 - Delta_a tau1 - c_a A tau0 (x) Y_a diagonalized, R = f(H), gaps
     self-consistent, j = sum_a (c_a/V_a) Tr(Y_a R_a) - chi_dia,a A;
 (b) Eilenberger theory with Born impurity scattering of rate gamma_a = hbar/tau_a and a uniform
     superfluid momentum (Doppler shift s = e v_a A cos(theta)), valid at any purity:
        g = w_s / sqrt(w_s^2 + D~^2), f = D~ / sqrt(w_s^2 + D~^2),  w_s = w~ + i s,
        w~ = w_n + (gamma/2) <g>,   D~ = Delta + (gamma/2) <f>,
        Delta_a = sum_b V_ab N_b [Delta_b I_b + 2 pi T sum_{n>=0} Re(<f_b> - Delta_b/sqrt(w_n^2 + Delta_b^2))],
        j_a = 2 e nu_a v_a 2 pi T sum_{n>=0} <cos(theta) Im g>,
     with the Fermi velocity fixed by the Drude weight, e^2 nu v^2 (2/3 in 3D, 1 in 2D) = sigma_n gamma/hbar.
     pi band: 3D (isotropic), sigma band: cylinder with A in plane (2D); both geometries are reported.
For the dirty pi band (b) reduces to Usadel/Abrikosov-Gor'kov depairing; for the nearly clean
sigma band it contains the Doppler (Landau) physics of the clean limit.
Prints Delta_pi, Delta_sigma and j versus A, and the depairing current density j_d = max j.
"""
import numpy as np
import twoband as tb
import realtime as rt

T = 7.0
E_CH = 1.602176634e-19
HBAR = 1.054571817e-34
VN = tb.V * tb.NR[None, :]
WC = rt.WC


# ------------------------------------------------------------------ (a) random-matrix BdG
def bdg_state(bands, A, D0, iters=200):
    D = np.array(D0, float)
    for _ in range(iters):
        F, J, Rs = [], 0.0, []
        for b in bands:
            N = b.N
            H = np.zeros((2 * N, 2 * N), complex)
            Y = 1j * b.K
            al = b.c * A
            H[:N, :N] = np.diag(b.eps) - al * Y
            H[N:, N:] = -np.diag(b.eps) - al * Y
            H[:N, N:] = -D[b.a] * np.eye(N)
            H[N:, :N] = -D[b.a] * np.eye(N)
            w, U = np.linalg.eigh(H)
            # R_ij = <psi_j^+ psi_i> = sum_k U_ik f(w_k) U_jk^*
            R = (U * rt.fermi(w, T)) @ U.conj().T
            Rs.append(R)
            F.append(np.sum(np.conj(np.diagonal(R[N:, :N]))[b.win]).real * b.de)
        Dn = VN @ np.array(F)
        if np.max(abs(Dn - D)) < 1e-6:
            D = Dn
            break
        D = 0.6 * Dn + 0.4 * D
    for b, R in zip(bands, Rs):
        J += b.current(R) - b.chi_dia * A
    return D, J


# ------------------------------------------------------------------ (b) Eilenberger, Born impurities
def eil_band(a, Delta, A, geom, nmax=3000, nang=48):
    """returns (delta F, j) for band a at gap Delta and vector potential A."""
    gam = tb.GAM[a]
    rhoD = rt.SIGN[a] * gam / tb.HBAR_meVs                     # S/(m s)
    nuJ = rt.NUD[a] / rt.MEV                                   # per J per m^3 per spin
    fac = 2.0 / 3.0 if geom == "3d" else 1.0
    v = np.sqrt(rhoD / (fac * E_CH ** 2 * nuJ))                # m/s
    p = E_CH * v * A / rt.MEV                                  # meV
    if geom == "3d":
        x, wq = np.polynomial.legendre.leggauss(nang); wq = wq / 2
    else:
        phi = (np.arange(nang) + 0.5) * 2 * np.pi / nang; x = np.cos(phi); wq = np.full(nang, 1.0 / nang)
    kT = tb.KB * T
    wn = np.pi * kT * (2 * np.arange(nmax) + 1)
    s = p * x[None, :]
    wt = wn.copy().astype(complex); Dt = np.full(nmax, Delta, complex)
    for _ in range(400):
        ws = wt[:, None] + 1j * s
        r = np.sqrt(ws ** 2 + Dt[:, None] ** 2)
        g = ws / r; f = Dt[:, None] / r
        gm = g @ wq; fm = f @ wq
        wt_n = wn + 0.5 * gam * gm; Dt_n = Delta + 0.5 * gam * fm
        if np.max(abs(wt_n - wt)) + np.max(abs(Dt_n - Dt)) < 1e-10:
            wt, Dt = wt_n, Dt_n
            break
        wt = 0.5 * wt + 0.5 * wt_n; Dt = 0.5 * Dt + 0.5 * Dt_n
    ws = wt[:, None] + 1j * s
    r = np.sqrt(ws ** 2 + Dt[:, None] ** 2)
    g = ws / r; f = Dt[:, None] / r
    dF = 2 * np.pi * kT * np.sum((f @ wq).real - Delta / np.sqrt(wn ** 2 + Delta ** 2))
    jj = 2 * E_CH * nuJ * v * 2 * np.pi * kT * rt.MEV * np.sum((g.imag * x[None, :]) @ wq)
    return dF, jj, v


def eil_state(A, D0, geoms=("3d", "2d"), iters=300):
    D = np.array(D0, float)
    for _ in range(iters):
        F, J = [], 0.0
        for a in range(2):
            dF, jj, _ = eil_band(a, D[a], A, geoms[a]) if D[a] > 1e-6 else (0.0, 0.0, 0)
            F.append(D[a] * tb._I(max(D[a], 1e-9), T, WC) + dF)
            J += jj
        Dn = np.maximum(VN @ np.array(F), 0)
        if np.max(abs(Dn - D)) < 1e-6:
            D = Dn
            break
        D = 0.6 * Dn + 0.4 * D
    return D, abs(J)


if __name__ == "__main__":
    import realtime_runs as rr
    bands = [rt.Band(0, rr.NB["N"][0], rr.NB["W"][0], rr.NB["seed"]), rt.Band(1, rr.NB["N"][1], rr.NB["W"][1], rr.NB["seed"] + 1)]
    D0 = rt.gap_eq_grid(bands, T)
    for a, g in ((0, "3d"), (1, "2d")):
        print(f"band {a}: Fermi velocity from the Drude weight ({g}) = {eil_band(a, 1.0, 1e-9, g)[2]:.3e} m/s")
    As = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30]) * 1e-8
    out = []
    De = np.array(tb.gaps(T, WC)); Db = D0.copy(); De3 = De.copy()
    print("   A (V s/m) | BdG: 2Dpi 2Dsig (THz)  j (A/cm^2) | Eilenberger pi 3D, sigma 2D | pi 3D, sigma 3D")
    for A in As:
        Db, jb = bdg_state(bands, A, np.maximum(Db, 0.05))
        De, je = eil_state(A, np.maximum(De, 0.05), ("3d", "2d"))
        De3, je3 = eil_state(A, np.maximum(De3, 0.05), ("3d", "3d"))
        out.append((A, *Db, abs(jb), *De, je, *De3, je3))
        print(f"  {A:9.2e} | {2*Db[0]/tb.THZ:.3f} {2*Db[1]/tb.THZ:.3f}  {abs(jb)/1e4:9.3e} | {2*De[0]/tb.THZ:.3f} {2*De[1]/tb.THZ:.3f}"
              f" {je/1e4:9.3e} | {2*De3[0]/tb.THZ:.3f} {2*De3[1]/tb.THZ:.3f} {je3/1e4:9.3e}", flush=True)
    out = np.array(out)
    np.save(rt.os.path.join(rt.HERE, "rt_out", "depairing.npy"), out)
    print(f"depairing current j_d = max j: BdG {out[:,3].max()/1e4:.2e} A/cm^2 at A = {out[np.argmax(out[:,3]),0]:.1e};"
          f" Eilenberger(3D/2D) {out[:,6].max()/1e4:.2e} at {out[np.argmax(out[:,6]),0]:.1e}; Eilenberger(3D/3D) {out[:,9].max()/1e4:.2e}")
    print(f"screening sheet current 2E/Z0 per film thickness: {2*150e5/376.73/20e-9/1e4:.2e} A/cm^2 at 150 kV/cm, "
          f"{2*21e5/376.73/20e-9/1e4:.2e} at 21 kV/cm")
