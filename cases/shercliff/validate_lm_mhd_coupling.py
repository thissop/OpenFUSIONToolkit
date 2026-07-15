"""
Gap-4 (coupling half): validate the LM_MHD plasma-blanket COUPLING formulas
against exact analytic references.

`coupling.py` (blanket -> plasma feedback: Biot-Savart loop field, current
integrals, regime metrics) and `plasma_interface.py` (plasma -> blanket sheath
loads + liquid-surface response) are formula-level, "not a solver." So the
correctness question is: do the formulas match their exact analytic references?
This harness answers it decisively. Each check has an independent closed form.

NOTE on scope: this closes "the coupling FORMULAS are correct." The full coupled
loop (a resolved MUG axisymmetric blanket current -> field perturbation ->
TokaMaker plasma equilibrium and back) is a separate integration test that needs
the coupled machinery + an axisymmetric MUG blanket run; it is NOT claimed here.

Run:  python3 validate_lm_mhd_coupling.py
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]


def _load(mod):
    p = _ROOT / "src/python/OpenFUSIONToolkit/LM_MHD" / f"{mod}.py"
    spec = importlib.util.spec_from_file_location(f"lm_{mod}", p)
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


MU0 = 4.0e-7 * np.pi
_results = []


def check(name, got, ref, rel=1e-6):
    got = np.asarray(got, float); ref = np.asarray(ref, float)
    err = np.max(np.abs(got - ref) / (np.abs(ref) + 1e-300))
    ok = err <= rel
    _results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: max rel.err {err:.2e} (tol {rel:.0e})")


def numerical_loop_field(a, z0, I, R, z, nseg=20000):
    """Independent Biot-Savart: discretize the loop, sum dB at probe (R,0,z)."""
    phi = np.linspace(0, 2 * np.pi, nseg, endpoint=False)
    dphi = 2 * np.pi / nseg
    lx, ly, lz = a * np.cos(phi), a * np.sin(phi), np.full_like(phi, z0)
    dlx, dly = -a * np.sin(phi) * dphi, a * np.cos(phi) * dphi   # dl vector, dlz=0
    sx, sy, sz = R - lx, -ly, z - lz                             # probe - loop
    s3 = (sx * sx + sy * sy + sz * sz) ** 1.5
    # dl x s
    cx = dly * sz - 0.0 * sy
    cy = 0.0 * sx - dlx * sz
    cz = dlx * sy - dly * sx
    k = MU0 * I / (4.0 * np.pi)
    Bx = k * np.sum(cx / s3); By = k * np.sum(cy / s3); Bz = k * np.sum(cz / s3)
    return Bx, By, Bz   # probe on x-axis: B_R = Bx, B_Z = Bz, By ~ 0


def main():
    cp = _load("coupling")
    pi = _load("plasma_interface")

    print("=== coupling.py: Biot-Savart circular loop field ===")
    a, z0, I = 1.3, 0.2, 1.0e5
    # (1) on-axis vs exact  B_z = mu0 I a^2 / (2 (a^2+zeta^2)^{3/2})
    zs = np.array([0.0, 0.5, 1.0, 2.0]) + z0
    probes = np.column_stack([np.zeros_like(zs), zs])
    fld = cp.circular_loop_poloidal_field(a, z0, I, probes)
    zeta = zs - z0
    exact_bz = MU0 * I * a**2 / (2.0 * (a**2 + zeta**2) ** 1.5)
    check("on-axis B_z vs exact", fld[:, 1], exact_bz, rel=1e-10)
    check("on-axis B_R == 0", fld[:, 0], np.zeros_like(zeta), rel=1e-12)
    check("loop-center B_z == mu0 I/(2a)", fld[0, 1], MU0 * I / (2 * a), rel=1e-10)
    # (2) off-axis vs independent numerical Biot-Savart
    offs = [(0.5, 0.0), (0.9, 0.4), (1.8, -0.6), (0.3, 1.1)]
    mods, nums = [], []
    for R, zrel in offs:
        z = z0 + zrel
        f = cp.circular_loop_poloidal_field(a, z0, I, np.array([[R, z]]))[0]
        Bx, By, Bz = numerical_loop_field(a, z0, I, R, z)
        mods += [f[0], f[1]]; nums += [Bx, Bz]
    check("off-axis (B_R,B_Z) vs numerical Biot-Savart", mods, nums, rel=2e-4)

    print("=== coupling.py: dimensionless regime metrics (PbLi-like) ===")
    sigma, rho, nu, B, U, L = 7.7e5, 9400.0, 1.1e-7, 4.0, 0.1, 0.1
    check("Hartmann number", cp.hartmann_number(B, L, sigma, rho, nu),
          B * L * np.sqrt(sigma / (rho * nu)))
    check("Alfven speed", cp.alfven_speed(B, rho), B / np.sqrt(MU0 * rho))
    check("magnetic Reynolds", cp.magnetic_reynolds_number(U, L, sigma), MU0 * sigma * U * L)
    check("interaction parameter", cp.interaction_parameter(B, L, U, sigma, rho),
          sigma * B**2 * L / (rho * U))

    print("=== plasma_interface.py: sheath + surface response ===")
    # Bohm sound speed: c_s = sqrt(e T_e / m_i), deuterium, T_e=100 eV
    e = 1.602176634e-19; mD = 2.0 * 1.67262192369e-27
    check("Bohm sound speed (100 eV D)", pi.bohm_sound_speed(100.0, ion_mass=mD),
          np.sqrt(e * 100.0 / mD), rel=1e-3)
    # semi-infinite constant-flux surface temperature: DeltaT = 2 q sqrt(alpha t/pi)/k
    q, t, k, rho_s, cp_s = 5.0e6, 2.0, 16.0, 9400.0, 190.0
    alpha = k / (rho_s * cp_s)
    check("semi-infinite const-flux surface DeltaT",
          pi.semi_infinite_constant_flux_temperature_rise(q, t, k, rho_s, cp_s),
          2.0 * q * np.sqrt(alpha * t / np.pi) / k)

    print("=== end-to-end sanity: blanket current -> field perturbation ===")
    # a 0.1 m radius blanket loop carrying 1 kA at 2 m minor-axis distance from a
    # B0=4 T plasma edge: feedback ratio delta_b/B0 should be tiny (weak coupling)
    f = cp.circular_loop_poloidal_field(0.1, 0.0, 1.0e3, np.array([[2.0, 0.0]]))[0]
    dB = np.hypot(*f)
    print(f"  |delta_B| at 2 m from a 1 kA / 0.1 m blanket loop = {dB:.3e} T "
          f"-> feedback delta_B/B0(4T) = {dB/4.0:.2e} (weak, as expected)")

    print(f"\nSUMMARY: {sum(_results)}/{len(_results)} analytic checks PASS")
    return 0 if all(_results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
