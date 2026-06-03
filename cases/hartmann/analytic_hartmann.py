"""
Analytic reference for fully-developed Hartmann channel flow (insulating walls).

Referee solution for the three-way verification in
docs/dev/lm_mhd_verification_results.md. PERSONAL VERIFICATION exercise reproducing
the classic Hartmann (1937) result; no capability claim about MUG.

Problem
-------
Fully-developed, pressure-driven, steady, incompressible flow between two parallel
insulating plates at y = -a and y = +a. Streamwise velocity u(y) in x, uniform
transverse applied field B0 in y, insulating walls (net streamwise current = 0).

    Ha = B0 * a * sqrt(sigma / (rho * nu))

Closed form, driven by constant pressure gradient G = -dp/dx, with amplitude
C = G / (sigma * B0**2):

    u(y) = C * [ 1 - cosh(Ha * y/a) / cosh(Ha) ]                               (1)

Centerline velocity:

    u_c = u(0) = C * [ 1 - 1/cosh(Ha) ]                                        (2)

Mean velocity (per unit span, integrating (1) over [-a, a]):

    u_mean = C * [ 1 - tanh(Ha)/Ha ]                                           (3)

EXACT mean/centerline ratio (this is the unambiguous referee scalar):

    u_mean/u_c = [ 1 - tanh(Ha)/Ha ] / [ 1 - 1/cosh(Ha) ]                      (4)

    limits:  Ha -> 0   =>  (Ha^2/3)/(Ha^2/2) = 2/3   (Poiseuille)
             Ha -> inf =>  1/1 = 1                    (plug core)

NOTE on the task wording: the prompt states "u_mean/u_c -> 1 - tanh(Ha)/Ha" and
"Ha large -> 1/Ha (plug core)". Strictly, 1 - tanh(Ha)/Ha is u_mean/C (mean over the
driving amplitude, eq. 3), NOT u_mean/u_c; and the "1/Ha" is the FLOW-RATE suppression,
not u_mean/u_c (which -> 1 in a plug core). This module implements the exact relations
(2)-(5) and verifies the two physically unambiguous limits: u_mean/u_c -> 2/3 (Ha->0)
and Q* ~ 1/Ha (Ha large). Logged in the overnight log.

Hartmann-normalized flow rate (the classic 1/Ha suppression for fixed G):

    Q*(Ha) = (1/Ha) * [ 1 - tanh(Ha)/Ha ]                                      (5)

References
----------
Hartmann, J. (1937). Hg-Dynamics I. Mat.-Fys. Medd. 15(6).
Mueller & Buehler, "Magnetofluiddynamics in Channels and Containers" (2001), Ch. 3.
"""

import numpy as np

try:
    _trapz = np.trapezoid          # numpy >= 2.0
except AttributeError:             # pragma: no cover
    _trapz = np.trapz


def u_profile(y, a, Ha, C=1.0):
    """Hartmann velocity profile u(y), eq. (1). C = G/(sigma*B0**2) sets the scale;
    the profile SHAPE is independent of C. Robust at Ha->0 (Poiseuille limit)."""
    y = np.asarray(y, dtype=float)
    Ha = float(Ha)
    if Ha < 1e-9:
        # Taylor of (1): 1 - cosh(Ha y/a)/cosh(Ha) -> (Ha^2/2)(1 - (y/a)^2).
        # Keep the Ha^2/2 factor so amplitude is consistent with eq.(2),(3) limits.
        return C * (Ha * Ha / 2.0) * (1.0 - (y / a) ** 2)
    return C * (1.0 - np.cosh(Ha * y / a) / np.cosh(Ha))


def u_centerline(a, Ha, C=1.0):
    """Centerline velocity, eq. (2)."""
    return u_profile(0.0, a, Ha, C)


def u_mean(Ha, a=1.0, C=1.0):
    """Mean velocity, eq. (3): C*[1 - tanh(Ha)/Ha]. Ha->0 limit -> C*Ha^2/3."""
    Ha = float(Ha)
    if Ha < 1e-6:
        return C * (Ha * Ha / 3.0)
    return C * (1.0 - np.tanh(Ha) / Ha)


def u_mean_over_uc(Ha):
    """EXACT mean/centerline ratio, eq. (4). -> 2/3 (Ha->0), -> 1 (Ha->inf)."""
    Ha = float(Ha)
    if Ha < 1e-4:
        # numerator -> Ha^2/3, denominator -> Ha^2/2  => 2/3
        return 2.0 / 3.0
    return (1.0 - np.tanh(Ha) / Ha) / (1.0 - 1.0 / np.cosh(Ha))


def u_mean_over_uc_numeric(a, Ha, npts=40001):
    """Mean/centerline ratio by direct numerical integration of the profile.
    Independent cross-check of eq. (4); used in the self-test."""
    y = np.linspace(-a, a, npts)
    um = _trapz(u_profile(y, a, Ha), y) / (2.0 * a)
    return um / u_centerline(a, Ha)


def flow_rate(Ha, a=1.0, C=1.0):
    """Dimensional flow rate per unit span Q = int_{-a}^{a} u dy = 2a*u_mean."""
    return 2.0 * a * u_mean(Ha, a, C)


def flow_rate_star(Ha):
    """Hartmann-normalized flow rate Q*, eq. (5). Q* ~ 1/Ha at large Ha."""
    Ha = float(Ha)
    if Ha < 1e-6:
        return Ha / 3.0
    return (1.0 / Ha) * (1.0 - np.tanh(Ha) / Ha)


def hartmann_layer_thickness(a, Ha):
    """Hartmann boundary-layer thickness ~ a/Ha."""
    return a / float(Ha)


def _selftest():
    ok = True
    print("Hartmann analytic self-test")
    print("-" * 64)

    # 1) Poiseuille limit (Ha->0): u_mean/u_c -> 2/3, by exact (4) and numeric
    r_exact = u_mean_over_uc(1e-5)
    r_num = u_mean_over_uc_numeric(1.0, 1e-3)
    print(f"  Ha->0 : u_mean/u_c exact={r_exact:.6f} numeric={r_num:.6f}  (expect 0.666667)")
    if abs(r_exact - 2/3) > 1e-3 or abs(r_num - 2/3) > 3e-3:
        print("    FAIL: Poiseuille 2/3 not recovered"); ok = False

    # 2) small-Ha profile is parabolic
    y = np.linspace(-1, 1, 201)
    u = u_profile(y, 1.0, 1e-3); par = u.max() * (1 - y**2)
    rel = np.linalg.norm(u - par) / np.linalg.norm(par)
    print(f"  Ha->0 : profile vs parabola rel L2 = {rel:.3e}  (expect < 1e-2)")
    if rel > 1e-2:
        print("    FAIL: small-Ha profile not parabolic"); ok = False

    # 3) plug-core limit (Ha large): u_mean/u_c -> 1
    r_big = u_mean_over_uc(200.0); r_big_n = u_mean_over_uc_numeric(1.0, 200.0)
    print(f"  Ha=200: u_mean/u_c exact={r_big:.6f} numeric={r_big_n:.6f}  (expect -> 1)")
    if r_big < 0.985 or r_big_n < 0.985:
        print("    FAIL: high-Ha core not plug-like"); ok = False

    # 4) flow-rate Q* ~ 1/Ha at large Ha
    for Ha in (50.0, 100.0, 200.0):
        q = flow_rate_star(Ha); approx = 1.0 / Ha
        rel = abs(q - approx) / approx
        print(f"  Ha={Ha:5.0f}: Q*={q:.6e}  1/Ha={approx:.6e}  rel={rel:.3e}")
        if rel > 0.05:
            print("    FAIL: Q* not ~ 1/Ha"); ok = False

    # 5) exact (4) vs numeric integration agree across the sweep
    print("  cross-check eq.(4) vs numeric over sweep:")
    for Ha in (1, 5, 10, 20, 50):
        e = u_mean_over_uc(Ha); n = u_mean_over_uc_numeric(1.0, Ha)
        rel = abs(e - n) / n
        flag = "ok" if rel < 1e-4 else "FAIL"
        if rel >= 1e-4:
            ok = False
        print(f"    Ha={Ha:3d}: exact={e:.6f} numeric={n:.6f} rel={rel:.2e} [{flag}]")

    print("-" * 64)
    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if _selftest() else 1)
