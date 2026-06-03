"""
Analyze a resolved-MUG Hartmann channel run.

Reads hartmann_mug.profile (rows "x y velx", with a '# Ha_pred = ...' header),
verifies the flow is developed (velx ~ independent of x), bins velx by the
wall-normal coordinate y, fits the analytic insulating-Hartmann shape to extract
an effective Ha, and compares to the predicted Ha = B0 a / sqrt(mu0 eta rho nu).

Usage: python analyze_mug.py <profile_file> [a=1.0] [label=Ha?]
Writes (next to the profile): <dir>/profile_binned.csv and overlay PNG; also
echoes Ha_fit, Ha_pred, and L2/Linf vs the fitted analytic profile.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analytic_hartmann as ah

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False

from scipy.optimize import curve_fit


def hartmann_shape(y, Ha, C, a=1.0):
    return C * (1.0 - np.cosh(Ha * y / a) / np.cosh(Ha))


def load_profile(path):
    Ha_pred = None
    with open(path) as f:
        for line in f:
            if "Ha_pred" in line:
                Ha_pred = float(line.split("=")[1])
                break
    data = np.loadtxt(path, comments="#")
    x, y, ux = data[:, 0], data[:, 1], data[:, 2]
    return x, y, ux, Ha_pred


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "mug", "hartmann_mug.profile")
    a = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    label = sys.argv[3] if len(sys.argv) > 3 else "MUG"
    x, y, ux, Ha_pred = load_profile(path)

    # developed-flow check: spread of velx across x at fixed y
    yb = np.round(y, 6)
    yu = np.unique(yb)
    u_of_y = np.array([ux[yb == yy].mean() for yy in yu])
    spread = np.array([ux[yb == yy].std() for yy in yu])
    umax = np.max(np.abs(u_of_y)) or 1.0
    dev = spread.max() / umax
    print(f"developed-flow check: max std(velx)/max|u| over x = {dev:.3e} "
          f"({'OK' if dev < 0.02 else 'NON-UNIFORM in x!'})")

    # sort by y
    o = np.argsort(yu); yu, u_of_y = yu[o], u_of_y[o]
    uc = np.interp(0.0, yu, u_of_y)
    un = u_of_y / uc

    # fit Ha (and amplitude) to the binned profile
    p0 = [Ha_pred if Ha_pred else 5.0, uc]
    try:
        popt, _ = curve_fit(lambda yy, Ha, C: hartmann_shape(yy, Ha, C, a),
                            yu, u_of_y, p0=p0, maxfev=20000)
        Ha_fit, C_fit = popt
        Ha_fit = abs(Ha_fit)
    except Exception as e:
        print("curve_fit failed:", e)
        Ha_fit, C_fit = float("nan"), uc

    # errors of the (normalized) profile vs analytic at Ha_fit and Ha_pred
    def errs(Ha):
        ana = ah.u_profile(yu, a, Ha); ana = ana / np.interp(0, yu, ana)
        d = un - ana
        l2 = np.sqrt(np.trapezoid(d**2, yu) / np.trapezoid(ana**2, yu))
        return l2, np.max(np.abs(d))//1e-12*1e-12 if False else np.max(np.abs(d))

    l2_fit, linf_fit = errs(Ha_fit)
    muc = np.trapezoid(u_of_y, yu) / (yu[-1] - yu[0]) / uc

    print(f"{label}: Ha_pred={Ha_pred:.4g}  Ha_fit={Ha_fit:.4g}  "
          f"ratio fit/pred={Ha_fit/Ha_pred if Ha_pred else float('nan'):.4g}")
    print(f"  normalized-profile vs analytic(Ha_fit): L2={l2_fit:.3e} Linf={linf_fit:.3e}")
    print(f"  u_mean/u_c: MUG={muc:.4f}  analytic(Ha_fit)={ah.u_mean_over_uc(Ha_fit):.4f}")

    # write binned profile CSV (for compare.py: Ha_<n>_profile.csv with y,ux)
    outdir = os.path.dirname(os.path.abspath(path))
    np.savetxt(os.path.join(outdir, "profile_binned.csv"),
               np.column_stack([yu, u_of_y]), delimiter=",", header="y,ux", comments="")
    # also write under the integer-Ha name compare.py looks for
    if Ha_pred:
        nm = os.path.join(outdir, f"Ha_{int(round(Ha_pred))}_profile.csv")
        np.savetxt(nm, np.column_stack([yu, u_of_y]), delimiter=",", header="y,ux", comments="")

    if HAVE_MPL:
        yy = np.linspace(-a, a, 400)
        ana = ah.u_profile(yy, a, Ha_fit); ana = ana / np.interp(0, yy, ana)
        plt.figure(figsize=(6, 4.2))
        plt.plot(ana, yy, "k-", lw=2, label=f"analytic (Ha={Ha_fit:.2f})")
        plt.plot(un, yu, "s", ms=3, color="C3", label=f"{label} (resolved)")
        plt.xlabel("u / u_centerline"); plt.ylabel("y / a")
        plt.title(f"Resolved MUG Hartmann channel, Ha_pred={Ha_pred:.2f}")
        plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
        png = os.path.join(outdir, "mug_overlay.png")
        plt.savefig(png, dpi=130); plt.close()
        print("  wrote", png)


if __name__ == "__main__":
    main()
