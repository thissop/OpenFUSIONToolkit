"""
Three-way Hartmann verification: analytic vs mhdFoam (vs MUG, if available).

Reads OpenFOAM sampled profiles, computes the analytic referee profile on the same
y-grid, normalizes by centerline velocity, and reports L2/Linf profile errors plus
the discriminating scalar u_mean/u_c and the flow-rate ratio. Writes per-Ha CSVs,
overlay PNGs, a flow-rate-vs-Ha PNG, and error_table.md.

Honest scope: personal verification exercise reproducing classic Hartmann flow; not a
MUG capability claim. Drag term OFF throughout (resolved channel).
"""
import os
import glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, HERE)
import analytic_hartmann as ah

RESULTS = os.path.join(HERE, "results")
OPENFOAM = os.path.join(HERE, "openfoam")
MUG = os.path.join(HERE, "mug")
HA_LIST = [1, 5, 10, 20, 50]
A = 1.0  # channel half-width

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception as e:  # pragma: no cover
    HAVE_MPL = False
    print("matplotlib unavailable:", e)


def trapz(y, x):
    f = getattr(np, "trapezoid", None) or np.trapz
    return f(y, x)


from extract_profile import get_profile as _vtk_profile


def read_openfoam_profile(ha):
    """Return (y, ux, label) for this Ha via foamToVTK + meshio cell-centroid
    extraction (the OF v1912 sampling functionObject is broken: 'sha1' IOstream)."""
    case = os.path.join(OPENFOAM, f"Ha_{ha}")
    if not os.path.isdir(case):
        return None
    res = _vtk_profile(case, xtarget=10.0)
    if res is None:
        return None
    y, ux = res
    return y, ux, f"{case}/VTK"


def read_mug_profile(ha):
    """Return (y, ux) from a MUG-exported profile CSV if present, else None.
    Expected file: cases/hartmann/mug/Ha_<ha>_profile.csv with columns y,ux (header allowed)."""
    f = os.path.join(MUG, f"Ha_{ha}_profile.csv")
    if not os.path.exists(f):
        return None
    try:
        data = np.loadtxt(f, delimiter=",", skiprows=1)
    except Exception:
        data = np.loadtxt(f)
    if data.ndim != 2 or data.shape[1] < 2:
        return None
    y, ux = data[:, 0], data[:, 1]
    order = np.argsort(y)
    return y[order], ux[order], f


def norm_by_centerline(y, u):
    """Normalize a profile by its centerline value u(y=0) (interpolated)."""
    uc = np.interp(0.0, y, u)
    if abs(uc) < 1e-30:
        return u.copy(), 1.0
    return u / uc, uc


def mean_over_uc(y, u):
    """u_mean/u_c using trapezoid integration over [-a, a]."""
    uc = np.interp(0.0, y, u)
    um = trapz(u, y) / (y[-1] - y[0])
    return um / uc


def main():
    os.makedirs(RESULTS, exist_ok=True)
    rows = []
    flow = {"Ha": [], "ana": [], "of": [], "mug": []}

    for ha in HA_LIST:
        # analytic on a fine grid for plotting, and on the OF grid for error
        of = read_openfoam_profile(ha)
        mug = read_mug_profile(ha)

        # reference grid: use OF grid if available else uniform
        if of is not None:
            y_of, ux_of, of_file = of
        else:
            y_of, ux_of, of_file = None, None, None

        y_plot = np.linspace(-A, A, 400)
        u_ana_plot = ah.u_profile(y_plot, A, ha)
        u_ana_plot_n, _ = norm_by_centerline(y_plot, u_ana_plot)

        # analytic scalars
        muc_ana = ah.u_mean_over_uc(ha)
        flow["Ha"].append(ha)
        flow["ana"].append(ah.flow_rate_star(ha))

        row = {"Ha": ha, "muc_ana": muc_ana}

        # --- OpenFOAM error ---
        if of is not None:
            u_ana_on_of = ah.u_profile(y_of, A, ha)
            ana_n, _ = norm_by_centerline(y_of, u_ana_on_of)
            of_n, uc_of = norm_by_centerline(y_of, ux_of)
            diff = of_n - ana_n
            l2 = np.sqrt(trapz(diff**2, y_of) / trapz(ana_n**2, y_of))
            linf = np.max(np.abs(diff))
            muc_of = mean_over_uc(y_of, ux_of)
            row.update(of_l2=l2, of_linf=linf, of_muc=muc_of,
                       of_muc_relerr=abs(muc_of - muc_ana) / muc_ana)
            flow["of"].append(muc_of)  # proxy; flow-rate ratio handled below
        else:
            row.update(of_l2=None, of_linf=None, of_muc=None, of_muc_relerr=None)
            flow["of"].append(None)

        # --- MUG error ---
        if mug is not None:
            y_m, ux_m, mug_file = mug
            u_ana_on_m = ah.u_profile(y_m, A, ha)
            ana_n_m, _ = norm_by_centerline(y_m, u_ana_on_m)
            m_n, _ = norm_by_centerline(y_m, ux_m)
            diffm = m_n - ana_n_m
            l2m = np.sqrt(trapz(diffm**2, y_m) / trapz(ana_n_m**2, y_m))
            linfm = np.max(np.abs(diffm))
            muc_m = mean_over_uc(y_m, ux_m)
            row.update(mug_l2=l2m, mug_linf=linfm, mug_muc=muc_m,
                       mug_muc_relerr=abs(muc_m - muc_ana) / muc_ana)
            flow["mug"].append(muc_m)
        else:
            row.update(mug_l2=None, mug_linf=None, mug_muc=None, mug_muc_relerr=None)
            flow["mug"].append(None)

        rows.append(row)

        # --- per-Ha CSV ---
        csv = os.path.join(RESULTS, f"profile_Ha{ha}.csv")
        with open(csv, "w") as fh:
            fh.write("y,u_analytic_norm,u_mhdfoam_norm,u_mug_norm\n")
            yg = y_of if y_of is not None else y_plot
            ana_g = ah.u_profile(yg, A, ha)
            ana_gn, _ = norm_by_centerline(yg, ana_g)
            of_gn = norm_by_centerline(yg, ux_of)[0] if of is not None else [None] * len(yg)
            if mug is not None:
                mug_gn = np.interp(yg, mug[0], norm_by_centerline(mug[0], mug[1])[0])
            else:
                mug_gn = [None] * len(yg)
            for i in range(len(yg)):
                a_ = f"{ana_gn[i]:.6e}"
                o_ = f"{of_gn[i]:.6e}" if of is not None else ""
                m_ = f"{mug_gn[i]:.6e}" if mug is not None else ""
                fh.write(f"{yg[i]:.6e},{a_},{o_},{m_}\n")

        # --- overlay plot ---
        if HAVE_MPL:
            plt.figure(figsize=(6, 4.2))
            plt.plot(u_ana_plot_n, y_plot, "k-", lw=2, label="analytic")
            if of is not None:
                of_pn, _ = norm_by_centerline(y_of, ux_of)
                plt.plot(of_pn, y_of, "o", ms=3, mfc="none", color="C0",
                         label="mhdFoam")
            if mug is not None:
                m_pn, _ = norm_by_centerline(mug[0], mug[1])
                plt.plot(m_pn, mug[0], "s", ms=3, color="C3", label="MUG")
            plt.xlabel("u / u_centerline")
            plt.ylabel("y / a")
            plt.title(f"Hartmann profile, Ha = {ha}")
            plt.legend()
            plt.grid(alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(RESULTS, f"overlay_Ha{ha}.png"), dpi=130)
            plt.close()

    # --- flow-rate vs Ha plot ---
    if HAVE_MPL:
        Ha = np.array(flow["Ha"], float)
        plt.figure(figsize=(6, 4.2))
        qstar = np.array([ah.flow_rate_star(h) for h in Ha])
        plt.loglog(Ha, qstar, "k-o", label="analytic Q*")
        plt.loglog(Ha, 1.0 / Ha, "k--", alpha=0.5, label="1/Ha reference")
        plt.xlabel("Ha")
        plt.ylabel("Q*  = (1/Ha)[1 - tanh(Ha)/Ha]")
        plt.title("Hartmann flow-rate suppression")
        plt.legend()
        plt.grid(alpha=0.3, which="both")
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS, "flowrate_vs_Ha.png"), dpi=130)
        plt.close()

    # --- error table ---
    with open(os.path.join(RESULTS, "error_table.md"), "w") as fh:
        fh.write("# Hartmann three-way verification — error table\n\n")
        fh.write("Profiles normalized by centerline velocity. "
                 "L2/Linf are on the normalized profile vs analytic; "
                 "u_mean/u_c is the discriminating scalar (2/3 at Ha=0, ->1 plug core).\n\n")
        fh.write("| Ha | u_mean/u_c (analytic) | mhdFoam L2 | mhdFoam Linf | mhdFoam u_mean/u_c | mhdFoam relerr | MUG L2 | MUG Linf | MUG u_mean/u_c | MUG relerr |\n")
        fh.write("|----|----|----|----|----|----|----|----|----|----|\n")
        def fmt(x, e=False):
            if x is None:
                return "—"
            return f"{x:.3e}" if e else f"{x:.4f}"
        for r in rows:
            fh.write("| {Ha} | {a:.4f} | {ol2} | {oli} | {omuc} | {ore} | {ml2} | {mli} | {mmuc} | {mre} |\n".format(
                Ha=r["Ha"], a=r["muc_ana"],
                ol2=fmt(r["of_l2"], True), oli=fmt(r["of_linf"], True),
                omuc=fmt(r["of_muc"]), ore=fmt(r["of_muc_relerr"], True),
                ml2=fmt(r["mug_l2"], True), mli=fmt(r["mug_linf"], True),
                mmuc=fmt(r["mug_muc"]), mre=fmt(r["mug_muc_relerr"], True)))
        fh.write("\n_Generated by cases/hartmann/compare.py_\n")

    print("Wrote results to", RESULTS)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
