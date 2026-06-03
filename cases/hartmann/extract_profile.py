"""
Extract a wall-normal velocity profile u_x(y) from an OpenFOAM case by reading the
foamToVTK output with meshio. This sidesteps the broken sampling functionObject
("sha1" IOstream error) in this OpenFOAM v1912 packaging.

Usage:
    python extract_profile.py <caseDir> [xtarget=10.0] [time=latest]
Writes <caseDir>/profile.csv with columns: y,ux  and prints summary.

Importable: get_profile(caseDir, xtarget=10.0) -> (y, ux) sorted by y.
"""
import os
import sys
import glob
import subprocess
import numpy as np
import meshio

OF_BASHRC = "/usr/share/openfoam/etc/bashrc"


def _ensure_vtk(case, time=None):
    """Run foamToVTK for the requested time if no VTK output is present."""
    vtk_dir = os.path.join(case, "VTK")
    have = glob.glob(os.path.join(vtk_dir, "*", "internal.vtu"))
    if not have:
        tflag = f"-time {time}" if time else "-latestTime"
        cmd = f"source {OF_BASHRC} >/dev/null 2>&1 && cd '{case}' && foamToVTK {tflag} >foamToVTK.log 2>&1"
        subprocess.run(["bash", "-lc", cmd], check=False)
        have = glob.glob(os.path.join(vtk_dir, "*", "internal.vtu"))
    return have


def get_profile(case, xtarget=10.0, time=None):
    """Return (y, ux) at the cell-column nearest x=xtarget, sorted by y."""
    files = _ensure_vtk(case, time)
    if not files:
        return None
    # choose the VTK dir with the largest trailing index (latest written)
    def idx(p):
        b = os.path.basename(os.path.dirname(p))
        try:
            return int(b.split("_")[-1])
        except ValueError:
            return -1
    f = max(files, key=idx)
    m = meshio.read(f)
    if "U" not in m.cell_data:
        return None
    U = np.vstack(m.cell_data["U"])  # (ncells, 3)
    conn = m.cells[0].data           # (ncells, 8)
    cent = m.points[conn].mean(axis=1)  # cell centroids (ncells,3)
    xc, yc = cent[:, 0], cent[:, 1]
    # pick the single column of cells nearest x=xtarget
    xs = np.unique(np.round(xc, 6))
    xcol = xs[np.argmin(np.abs(xs - xtarget))]
    sel = np.abs(xc - xcol) < 1e-6
    y = yc[sel]
    ux = U[sel, 0]
    order = np.argsort(y)
    return y[order], ux[order]


def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    case = sys.argv[1]
    xtarget = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
    time = sys.argv[3] if len(sys.argv) > 3 else None
    res = get_profile(case, xtarget, time)
    if res is None:
        print("FAILED to extract profile from", case); sys.exit(2)
    y, ux = res
    out = os.path.join(case, "profile.csv")
    np.savetxt(out, np.column_stack([y, ux]), delimiter=",", header="y,ux", comments="")
    print(f"wrote {out}: {len(y)} points, x~{xtarget}, "
          f"u_c={np.interp(0,y,ux):.4f}, y in [{y.min():.3f},{y.max():.3f}]")


if __name__ == "__main__":
    main()
