"""
V2 transient-conjugate waveform extractor — turns a chain of MUG restart
checkpoints (xmhd2d_NNNNN.rst, written by test_hunt_mug at rst_freq) into the
spatiotemporal fields v(t,y) and by(t,y) for the V2 capstone comparison against
COMSOL's Tier6 transient conjugate reference.

Each .rst is an HDF5 checkpoint holding U_vely, U_by (+ all fields), plus scalars
t and dt. This reads a directory of them in time order and extracts, along the
Hartmann centerline cut (x=0), the axial velocity vely(t,z) and induced field
by(t,z) — the two waveforms whose spatiotemporal L2 + Alfven-transit timing is
the V2 metric.

Requires the mesh node coordinates (same order as the vector DOFs). We read them
from the companion profile the driver writes, or pass --coords a "x z ..." file.

Usage: python3 v2_waveform.py --rstdir <case>/ --coords <case>/hunt_mug.profile
       [--xcut 0.0] [--out v2_waveform.npz]
"""
import argparse
import glob
import os

import numpy as np

try:
    import h5py
except ImportError:
    h5py = None


def load_rst(path):
    f = h5py.File(path, "r")
    t = float(np.array(f["t"]).ravel()[0])
    vely = np.array(f["U_vely"]).ravel()
    by = np.array(f["U_by"]).ravel()
    f.close()
    return t, vely, by


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rstdir", required=True)
    ap.add_argument("--coords", required=True,
                    help="profile file whose rows are 'x z ...' in DOF order")
    ap.add_argument("--xcut", type=float, default=0.0)
    ap.add_argument("--out", default="v2_waveform.npz")
    args = ap.parse_args()
    if h5py is None:
        raise SystemExit("h5py required: pip install h5py")

    coords = np.loadtxt(args.coords)
    xs, zs = coords[:, 0], coords[:, 1]
    # nodes on the Hartmann centerline cut x = xcut
    cut = np.where(np.abs(xs - args.xcut) < 1e-6)[0]
    if cut.size == 0:
        # fall back to nearest x-column
        xu = np.unique(np.round(xs, 6))
        xnear = xu[np.argmin(np.abs(xu - args.xcut))]
        cut = np.where(np.abs(xs - xnear) < 1e-6)[0]
        print(f"[note] no nodes exactly at x={args.xcut}; using x={xnear}")
    order = np.argsort(zs[cut])
    cut = cut[order]
    zcut = zs[cut]

    files = sorted(glob.glob(os.path.join(args.rstdir, "xmhd2d_*.rst")))
    rows = []
    for fp in files:
        t, vely, by = load_rst(fp)
        rows.append((t, vely[cut], by[cut]))
    rows.sort(key=lambda r: r[0])
    tarr = np.array([r[0] for r in rows])
    vgrid = np.array([r[1] for r in rows])   # shape (nt, nz)
    bgrid = np.array([r[2] for r in rows])

    np.savez(args.out, t=tarr, z=zcut, vely=vgrid, by=bgrid)
    print(f"V2 waveform: {tarr.size} times t=[{tarr[0]:.4g},{tarr[-1]:.4g}], "
          f"{zcut.size} points on x={args.xcut} cut")
    print(f"  vely range [{vgrid.min():.3e},{vgrid.max():.3e}], "
          f"by range [{bgrid.min():.3e},{bgrid.max():.3e}]")
    print(f"  wrote {args.out} (t, z, vely[t,z], by[t,z]) — hand to the COMSOL overlay")


if __name__ == "__main__":
    main()
