"""
Pytest driver for test_drag_decay: validates the reduced Hartmann/wall-drag
source term added to xmhd_2d.F90.

The Fortran program writes drag_decay.results with two lines:
  Line 1: relative error of velx vs backward-Euler prediction
  Line 2: relative error of vely (should be undamped)

Pass criteria:
  - velx relative error < 2%   (backward-Euler prediction)
  - vely relative error < 1e-6  (drag_bhat direction, should not decay)
"""
from __future__ import print_function
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

# Input file template — periodic domain, 2D surface mesh
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
 use_petsc=F
/

&mesh_options
 meshname='test'
 cad_type=92
 nlevels={1}
 nbase={0}
 grid_order=2
 fix_boundary=T
/

&cube_options
 mesh_type=1
 ni=4,1,0
 rscale=1.0,1.0,0.
 shift=-0.5,-0.5,0.
 ref_per=T,T,F
/

&drag_test_options
  order={2}
  nsteps={3}
  dt={4}
  n0=1.E19
  velx0=1.0
  vely0=0.75
  t0=1.E3
  chi=1.E-30
  eta=1.E-30
  nu=1.E-30
  D_diff=1.E-30
  gamma=1.67
  den_scale=1.E19
  alpha_drag={5}
  lin_tol=1.E-10
  nl_tol=1.E-8
  pm=F
/
"""


def drag_decay_setup(nbase, nlevels, order=2, nsteps=10, dt=0.05, alpha_drag=2.0):
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(
            nbase, nlevels, order, nsteps, dt, alpha_drag))
    return run_OFT("./test_drag_decay", 1, 120)


def validate_drag_result(velx_err_tol=0.02, vely_err_tol=1.e-6):
    """Read drag_decay.results and check tolerances."""
    retval = True
    with open('drag_decay.results', 'r') as fid:
        velx_err = float(fid.readline())
        vely_err = float(fid.readline())

    print("  velx rel. error vs backward-Euler: {:.4e}  (tol={:.4e})".format(
        velx_err, velx_err_tol))
    print("  vely rel. error (should be ~0):     {:.4e}  (tol={:.4e})".format(
        vely_err, vely_err_tol))

    if velx_err > velx_err_tol:
        print("FAILED: velx decay error too large!")
        retval = False
    if vely_err > vely_err_tol:
        print("FAILED: vely changed (should be undamped by drag)!")
        retval = False
    return retval


# ---------------------------------------------------------------------------
# Test: basic drag decay — order 2, single level
# ---------------------------------------------------------------------------
def test_drag_p2_r1():
    """Drag decay test: p=2, 10 steps, alpha=2, dt=0.05"""
    assert drag_decay_setup(1, 1, order=2, nsteps=10, dt=0.05, alpha_drag=2.0)
    assert validate_drag_result(velx_err_tol=0.02, vely_err_tol=1.e-6)


# ---------------------------------------------------------------------------
# Test: higher drag coefficient (alpha=5) — faster decay, same structure
# ---------------------------------------------------------------------------
def test_drag_p2_r1_high_alpha():
    """Drag decay test: p=2, 10 steps, alpha=5, dt=0.05"""
    assert drag_decay_setup(1, 1, order=2, nsteps=10, dt=0.05, alpha_drag=5.0)
    assert validate_drag_result(velx_err_tol=0.02, vely_err_tol=1.e-6)


# ---------------------------------------------------------------------------
# Test: higher polynomial order
# ---------------------------------------------------------------------------
def test_drag_p3_r1():
    """Drag decay test: p=3, 10 steps, alpha=2, dt=0.05"""
    assert drag_decay_setup(1, 1, order=3, nsteps=10, dt=0.05, alpha_drag=2.0)
    assert validate_drag_result(velx_err_tol=0.02, vely_err_tol=1.e-6)
