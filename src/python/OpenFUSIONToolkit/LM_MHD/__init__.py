#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''Utilities for liquid-metal MHD regime triage and coupling data exchange.'''

from .coupling import (
    MU0,
    CouplingSnapshot,
    LiquidMetalProperties,
    alfven_mach_number,
    alfven_speed,
    blanket_feedback_ratio_sheet,
    classify_coupling,
    feedback_ratio,
    flow_time,
    hartmann_number,
    interaction_parameter,
    magnetic_diffusion_time,
    magnetic_reynolds_number,
    motional_current_density,
    read_coupling_snapshot,
    representative_liquid_metals,
    reynolds_number,
    sheet_delta_b,
    write_coupling_snapshot,
)
from .closures import (
    cells_per_hartmann_layer,
    classify_hartmann_resolution,
    drag_power_density,
    hartmann_channel_effective_drag,
    hartmann_layer_thickness,
    hartmann_mean_over_centerline,
    implicit_drag_update,
    parallel_velocity,
    perpendicular_velocity,
)

__all__ = [
    "MU0",
    "CouplingSnapshot",
    "LiquidMetalProperties",
    "alfven_mach_number",
    "alfven_speed",
    "blanket_feedback_ratio_sheet",
    "classify_coupling",
    "feedback_ratio",
    "flow_time",
    "hartmann_number",
    "interaction_parameter",
    "magnetic_diffusion_time",
    "magnetic_reynolds_number",
    "motional_current_density",
    "read_coupling_snapshot",
    "representative_liquid_metals",
    "reynolds_number",
    "sheet_delta_b",
    "write_coupling_snapshot",
    "cells_per_hartmann_layer",
    "classify_hartmann_resolution",
    "drag_power_density",
    "hartmann_channel_effective_drag",
    "hartmann_layer_thickness",
    "hartmann_mean_over_centerline",
    "implicit_drag_update",
    "parallel_velocity",
    "perpendicular_velocity",
]
