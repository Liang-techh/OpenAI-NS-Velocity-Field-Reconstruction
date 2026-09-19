# MATLAB visualization

## Launch

From the repository root in MATLAB:

```matlab
start_here
start_here('core')
```

The default presentation preset opens `ns_explorer` with the selected seed radius `0.621`, seed-height center `-0.221`, radial view `0.938`, color multiplier `0.932` and time `0.5`. The slice is intentionally placed on the axis (`y=0`), rather than at the off-axis value in the earlier screenshot. These are viewing choices, not modifications to the field.

The direct entry points remain available:

```matlab
addpath('visualization/matlab')
ns_explorer
ns_core_explorer
```

Normal operation needs base MATLAB and the included MAT data, not Python. The original viewer was natively tested on MATLAB R2024b. Release-specific tests, including corrected core tools, are recorded separately by the release workflow; do not infer new GUI testing solely from old receipts. R2021a+ is a compatibility target, not an exhaustive version certification.

## Full-field controls

The bottom physical-time slider and Play button evaluate the continuous finite-basis field in `[0.25,0.75]`. Other controls select the candidate, instantaneous streamlines/vorticity surface/arrows, scalar and line coloring, XZ/XY/YZ slices and their offset, seed radius and height, seed count, arc cap, isovalue, opacity, radial display range, cutaway and fixed color multiplier.

Color and arrow gains are explicit. A full-field radial zoom changes the displayed limits, not the underlying sampling resolution. Drag/play uses a preview grid; releasing the slider restores the selected detail level. Full-field seed count is the number of starting seeds, not a guarantee that every seed produces a long line.

## Corrected central-core tools

`ns_core_explorer` evaluates a NEW local box when radius, half-height, time or candidate changes. It uses the separable evaluator rather than building huge per-point basis tensors. Seed radius and vertical span are independent and clamped inside the box. The time slider is below the plots, supports preview and playback, and all controls retain the frozen coefficients.

The release repairs the earlier supplement's invalid temporary-result indexing, mistaken struct-method call, UI-layout construction and streamline vertex-budget calculation. It also removes redundant full-volume evaluation during every metric update. Imported sampled grids stay in the full viewer; the core viewer refuses them rather than inventing pressure/residual arrays or presenting interpolated data as spectral evaluation.

The geometry table uses a cylinder inside the local square box, Cartesian trapezoid volume weights and actual vorticity squared. It reports radial/axial RMS extents, their ratio, flow/pressure directions, component ratios and forward instantaneous-streamline summaries. Axis points are excluded from radial-sign/ratio tests; bipolar/axial-pressure directions exclude `|z|<=0.03`; component ratios exclude radial inflow weaker than `1e-5`. Missing/undefined measurements are not silently turned into successes.

These are local, observation-window-dependent statistics. They are NOT a detected physical core boundary or a numerical similarity score against the original schematic. Finite streamline turns depend on seeding, arc cap and boundary termination, and are NOT time-integrated particle rotations.

## Time-series diagnostics

```matlab
T = ns_core_timeseries([], 'ST054-Q2', ...
    'Times', linspace(0.25,0.75,13), ...
    'RadiusMax', 0.35, 'ZHalfSpan', 0.60);
```

The window is fixed in physical coordinates. The output table contains the time and local metrics. Separate figures show extents, anisotropy, direction fractions and component ratios. Do not interpret window-truncated extent changes as a fitted singularity exponent.

## Recommended comparison discipline

Keep the camera, physical axis ratio, seeds, color scale and observation window fixed when comparing times or candidates. Move the XZ slice to `y=0` when inspecting the axis. Use a narrower seed-height span for the central core; use the full viewer for return flow and support-edge residuals. Never stretch one axis or alter a coefficient merely to mimic a schematic.

The MAT payload is identical to the pinned natively tested upstream export. Its velocity, pressure and residual reference arrays are checked against the Python package. Evaluator agreement near machine precision must not be confused with the much larger NS residual.
