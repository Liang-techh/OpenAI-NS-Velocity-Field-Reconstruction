# MATLAB visualization

From the repository root:

```matlab
addpath('visualization/matlab');
ns_explorer;       % full-field explorer
ns_core_explorer;  % dense central-core view
```

ST054-Q2 is the lower-L2 default; ST054-M3 is the lower-peak alternative. data/st054_models.mat contains both complete finite-basis models. Neither meets the full NS momentum gates.

## Full-field explorer

The bottom slider evaluates continuous time over [0.25,0.75] from the original eight-term Chebyshev expansion. Select instantaneous 3-D streamlines, vorticity surfaces or velocity arrows, with spectral slices of velocity, vorticity, pressure, pressure force and full residual. Seed placement, line count/length, opacity, cutaway and color limits change the view, not the field. Coordinates are not stretched to resemble an illustration. Radial view crops the axes; use the core app for local resampling.

A streamline is not a particle trajectory in a time-dependent field. Only time changes the evaluated field; other controls never change coefficients, forcing or viscosity. Keep camera, scales and seeds fixed for controlled comparisons. Use zero slice offset to see the axis.

## Central-core explorer

The core app re-evaluates a local x/y box and independent z extent directly from the spectral field. It does not enlarge a few coarse cells. Core radius, half-height, seed radius, seed height span, count, arc-length and slice offset are independent controls. Dragging gives a lower-detail preview; release rebuilds detail. Time and Play are at the bottom. Seeds are clipped inside the sampled box.

Metrics use actual enstrophy with cylindrical volume weights, not unweighted samples or substituted speed. Extents depend on the selected observation cylinder. The excluded midplane strip and inflow cutoff are explicit; undefined near-zero-inflow ratios are excluded. These are geometry diagnostics, not source-calibrated optimization constraints or proof of a matching mechanism.

```matlab
T = ns_core_timeseries([], 'ST054-Q2', 'Times', linspace(0.25,0.75,13));
```

This samples a fixed physical window over time. It does not fit shrinking scales, infer singularity or measure target-image similarity. The core tools replace the earlier untested supplement with bounded tensor sampling, valid MATLAB indexing and explicit volume weighting.

The full app can import [time,x,y,z] velocity MAT grids. They use interpolation and approximate curl; missing pressure/residual is not invented. Focused core tools require the spectral bundle. Do not extrapolate beyond the declared time window. Performance depends on hardware and detail.
