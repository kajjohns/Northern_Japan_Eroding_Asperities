# NoJapan Eroding Asperities Model Result Tables

This directory contains tab-delimited text files derived from the Northern Japan eroding-asperities MCMC model. The files are intended for archiving with the code repository `NoJapan_eroding_asperities` and for use in reproducing the model-result panels plotted by `plot_fancy_MCMC_inversion_joint.m`.

The files were generated with:

```matlab
setup_NorthernJapan
export_model_results_for_zenodo
```

## Files

Patch-based model-result files:

```text
coupling_ratio_1998.txt
coupling_ratio_2009.txt
locking_probability_1998.txt
locking_probability_2009.txt
process_zone_stressing_rate_1998.txt
process_zone_stressing_rate_2009.txt
```

Velocity-fit files:

```text
velocity_fit_1998.txt
velocity_fit_2009.txt
```

## Patch-Based Files

Each row describes one triangular fault patch. These files correspond to quantities plotted with `trisurf` in `plot_fancy_MCMC_inversion_joint.m`.

Columns:

```text
patch_id
centroid_lon
centroid_lat
node1_id
node1_lon
node1_lat
node2_id
node2_lon
node2_lat
node3_id
node3_lon
node3_lat
<quantity>_mean
<quantity>_std
```

Column notes:

- `patch_id`: row number of the triangular patch in the model mesh.
- `centroid_lon`, `centroid_lat`: longitude and latitude of the patch centroid.
- `node*_id`: index of the mesh node used by the patch.
- `node*_lon`, `node*_lat`: longitude and latitude of each patch node.
- `<quantity>_mean`: posterior mean of the exported model quantity.
- `<quantity>_std`: posterior standard deviation of the exported model quantity.

Patch quantities:

- `coupling_ratio`: dimensionless coupling ratio, computed as `1 - creep_rate / plate_rate`.
- `locking_probability`: posterior probability that a patch is locked, from the mean of the binary locked-patch indicator.
- `process_zone_stressing_rate`: reconstructed process-zone stressing-rate quantity plotted by the MATLAB scripts, using the same scaling as `plot_fancy_MCMC_inversion_joint.m`.

## Velocity-Fit Files

Each row describes one observation site.

Columns:

```text
site_id
lon
lat
observed_east_mm_per_yr
observed_north_mm_per_yr
modeled_east_mean_mm_per_yr
modeled_north_mean_mm_per_yr
modeled_east_std_mm_per_yr
modeled_north_std_mm_per_yr
```

Column notes:

- `site_id`: row number of the observation site.
- `lon`, `lat`: longitude and latitude of the observation site.
- `observed_east_mm_per_yr`, `observed_north_mm_per_yr`: observed horizontal velocity components.
- `modeled_east_mean_mm_per_yr`, `modeled_north_mean_mm_per_yr`: posterior mean modeled horizontal velocity components.
- `modeled_east_std_mm_per_yr`, `modeled_north_std_mm_per_yr`: posterior standard deviation of modeled horizontal velocity components.

## Coordinates and Units

- Geographic coordinates are longitude and latitude in decimal degrees.
- Velocity components are in millimeters per year.
- Coupling ratio and locking probability are dimensionless.
- The patch values are associated with triangular fault patches; use the centroid columns for point-style plotting or the node columns to reconstruct patch geometry.

## Provenance

These tables are derived from the MCMC output files in `NoJapan_outputs_M9cycle_joint/` using `export_model_results_for_zenodo.m`. The export script loads the MCMC chain outputs, removes the configured number of burn-in samples, computes posterior means and standard deviations, and writes the result tables in this directory.
