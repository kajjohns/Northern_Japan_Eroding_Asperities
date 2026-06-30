[![DOI](https://zenodo.org/badge/1284458728.svg)](https://doi.org/10.5281/zenodo.21073351)

# NoJapan_eroding_asperities

MATLAB code for conducting inversions for manuscript submitted to GRL:  Decadal Creep Acceleration on the Northern Japan Trench Driven by Erosion of Heterogeneously Locked Asperities, 
Kaj M. Johnson and Eric D. Burton, Department of Earth and Atmospheric Sciences, Indiana University. 

The code models eroding asperities and interseismic coupling along the northern Japan subduction margin. The main workflow estimates asperity radius changes and process-zone stressing using two interseismic velocity fields, Sanriku postseismic corrections, elastic Green's functions, and a viscoelastic M9-cycle contribution.

The intended run order is:

```matlab
setup_NorthernJapan
mcmc_inversion_M9cycle_joint
plot_fancy_MCMC_inversion_joint
```

`setup_NorthernJapan.m` loads and prepares the shared model inputs, `mcmc_inversion_M9cycle_joint.m` runs the MCMC inversion, and `plot_fancy_MCMC_inversion_joint.m` loads the output files and makes summary figures.

## Repository Layout

```text
.
├── setup_NorthernJapan.m
├── mcmc_inversion_M9cycle_joint.m
├── plot_fancy_MCMC_inversion_joint.m
├── data/
├── tools/
├── hmmvp0.16/
├── NoJapan_outputs_M9cycle_joint/
└── Hmat_ne15326hmat1rerr-3.00.dat
```

Key files and folders:

- `setup_NorthernJapan.m`: loads mesh, Green's functions, velocity data, earthquake catalog, asperity geometry, rupture scaling, and Sanriku postseismic correction fields.
- `mcmc_inversion_M9cycle_joint.m`: runs the MCMC sampler for asperity radius scales and process-zone ring stresses.
- `plot_fancy_MCMC_inversion_joint.m`: reconstructs model fields from MCMC output and makes publication-style figures.
- `data/`: required model input data.
- `tools/`: helper functions for asperity geometry, coordinate conversion, plotting, DEM handling, and forward modeling.
- `hmmvp0.16/`: H-matrix matrix-vector product code used by the mechanical forward model.
- `NoJapan_outputs_M9cycle_joint/`: MCMC output folder. Existing output can be used for plotting or continuing a chain.

## Requirements

This package is written for MATLAB.

Expected MATLAB functionality includes:

- Basic MATLAB numerical and plotting functions.
- `gmres` for solving the creep-rate system.
- `griddata`, `trisurf`, `tiledlayout`, and standard graphics functions.
- Mapping/geospatial support for `geotiffread` and `shaperead`, used by the DEM and coastline plotting helpers.
- A working `hm_mvp` MEX binary compatible with your platform, or the ability to rebuild it from `hmmvp0.16/`.

The repository currently includes precompiled `hm_mvp` binaries for macOS architectures. If MATLAB cannot run `hm_mvp`, rebuild the MEX files from `hmmvp0.16/` for your system.

## Data Inputs

The active workflow reads these files from `data/`:

```text
data/setup_mesh_minimal.mat
data/interseismic_velocities_1998_2009.txt
data/earthquake_catalog.txt
data/Sanriku_sep19_2022.txt
data/M9_points.mat
data/visco_cycle_M9_contribution.mat
```

The topographic and coastline plotting utilities use:

```text
tools/dem/ETOPO1_Ice_g_geotiff.tif
tools/ne_10m_coastline/ne_10m_coastline.*
```

## Workflow

### 1. Start MATLAB in the repository root

Set MATLAB's current folder to the repository root, the folder containing the three top-level scripts.

```matlab
cd path/to/NoJapan_eroding_asperities
```

### 2. Run setup

```matlab
setup_NorthernJapan
```

This script adds required paths and creates the workspace variables used by the MCMC and plotting scripts, including:

- mesh and observation-site variables,
- horizontal Green's function matrices,
- observed 1998 and 2009 velocity fields,
- asperity geometry and rupture scaling,
- Sanriku postseismic correction fields.

### 3. Run the MCMC inversion

```matlab
mcmc_inversion_M9cycle_joint
```

The main user-editable settings are near the top of the script:

- `continuing`: if `true`, resume from the last row of the existing output files.
- `folder_name`: output folder for the chain.
- `startval` and `starting_folder`: optionally initialize from a different output folder.
- `radii`, `stepsize_radii`: initial radius-scale values and proposal step sizes.
- `ring_taus`, `stepsize_ring_tau`: initial process-zone stress values and proposal step sizes.
- `D`: process-zone ring width in km.
- `mu`: elastic shear modulus in Pa.

The script writes chain state to text files in `folder_name`.

### 4. Plot results

After running `setup_NorthernJapan.m` and after MCMC outputs exist, run:

```matlab
plot_fancy_MCMC_inversion_joint
```

The plotting script expects the setup variables to still be available in the MATLAB workspace. It loads the MCMC output files, reconstructs ring stressing rates, and makes figures for coupling, locking probability, process-zone stressing, velocity fits, and moment-rate summaries.

### 5. Export model result tables for Zenodo

After running `setup_NorthernJapan.m` and after MCMC outputs exist, run:

```matlab
export_model_results_for_zenodo
```

This writes tab-delimited text files to `Zenodo_model_results/`. The patch-based files include one row per triangular mesh patch, with patch centroid lon/lat, the lon/lat of the three mesh nodes, posterior mean, and posterior standard deviation. The velocity files include observation-site lon/lat plus observed and posterior-mean modeled east/north velocity components for 1998 and 2009.

## MCMC Outputs

The MCMC script writes these files to `folder_name`:

```text
M_radii.txt
M_ring_tau.txt
logrho.txt
dhat.txt
locked_index.txt
creep_rates.txt
visco_scale.txt
```

These are plain text chain outputs used by `plot_fancy_MCMC_inversion_joint.m`.

## Notes for Users

- Run scripts from the repository root so relative paths resolve correctly.
- The setup script must be run before the MCMC script.
- The setup script must also be run before the plotting script unless the required variables are already in the workspace.
- Long MCMC runs append to output files when `continuing = true`.
- Set `continuing = false` only when you want to overwrite the output files in `folder_name`.
- If `hm_mvp` fails to initialize, check that the MEX binary matches your MATLAB platform.

## Helper Code

Important helper files in `tools/` include:

- `calculate_M9cycle_joint_forward_model.m`: evaluates one MCMC proposal.
- `define_asperities_seismic_M9.m`: defines polygon and circular asperity geometry.
- `make_rupture_dimension_M9.m`: computes rupture areas and scale factors.
- `get_locked_indices_asperities.m`: identifies locked mesh patches for a proposal.
- `get_ringtau_asperities.m`: builds process-zone ring stress fields.
- `mvp.m`: wrapper around `hm_mvp`.
- `plot_asperities_llh.m`, `plot_coast.m`: plotting helpers.

## Citation and License

If you use this package in a publication, cite the associated study and any external datasets or utilities used by your analysis. The DEM helper includes its own license file in `tools/dem/license.txt`; the coastline data include Natural Earth metadata in `tools/ne_10m_coastline/`.
