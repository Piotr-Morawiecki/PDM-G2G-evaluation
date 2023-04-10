# Evaluation of PDM and G2G rainfall-runoff models
Scripts used to evaluate PDM model (R. Lamb, 1999) and Grid-to-Grid model (V. Bell, 2007)

## About the repository

This GitHub repository includes scripts, which were used in two-part paper "Asymptotic analysis and validation of conceptual rainfall-runoff models" by P. Morawiecki and P.H. Trinh [1,2] for comparing Probability Distributed Model (PDM) implementation by R. Lamb [3] and Grid-to-Grid (G2G) model [4] to the 1D physical benchmark model developed by us in [5].

The repository includes MODELS/ directory with MatLAB implementations of PDM, G2G and 1D physical benchmark model, all applied for a simple 1D benchmark geometry representing a single hillslope. The scripts in the main directory recreate figures from papers [1,2]. The visual appearance of the figures is different than in paper, since for the published version the output data were illustrated using LaTeX TikZ library.

## Repository content

The repository includes the following directories:

* DATA - directory for output data generated by the scripts
* FIGURES - directory for storing output figures
* MODELS - directory with implementations of PDM, G2G and 1D physical benchmark model

The main directory includes the following scripts:

| Filename | Short description | Figures recreated |
| ---      | ---       | ---       |
| MODELS/PDM.m | Class implementing the PDM implementation from [3]       | - |
| MODELS/G2G.m | Class implementing the Grid-to-Grid (G2G) model from [4] | - |
| MODELS/Catchment1D.m | Class implementing the 1D physical benchmark model from [5] | - |
| MODELS/MvG_model.m | Class implementing the Mualem-Van Genuchten model used in the 'Catchment1D' class | - |
| Physical_benchmarks.m | Generates physical benchmark hydrographs with varying P and P0 | - |
| Physical_drying_scenario.m | Generates physical benchmark hydrograph for the drying scenario | Fig. 9 in [2] |
| PDM_regimes.m | Illustrates three regimesof the PDM | Fig. 5 in [1] |
| PDM_timescales.m | Illustrates three distinct timescales characterising the PDM | Fig. 6 in [1] |
| PDM_P0_dependance.m | Compares PDM results with physical benchmarks for varying mean rainfall P0 | Fig. 8, 10 in [1] |
| PDM_P_dependance.m | Compares PDM results with physical benchmarks for varying simulated rainfall P  | Fig. 9, 11 in [1] |
| G2G_P_dependance.m | Compares G2G results with physical benchmarks for varying simulated rainfall P | Fig. 3 in [2] |
| G2G_P0_depandance.m | Compares G2G results with physical benchmarks for varying mean rainfall P0 | Fig. 4 in [2] |
| G2G_drying_scenario.m | Compares G2G results with physical benchmarks for the drying scenario | Fig. 5, 6 in [2] |
| G2G_upstream_rainfall.m | Illustrates differences between the surface flow model used G2G and physical models when applied to an upstream rainfall scenario | Fig. 7, 8 in [2] |

## References

<a id="1">[1]</a> P. W. Morawiecki and P. H. Trinh. "Asymptotic analysis and validation of conceptual rainfall-runoff models. Part 1. Evaluation of a lumped Probability-Distributed Model." T.B.C. (2023).

<a id="1">[1]</a> P. W. Morawiecki and P. H. Trinh. "Asymptotic analysis and validation of conceptual rainfall-runoff models. Part 2. Evaluation of a Grid and Grid-to-Grid models." T.B.C. (2023).

<a id="1">[3]</a> R. Lamb. "Calibration of a conceptual rainfall‐runoff model for flood frequency estimation by continuous simulation." Water Resources Research 35.10 (1999): 3103-3114.

<a id="1">[4]</a> V. A. Bell et al. "Development of a high resolution grid-based river flow model for use with regional climate model output." Hydrology and Earth System Sciences 11.1 (2007): 532-549.

<a id="1">[5]</a> P. W. Morawiecki and P. H. Trinh. "On the development and analysis of coupled surface-subsurface models of catchments. Part 3. Analytical solutions and scaling laws." T.B.C. (2022).
