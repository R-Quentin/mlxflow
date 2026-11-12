# mlxflow

Model building and diagnostics workflow for population PK models run with
[Monolix](https://monolixsuite.slp-software.com/).

> **Unofficial.** `mlxflow` is not affiliated with, nor endorsed by, Lixoft or
> Simulations Plus. "Monolix" is their trademark; the `mlx` prefix only says that
> the package works on Monolix files.

Work in progress: version `0.0.1`.

## What it does

| Layer | Needs Monolix running? | Main functions |
|---|---|---|
| **Engine** | yes (through `lixoftConnectors`) | `backward_elimination()`, `run_cossac()`, `init_monolix()` |
| **Reading** | no, only files | `read_summary()`, `read_models()`, `read_model_building_report()`, `compare_two_models()`, `compare_models()` |
| **Reporting** | no, only files | `run_diagnostics()`, `run_comparison()`, `build_param_table()`, `build_*_data()`, `plot_flowchart()`, `plot_flowchart_combined()` |

Only the engine layer needs `lixoftConnectors`, which ships with MonolixSuite (it is not on CRAN):
it is listed in `Suggests` and checked when a function really talks to Monolix.

## Installation

```r
# install.packages("remotes")
remotes::install_github("R-Quentin/mlxflow")
```

`ggplot2 >= 3.5.0` is required.

## Quick tour

```r
library(mlxflow)

# Backward elimination of the covariates of a model (LRT, alpha = 0.01)
init_monolix("warfarin_final.mlxtran")
res <- backward_elimination("warfarin_final.mlxtran", alpha = 0.01,
                            output_dir = "warfarin_final/Backward")

# Every model of the run, then the flowchart of the search
models <- read_models("warfarin_final/Backward") |>
  dplyr::mutate(covariates = purrr::map(covariates, clean_cov_labels))
plot_flowchart(build_backward_data(models, alpha = 0.01), title = "Backward elimination")

# Diagnostic report of one model
run_diagnostics("warfarin_final", output_dir = "warfarin_final/Diagnostics", save_output = TRUE)

# Two models side by side: parameters and criteria in common tables, paired performance plots
run_comparison("warfarin_example", "warfarin_final", label_1 = "Base", label_2 = "Final")
```

A complete walk-through on the Monolix warfarin demo is in `inst/tutorials/warfarin_tutorial.R`
(`system.file("tutorials", "warfarin_tutorial.R", package = "mlxflow")`).

Backward runs write one folder and one project per model, named
`run<NNN>_I<k>_<parameter>_<covariate>`: `NNN` is a run number that is never reused (like
Monolix's `ModelBuilding_NNN`), `k` the iteration, and the name gives the relation removed.
The verdict of each test is in `backward_log.csv`, never in a file name.

`alpha` is mandatory in the flowchart builders: it must be the level of the run being drawn
(for a Monolix search, the threshold of its settings), because the builders *replay* the
decisions, they do not take them.

## Help and arguments

Every exported function documents each of its arguments (`?backward_elimination`). In RStudio
the description of an argument appears in the completion pop-up while you type.

## Tests

`devtools::test()`. Monolix is replaced by a small fake engine, so the tests run anywhere.
The test fixtures are summaries of the Monolix warfarin demo project.
