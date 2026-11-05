# mlxflow 0.0.0.9000

First assembly of the package from three scripts (engine, readers/diagnostics, flowcharts).

## New

* `run_comparison()` and the `model_comparison` report template: two models side by side
  (information criteria and parameters in common tables; prediction-error metrics,
  observations vs predictions, VPC and NPDE as paired plots).
* `compare_criteria()` and `compare_obspred_metrics()`, the building blocks of that report.
* `read_summary()` and `read_models()` expose `ll_method` (importance sampling or linearization);
  `read_models()` and `compare_criteria()` warn when models computed differently are mixed, since
  their OFV, AIC and BIC are not comparable.
* `read_model_building_report()` and `get_monolix_best_model()`: the model retained by Monolix,
  read from `modelBuilding.txt`.
* `auto_pack_rows()`, `check_chartsdata()`, `resolve_dv()` and `read_iparam_pdf()` are exported
  (the report templates use them).
* `inst/tutorials/warfarin_tutorial.R`: the tutorial of the package, on the Monolix warfarin demo.

## Renamed functions

| Before | Now |
|---|---|
| `read_all_models_1()`, `read_all_models_2()` | `read_models(layout = "auto")` |
| `compare_2_summaries_1()`, `compare_2_summaries_2()` | `compare_two_models(values = )` |
| `compare_all_summaries_1()`, `compare_all_summaries_2()` | `compare_models(values = )` |
| `plot_scm_flowchart()` | `plot_flowchart()` |
| `plot_combined_flowchart()` | `plot_flowchart_combined()` |
| `run_my_diagnostics()` | `run_diagnostics()` |
| `covariate_effects_poppk()` | `compute_covariate_effects()` |
| `rename_folder_b()` | `renumber_model_folder()` |
| `ofv_threshold()` | `critical_ofv()` |
| `match_monolix_best()` | `read_model_building_report()` + `get_monolix_best_model()` |

Removed: `compare_summaries_1/2/3`, `build_scm_data_V1/V2` (duplicates of newer functions).

## Behaviour changes

* `alpha` is required by `build_scm_data()`, `build_cossac_data()` and `build_backward_data()`;
  `plot_flowchart()` takes the level from the data (its `alpha` argument is gone).
* Significance uses the degrees of freedom of each test (coefficients added or removed), not 1.
* Backward runs: new folder naming (`run<NNN>_I<k>_...`), `save_all = TRUE` by default,
  `max_steps` is now `max_iter`, `backward_steps.csv` is now `backward_iterations.csv`.
* `compare_two_models()`: the one-value case now returns `delta_<value>` / `delta_rel_<value>`.
* `run_diagnostics()`: the html report is named after the model (it used to be a generic
  `model_checkup.html`) and is written in the working directory unless `knit_dir` is given.
* `model_name` of `read_models()` is the folder name.

## Bug fixes

* `read_summary()`: covariates are read from the population block only; the first column of the
  Fisher matrix was read as an extra, glued covariate when a beta had the widest label.
* Comparisons work on nested (ModelBuilding) folders and on models without standard errors.
* Functions no longer depend on which packages the session has attached (`fixed()`, `union()`,
  `footnote()`, igraph functions).
* Single-node flowcharts are drawn.
