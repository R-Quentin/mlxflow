# mlxflow 0.1.1

First assembly of the package from three scripts (engine, readers/diagnostics, flowcharts).

## Changes (this round)

* `test_variability()` (replaces the tutorial helper `test_IOV()`): tests the addition of a random
  effect on each parameter, at an occasion level (IOV) or at `"id"` (IIV), by LRT. `forward = TRUE`
  selects the effects stepwise like `backward_elimination()`. Every model has its own result folder
  (the starting project is no longer overwritten, which caused the "old inconsistent project" warnings),
  the default level is the first one that is not `"id"`, and a model with `"id"` only asks for
  `level = "id"` explicitly.
* `backward_elimination()` and `test_variability()` take a `scenario`: the tasks of a Monolix scenario
  are run for every model (fast, or complete with standard errors and plots) instead of the fixed chain.
* `load_project_quiet()` is exported: loads a project without the Monolix notices, checks the path first
  and shows what Monolix printed if it fails. Internal callers are unchanged.
* `compare_two_models()` and `compare_models()` take `delta` (`TRUE`, `FALSE`, `"abs"`, `"rel"`) to choose
  the differences shown in the table; `compare_models()` can now add the differences to a `reference` model.
* `kable_param_comparison()`: parameter comparison table with estimates, RSE and shrinkage, their change,
  and the bold / colour thresholds of `kable_param_table()`. Used by the comparison report.
* `compare_obspred_metrics()` returns the metrics in columns by default (`format = "wide"`, `"long"` gives
  the previous layout) and handles strata (`split`) correctly instead of listing them as a metric; strata
  or observed variables that only one model has are reported, and a single observed variable named
  differently in the two models is matched.
* The backward flowchart: kept relations in grey, non-significant ones in white with a dark red outline, the
  dropped relation in red, with new legend texts (`not_sig_label`, `{alpha}` in labels). The combined chart
  no longer draws the border of significant nodes with the non-significant colour.
* The diagnostic report is now written as `<model>_report.html` (it was `<model>_model_checkup.html`), the
  comparison report as `<prefix>_comparison.html` (it was `<prefix>_model_comparison.html`).
* `test_variability()` can be resumed (`resume = TRUE`): the state is saved after every estimated model
  (`variability_state.rds`) and the settings of the interrupted run are kept. Like `backward_elimination()`,
  everything is written in `output_dir`.
* `test_variability()`: an OFV that is not a finite number (`-Inf`, `NaN`: importance sampling that breaks
  down) is flagged `invalid`, kept in the table and never adopted, instead of stopping the search with an
  error. `max_rse` requires the new random effect to be well estimated to be adopted (needs a scenario with
  standard errors), and `prune_alpha` drops, after the first iteration, the parameters that showed nothing
  against the base model (pre-selection). New columns `rse_new` and `eligible` in the runs table and the log.
* Convergence assessment: the standard errors of `getAssessmentResults()` come as a vector sorted alphabetically
  by parameter name, not in the order of the convergence columns; they were attributed to the wrong
  parameters in `tidy_convergence_assessment()` / `plot_convergence_assessment(ci = TRUE)`, and are now matched
  by name. `read_assessment()` only reads the `Assessment_NNN` runs (the model's own results were counted as an
  extra run when the model folder was given), accepts the model folder, and uses an `AssessmentResults.rds` found
  next to the runs or in the model folder (it holds the SAEM trajectories). `run_diagnostics()` also looks for
  `<model_dir>/AssessmentResults.rds`.
* `build_variability_data()` turns a `test_variability()` result (or its output folder, even for an
  interrupted run) into the data of `plot_flowchart()`: the flowchart of an IOV / IIV test, with the adopted
  model of each iteration highlighted and invalid, failed or imprecise models labelled.
* Charts data (`test_variability()`, `backward_elimination(export_charts = TRUE)`): `computeChartsData()` is
  called once the project is saved (the charts are only written to a saved project).
* `read_assessment()` reads the `Assessment` folder written by Monolix (estimates and standard errors of each
  run, SAEM trajectories when the charts of the runs exist) into the structure of `getAssessmentResults()`.
  `plot_convergence_assessment()` and `tidy_convergence_assessment()` accept that folder, and
  `run_diagnostics()` uses `<model_dir>/Assessment` by default when it exists.
* `compare_obspred_metrics()`: `match_strata = TRUE` compares the strata of one model (Monolix `split`, e.g. by
  route) with the observed variables of the other (one error model per route: y1 / y2), after renaming them
  alike with `split_labels` / `y_type_labels`.
* `kable_param_comparison()` shows absolute changes for the estimates, RSE and shrinkage.
* `read_covariate_assessment()`, `plot_covariate_search()` and `covariate_tests_table()` read bounded p-values
  such as `<2.2e-16` in the random effect x covariate tests (the column was text: error in `log10()`).
* The tutorial is now `inst/tutorials/mlxflow_tutorial.R` (its IOV part uses another demo model).

## New

* `compare_two_models()` now returns a list like `compare_models()`: `data` (the tibble it used to
  return) and `kable` (a ready-made table). Both accept `bold_above = c(rse_pct = 30)` to print the
  values above a threshold in bold.
* `compare_obspred_metrics()` accepts model folders and `read_summary()` objects (`read_summary()`
  now keeps `summary_path`), and explains what it cannot use.
* `backward_monitor()`: a ready-made `diagnostic_fun` for `backward_elimination()` (relations, OFV,
  BICc and its change, largest eta-shrinkage).
* Convergence assessment: `tidy_convergence_assessment()`, `plot_convergence_assessment()` and an
  optional section of the check-up report (`run_diagnostics(assessment = )`).
* The check-up report shows its run configuration as a table and no longer prints `NULL` for
  analyses that were not run; the comparison report fixes its metrics table.
* `tests/testthat/test-templates.R` renders both reports end to end when a model with ChartsData is
  placed in `tests/testthat/fixtures/warfarin_final/` (skipped otherwise).

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
  `<model>_model_checkup.html`) and is written in the working directory unless `knit_dir` is given.
* `model_name` of `read_models()` is the folder name.

## Bug fixes

* `read_summary()`: covariates are read from the population block only; the first column of the
  Fisher matrix was read as an extra, glued covariate when a beta had the widest label.
* Comparisons work on nested (ModelBuilding) folders and on models without standard errors.
* Functions no longer depend on which packages the session has attached (`fixed()`, `union()`,
  `footnote()`, igraph functions).
* Single-node flowcharts are drawn.
