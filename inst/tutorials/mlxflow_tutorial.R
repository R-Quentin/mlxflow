
#' @----------------------------------------------------------------------------
#' Tutorial for model building using R only (no GUI)
#' Combining LixofConnector with personal functions
#' Using a Monolix demos as examples
#' @----------------------------------------------------------------------------


# 1. Set-up: packages & helpers ================================================

# Lixoft connectors (make sure to restart session before running)
library(lixoftConnectors) 
initializeLixoftConnectors(software = "monolix") 

# Useful packages
# install.packages("kableExtra", "rmarkdown", "tidyverse")
library(kableExtra)
library(rmarkdown)
library(tidyverse)


# The package:
# install.packages("remotes")
# remotes::install_github("R-Quentin/mlxflow")
library(mlxflow)


# Directory
cur_dir <- this.path::here()
setwd(cur_dir)



# 2. Load, save and run project ================================================
# Load demo
project <- paste0(getDemoPath(), "/5.models_for_individual_parameters/5.2.covariate_model/warfarin_covariate1_project.mlxtran")
loadProject(projectFile = project)
# load_project_quiet(project)  # Same, without the Monolix notices on the console (checks the path first)

# Copy in working directory
saveProject(projectFile = "warfarin_example")
# getPopulationParameterInformation()


# Run estimations
# -> One by one run
# runPopulationParameterEstimation()
# runConditionalModeEstimation()
# runConditionalDistributionSampling()
# runStandardErrorEstimation()
# runLogLikelihoodEstimation()

# -> Batch run
getScenario() # Not all task enabled
scenario <- getScenario()
scenario$tasks = c(populationParameterEstimation   = T, 
                   conditionalModeEstimation       = T, 
                   conditionalDistributionSampling = T, 
                   standardErrorEstimation         = T, 
                   logLikelihoodEstimation         = T)
# scenario$linearization = TRUE
setScenario(scenario)
runScenario()

saveProject()
computeChartsData()



# 3. Explore model diagnostics  ================================================

#' @If not done already, run: loadProject(projectFile = "warfarin_final.mlxtran")

  ## Basic diagnostics    ------------------------------------------------------

# Observations
plotObservedData()

plotObservedData(
  # obsName = "concentration",
  settings = list(
    mean = T,
    error = T,
    meanMethod = "geometric",
    # timeAfterLastDose = T,
    # ylog = T,
    # lines = F,
    scales = "free_x",
    ncol = 2
  ),
  stratify = list(
    split = c("sex"),
    color = c("sex") #, mergedSplits = TRUE 
  )
)


# Prediction distribution
plotPredictionDistribution(settings = list(obs = TRUE))


# Individual fits
plotIndividualFits()

plotIndividualFits(
  stratify = list(individualSelection = list(indices = c(1, 8), isRange = TRUE)),
  settings = list(cens   = TRUE, 
                  popCov = TRUE,
                  legend = FALSE,
                  ylog   = TRUE,
                  xlab   = "Time (hr)",
                  ylab = "Concentration (ng/mL)"))

# Obs vs. Pred
plotObservationsVsPredictions()

plotObservationsVsPredictions(
  predictions = c("pop", "indiv"),
  stratify    = list(split = c("sex"), color = c("sex")),
  settings = list(ncol  = 1)
)


# VPC
plotVpc()

plotVpc(
  stratify = list(split = c("wt")),
  settings = list(timeAfterLastDose  = TRUE)
  )



# Population parameters 
getPopulationParameterInformation()


# Tests
getTests()



  ## Nicer diagnostics    ------------------------------------------------------

# -> Summary based diagnostics
ex_sum <- read_summary("warfarin_example/summary.txt", full = TRUE) 

ex_sum_tbl <- build_param_table(ex_sum, shrinkage_method = "CD", add_criteria = TRUE)

ex_sum_tbl |>
  select(-parameter) |> # = Use label
  kable_param_table(
    caption = sprintf("Parameter summary \u2014 %s", ex_sum$model_name),
    cols_to_remove = c("se_sa", "ci_low", "ci_high")
  ) |>
  # Post-processing aesthetics
  kable_styling(full_width = FALSE, font_size = 11, bootstrap_options = c("striped", "hover")) |>
  footnote(symbol = c("Conditional distribution")) |>
  column_spec(1, bold = TRUE)

# NB: Monolix RSE thresholds are: >50%, >100% and >200%
# The function supports full and partial estimation (eg no standard error estimation)


# -> FIM
plot_fisher_matrix(ex_sum, threshold = 0.6, ratio_warn = 20)

# -> Tests
plot_covariate_search("warfarin_example")


# -> Chats based diagnostics
# Needs ChartsData saved using computeChartsData()
obspred_res <- obspred_diagnostics("warfarin_example")
obspred_res$obspred # obspred_res$data / $metrics / $obspred / $aipe

plot_vpc("warfarin_example", facet_scales = "free", show_observations = TRUE)





  ## Generate full report ------------------------------------------------------

# Don't forget to save ChartsData using computeChartsData()
# Let's check:
fs::dir_ls("warfarin_example/ChartsData")                     # Short path if cur_dir = working directory
fs::dir_ls(file.path(cur_dir, "warfarin_example/ChartsData")) # Full path in case its afar


#' ----- run_diagnostics()
#' Function allowing you to run an Rmarkdown report template from the console.
#' @param workflow    : MANDATORY input for function. Diagnostic workflow (.Rmd)
#' @param knit_prefix : Optional input for function. Prefix used for html output (default is prefix)
#' @param knit_dir    : Optional input for function. Directory used for html output (default is directory of the input file)

#' @param model_dir   : MANDATORY input for Rmd. Desired model Monolix output directory (don't forget to export ChartsData)
#' @param save_output : Optional input for Rmd. When TRUE, save outputs in addition to knitted report
#' @param output_dir  : Optional input for Rmd. Directory used for saved outputs (plots & tables - defaults to model_dir)
#' @param prefix      : Optional input for Rmd. Prefix used for saved outputs names


# Most basic use:
run_diagnostics(model_dir = "warfarin_example")   # -> warfarin_example_report.html


# Tests (section "2.5 Covariate assessment" of the report) indicates that wt
# could be a covariate of V and CL, and sex a covariate of V as well.
# Let's build the covariate model.



# 4. Covariate model building   =================================================

#' @If not done already, run: loadProject(projectFile = "warfarin_example.mlxtran")

  ## a) Transform covs   -------------------------------------------------------

plotCovariates()
getCovariateInformation()$covariate |> 
  summary()

addContinuousTransformedCovariate(logtWt = "log(wt/70)")
addCategoricalTransformedCovariate(
  tSEX = list(
    from = "sex",
    transformed = list(G_0 = c("0"), G_1 = c("1")),
    reference = "G_0"))

getCovariateInformation()$covariate |> colnames()
saveProject()



  ## b1) SCM             -------------------------------------------------------

# -> Modify the settings
set_SCM <- getModelBuildingSettings()

{
  # Set method and thresholds
  set_SCM$strategy             <- 'scm' # "cossac" or "scm"
  set_SCM$criterion            <- "LRT" # "BIC" or "LRT"
  set_SCM$threshold$lrt[1]     <- 0.01  # forward p-value threshold for LRT
  set_SCM$threshold$lrt[2]     <- 0.001 # backward p-value threshold for LRT
  
  # Set covariates and parameters
  set_SCM$covariates <- c("sex", "logtWt")
  set_SCM$parameters <- c("Tlag", "ka", "V", "Cl")
  
  # Disable linearisation
  set_SCM$useLin = FALSE
}

# -> Run
runModelBuilding(settings = set_SCM)


# -> Get results
SCM_results <- getModelBuildingResults()
SCM_results |> saveRDS("SCM_results.rds")
# SCM_results <- readRDS("SCM_results.rds") # Reload without rerun

# Best model number
which(map_lgl(SCM_results, ~ .x$bestModel %in% TRUE)) 
# Gives either 7 or 23. 
# If you get 7, SCM stopped early. 
# If you get 23 SCM went further before coming back (backward)
# Anyway best relationship is : logtWt ~ V.


# Build table
SCM_results_table <- imap_dfr(SCM_results, \(res, idx) {
  
  ind_models <- if (!is.null(res$individualModels)) {
    list(res$individualModels)
  } else {
    list(NULL)
  }
  
  tibble(
    run              = idx,
    BICc             = res$BICc %||% NA_real_,
    LL               = res$LL %||% NA_real_,
    testedParameter  = res$tested[1] %||% NA_character_,
    testedCovariate  = res$tested[2] %||% NA_character_,
    bestModel        = res$bestModel %||% FALSE,
    individualModels = ind_models  # List
  )
})
head(SCM_results_table)


# -> Rename folder for future use
file.rename("warfarin_example/ModelBuilding", "warfarin_example/SCM")




  ## b2) COSSAC          -------------------------------------------------------

# -> Modify the settings
set_COSSAC <- getModelBuildingSettings()

{
  # Set method and thresholds
  set_COSSAC$strategy             <- 'cossac' # "cossac" or "scm"
  set_COSSAC$useSambaBeforeCossac <- FALSE    # if true, corresponds to "covSAMBA-COSSAC" in the GUI
  set_COSSAC$criterion            <- "LRT"    # "BIC" or "LRT"
  set_COSSAC$threshold$lrt[1]     <- 0.01     # forward p-value threshold for LRT
  set_COSSAC$threshold$lrt[2]     <- 0.001    # backward p-value threshold for LRT
  
  # Set covariates and parameters
  set_COSSAC$covariates <- c("sex", "logtWt")
  set_COSSAC$parameters <- c("Tlag", "ka", "V", "Cl")
  
  # Disable linearisation
  set_COSSAC$useLin = FALSE
}

# -> Run
runModelBuilding(settings = set_COSSAC)


# -> Get results
COSSAC_results <- getModelBuildingResults()
COSSAC_results |> saveRDS("COSSAC_results.rds")
# COSSAC_results <- readRDS("COSSAC_results.rds") # Reload without rerun


# Best model number
which(map_lgl(COSSAC_results, ~ .x$bestModel %in% TRUE)) # 2: logtWt ~ V

# Build table
COSSAC_results_table <- imap_dfr(COSSAC_results, \(res, idx) {
  
  ind_models <- if (!is.null(res$individualModels)) {
    list(res$individualModels)
  } else {
    list(NULL)
  }
  
  tibble(
    run              = idx,
    BICc             = res$BICc %||% NA_real_,
    LL               = res$LL %||% NA_real_,
    testedParameter  = res$tested[1] %||% NA_character_,
    testedCovariate  = res$tested[2] %||% NA_character_,
    bestModel        = res$bestModel %||% FALSE,
    individualModels = ind_models  # List
  )
})
head(COSSAC_results_table)

# -> Rename folder for future use
file.rename("warfarin_example/ModelBuilding", "warfarin_example/COSSAC")




  ## c) Plot Flowchart   -------------------------------------------------------
    ### 1/ SCM    --------------------------------------------------------------
forward_df <- read_models(file.path("warfarin_example/SCM"), full = FALSE)

# Strip categorical suffixes from covariate labels before building graph
forward_df <- forward_df |>
  mutate(covariates = map(covariates, clean_cov_labels))

# Construct SCM data
# Alpha MUST be the level used by the Monolix run to reflect the Monolix's decisions. 
# Builder can choose a different value above (via set_SCM$threshold$lrt[1] = 0.01).
forward_data <- build_scm_data(forward_df, alpha = 0.01)

# Get full model
full_model <- forward_df |>
  filter(model_name == forward_data$last_best_name) |>
  pull(covariates)

message("Last selected model: ", forward_data$last_best_name, 
        "\nFull covariate model: ", paste(unlist(full_model), collapse = ", "))

# If alpha coincides, mlxflow "$last_best_name" is the same as Monolix "$bestModel"

# However you can alaways cross-check Monolix's own report 
# (modelBuilding.txt in the model building folder)
mb_report <- read_model_building_report("warfarin_example/SCM")
mb_report$best_model_id

# If you had "23" as best model, this function find the "true" corresponding model
get_monolix_best_model(mb_report)
get_monolix_best_model(mb_report, forward_df)


# Visualisation
plot_flowchart(
  forward_data,
  title         = "Forward SCM - Warfarin PopPK",
  iter_label    = "Iteration",
  best_label    = "Best of iteration",
  node_size     = 30
)



    ### 2/ COSSAC --------------------------------------------------------------
# Read all models
cossac_df <- read_models(file.path("warfarin_example/COSSAC"), full = FALSE) |> 
  mutate(covariates = map(covariates, clean_cov_labels))

# Construct COSSAC data
cossac_data <- build_cossac_data(cossac_df, alpha = 0.01) # = set_COSSAC$threshold$lrt[1]


# Get full model
full_model <- cossac_df |>
  filter(model_name == cossac_data$last_best_name) |>
  pull(covariates)

message("Last selected model: ", cossac_data$last_best_name,
        "\nFull covariate model: ", paste(unlist(full_model), collapse = ", "))


# Visualisation + Short colour management (on the fly)
plot_flowchart(
  cossac_data,
  title = "COSSAC - Warfarin PopPK",
  node_size     = 30,
  colors = modifyList(
    eval(formals(plot_flowchart)$colors),  
    list(
      best  = list(fill = "palegreen3", border = "palegreen4", font = "black", arrow = "palegreen4")
    )
  ) 
)


# Alternative linear version
# plot_flowchart(scm_data = build_cossac_data_linear(cossac_df, alpha = 0.01))





# 5. Optional Backward deletion ================================================
  ## a) Automated function -----------------------------------------------------

# Backward is automated in Monolix SCM/COSSAC but you might want to check on your own.
# For example to :
# - Calculate LRT in case of COSSAC or for plotting reasons
# - Use a different elimination threshold
# - After manually deleting or adding a covariate


#' ----- backward_elimination():
#' Arguments and output files are documented in the header of backward_elimination()
#' In short: every model is saved as run<NNN>_I<k>_<parameter>_<covariate>   
#' (NNN = run number, never repeated; k = iteration)
#' and backward_log.csv holds every test (delta OFV, df, p-value, verdict).



# -> Best model selection
# Those are the best model numbers depending on model building strategy
# full_mod_path <- "warfarin_example/SCM/ModelBuilding_007/ModelBuilding_007.mlxtran"    # If SCM run
# full_mod_path <- "warfarin_example/COSSAC/ModelBuilding_002/ModelBuilding_002.mlxtran" # If COSSAC run

# But for demonstration purposes we will choose an over-parametrized model (with logtWt ~ Cl)
full_mod_path <- "warfarin_example/SCM/ModelBuilding_016/ModelBuilding_016.mlxtran"    # If SCM run
full_mod_path <- "warfarin_example/COSSAC/ModelBuilding_003/ModelBuilding_003.mlxtran" # If COSSAC run


# -> Output directory
backward_path <- "warfarin_example/Backward"


# -> Optional personalized scenario
# By default every model is estimated with the shortest valid chain (SAEM + CD
# sampling + log-likelihood) when running Lixoft model building.
# A scenario (same list format as for setScenario()) can be run instead for
# every model of the search: use it to get richer estimations for each one.
full_scenario <- getScenario()
full_scenario$tasks[c("populationParameterEstimation", "conditionalDistributionSampling",
                      "standardErrorEstimation", "logLikelihoodEstimation", "plots")] <- TRUE
 
# mlxflow backward elimination function uses the scenario of the model from
# which it originates: the shortest here as it comes from Monolix model building.
# However, it could have been a model with a "full" estimation scenario. In this
# case, you might want to reduce it to save computation time.



# -> Optional personalized diagnostic
# Any function(iter, project_file) can be used, e.g. to launch run_diagnostics() on each model.

# I've integrated this progress monitor as an example : 
# backward_monitor() prints one line after the base model and after each adopted
# removal: number of relations, OFV, BICc (and its change), largest
# eta-shrinkage (flagged above shrinkage_warn). 
# Run "backward_monitor" to print the function.
my_diag <- backward_monitor(shrinkage_warn = 30)


# -> Run backward check
result <- backward_elimination(
  project_file   = full_mod_path,
  alpha          = 0.01,          #' = OFV +6.63 for 1 df (see below)
  resume         = TRUE,          #' safe to re-launch after an interruption
  fast_screening = FALSE,         #' explicit search (SCM), change to TRUE for wald pre-ordering
  save_all       = TRUE,          #' save all model run, change to FALSE to keep best only
  # resume         = TRUE,          #' resume from last finished run in case of failure
  overwrite      = TRUE,          #' overwrite existing output_dir (did parameter changed?)
  # scenario       = full_scenario, #' optional: tasks run for every model (default: fastest valid chain)
  # diagnostic_fun = my_diag,       #' optional personalized diagnostic (verbose)
  output_dir     = backward_path
)

result$history       # step-by-step LRT summary (df, delta OFV, p-value) - step 0 = base model
result$stop_reason
result$final_project


# Just so you know:
# critical_ofv(alpha, df) = qchisq(1 - alpha, df)
# critical_ofv(0.01, 1) # 6.63
# critical_ofv(0.01, 2) # 9.21



  ## b) Plot Flowchart     -----------------------------------------------------

# ----- Backward Flowchart
backward_df <- read_models("warfarin_example/Backward", full = FALSE)

# Strip categorical suffixes from covariate labels before building graph
backward_df <- backward_df |>
  mutate(covariates = map(covariates, clean_cov_labels))

backward_data <- build_backward_data(backward_df, alpha = 0.01) # = alpha given to backward_elimination()

last_model <- backward_df |>
  filter(model_name == backward_data$last_best_name) |>
  pull(covariates)

message("Last selected model: ", backward_data$last_best_name,
        "\nLast covariate model: ", paste(unlist(last_model), collapse = ", "))


# Default look: relations kept (significant) are grey, those that are not significant are white
# with a dark red outline, and the relation dropped at each iteration is red.
# Labels and colours follow the type of the search, so the call stays short:
plot_flowchart(
  backward_data,
  title      = "Backward Elimination - Warfarin PopPK"
)


# -> Example personalized plot with reusable colour management: 
# Start from the default palette and change one category
my_backward_colors <- eval(formals(plot_flowchart_combined)$colors_backward)
my_backward_colors$sig <- list(fill = "skyblue3", border = "skyblue4", 
                               font = "black",    arrow  = "skyblue4") # Colours significant (kept) variables

# Change labels
plot_flowchart(
  backward_data,
  title         = "Backward Covariates Deletion for Warfarin Model",
  iter_label    = "Run",
  base_label    = "Full model \n(V ~ Wt, Cl ~ Wt)",
  best_label    = "Not significant (dropped)",        # Legend labels can be changed as well
  sig_label     = "Significant: p < {alpha}) (kept)", # {alpha} is replaced by the level of the test
  colors        = my_backward_colors
)




# ----- Combined Flowchart
# Get a full forward + backward flowchart plot

# SCM
merged_data_SCM <- merge_scm_data(forward_data, backward_data, gap = 2)
# str(merged_data_SCM, max.level = 2) # Check 

plot_flowchart_combined(merged_data_SCM, title = "Covariate Model Building - Warfarin PopPK")


# COSSAC + Reuse colour management (optional)
merged_data_COSSAC <- merge_scm_data(cossac_data, backward_data, gap = 2)

plot_flowchart_combined(
  merged_data_COSSAC,
  title              = "Covariate Model Building - Warfarin PopPK",
  forward_title      = "Forward inclusion (COSSAC)",
  colors_backward    = my_backward_colors,
  node_size          = 30
)




# 6. Final covariate model      ================================================

# Final model is V ~logWT 
final_mod_path <- file.path(cur_dir, "warfarin_example/Backward/run002_I1_Cl_logtWt.mlxtran")

# Let's run the full estimation and full report
loadProject(projectFile = final_mod_path)
saveProject(projectFile = "warfarin_final")

# Run missing estimations
getScenario() # SE and plots disabled
setScenario(full_scenario) # See above in 6.a)
runScenario()
saveProject()
computeChartsData()

# Final diagnostics - Let's save all tables and figures to a designated folder
run_diagnostics(
  model_dir   = "warfarin_final",
  output_dir  = "warfarin_final/Diagnostics",
  save_output = TRUE, 
  prefix      = "warfarin_FINAL"
)


# Let's compare our final model to where we started
# Side by side: parameters, information criteria (common tables) 
# and predictive performance (paired plots).
# Covariates and convergence stay in run_diagnostics().
run_comparison(
  "warfarin_example", "warfarin_final", 
  label_1 = "Base",   label_2 = "Final",            # Labels used in tables and plots
  knit_prefix = "Warfarin_structural_VS_covariate") # Report title ("/Warfarin_structural_VS_covariate_comparison.html")


# Templates for this functions are included in the package but you can use your
# own via the "workflow argument".
# /!\ A template must declare the parameters 'model_dir', 'save_output',
# 'output_dir' and 'prefix'.




# 7. Model comparison           ================================================

# Apart from the comparison template provided (see above), you might want to use
# the functions independently or compare more than two models at once.

  ## a) Two models comparison      ---------------------------------------------
    ### 1/ Situation I  --------------------------------------------------------
# Backward deletion found logtWt ~ Cl should be deleted but you have a doubt.
# Let's compare the model with and without

# If SCM run
model_with    <- read_summary("warfarin_example/SCM/ModelBuilding_016/ModelBuilding_016/summary.txt", full = TRUE)
model_without <- read_summary("warfarin_example/SCM/ModelBuilding_007/ModelBuilding_007/summary.txt", full = TRUE)
# If COSSAC run
model_with    <- read_summary("warfarin_example/COSSAC/ModelBuilding_003/ModelBuilding_003/summary.txt", full = TRUE)
model_without <- read_summary("warfarin_example/COSSAC/ModelBuilding_002/ModelBuilding_002/summary.txt", full = TRUE)



# -> Here we compare only 1 column (estimate) and compute the differences 

# compare_two_models():
# Fully join TWO models and compare one or several columns (delta columns computed)
#     - first and second : models (read_summary(..., full = TRUE) objects, summary.txt paths or model folders)
#     - values           : column(s) to compare (estimate, rse_pct etc.). Comparison following : second - first
#     - clean            : if TRUE (default), clean covariates label (beta_CL_WT -> ~WT)
#     - delta            : differences shown in the table: c("abs", "rel") (default, = TRUE), "abs", "rel" or FALSE (none)
#                          Only the table changes: $data always holds both.

one_param_two_mod <- compare_two_models(
  first  = model_with, 
  second = model_without, 
  values = "estimate", 
  clean  = TRUE
  )

# A list: the data and a ready-made table
one_param_two_mod$data    # raw table
one_param_two_mod$kable   # styled table (the footnote gives the two models)
one_param_two_mod         # both are printed

# The table is a kableExtra object, so it can be styled further
one_param_two_mod$kable |> 
  column_spec(5, italic = TRUE) |> 
  kable_styling(font_size = 9)

# Without the differences (or with only one of them)
compare_two_models(model_with, model_without, delta = FALSE)$kable
compare_two_models(model_with, model_without, delta = "rel")$kable


# You can also use:
compare_criteria(first = model_with, second = model_without)
# (compare_obspred_metrics() compares the prediction errors of models: see situation B)



    ### 2/ Situation II --------------------------------------------------------
# What changed from the initial and final model ?

base  <- read_summary("warfarin_example/summary.txt", full = TRUE)
final <- read_summary("warfarin_final/summary.txt", full = TRUE)

# -> Now we look at more thant one column (and compute the differences)

#' ----- compare_two_models():
#' Fully join TWO models and compare one or several columns (delta columns computed)
#'     - first and second : models (read_summary(..., full = TRUE) objects, summary.txt paths or model folders)
#'     - values           : column(s) to compare (estimate, rse_pct etc.). Comparison following : second - first
#'     - clean            : if TRUE (default), clean covariates label (beta_CL_WT -> ~WT)
#'     - bold_above       : values above a threshold are printed in bold in the table (here RSE > 30 %)
#'     - delta            : see situation A


# -> Estimated parameters
all_param_two_mod <- compare_two_models(
  first      = base, 
  second     = final, 
  values     = c("estimate", "rse_pct"), 
  clean      = TRUE,
  bold_above = c(rse_pct = 30)
  )

all_param_two_mod


# -> Prediction errors (read from their ChartsData)
# compare_obspred_metrics():
#     - first and second : model folders, read_summary() objects or obspred_diagnostics() outputs
#     - labels           : names of the two models
#     - format           : "wide" (default): one row per model (and "delta") with the metrics in columns,
#                          "long": one row per metric, one column per model
#     - split            : keep the strata of the Monolix `split` when both models have them (default).
#                          If only one has, the observations are pooled for both (a message tells so)
#     - y_type_labels, split_labels : rename observed variables / strata so that the two models agree
two_mod_metrics <- compare_obspred_metrics(first = base, second = final, labels = c("Base", "Final"))
two_mod_metrics

# Can also be formatted through kableExtra:
two_mod_metrics |>
  select(-DV, -Y_TYPE) |> 
  select(-Prediction) |> 
  kable(
    digits   = 2,
    caption  = "Performance changes",
    booktabs = TRUE
  ) |> 
  row_spec(which(two_mod_metrics$Model == "delta"), bold = TRUE) |> 
  auto_pack_rows(two_mod_metrics, group_vars = "Prediction") |> 
  kable_styling(font_size = 10, bootstrap_options = c("striped", "bordered", "hover"))

# Long format, if you prefer to reshape it yourself
# compare_obspred_metrics(first = base, second = final, labels = c("Base", "Final"), format = "long")



  ## b) Multiple models comparison ---------------------------------------------

# Just for example : COSSAC workflow models comparison 

# A real situation could be : Is IOV useful? Try IOV implementation on each
# parameter and compare. See "Testing variability (IOV / IIV)" section.


# One column
one_param_all_mod <- compare_models(dir = "warfarin_example/COSSAC", values = "estimate", clean = TRUE)
one_param_all_mod


# Multiple columns
all_param_all_mod <- compare_models(
  dir        = "warfarin_example/COSSAC", 
  values     = c("estimate", "shrinkage"), 
  clean      = TRUE,
  bold_above = c(shrinkage = 50))
all_param_all_mod

# Differences to a reference model (the first one by default, or `reference = "<folder name>"`)
# delta = TRUE / "abs" / "rel" adds them next to each model (FALSE by default)
compare_models(dir = "warfarin_example/COSSAC", values = "estimate", clean = TRUE, delta = "abs")


# 8. Convergence assessment     ================================================

#' @If not done already, run: loadProject(projectFile = "warfarin_final.mlxtran")

# runAssessment()

# -> Personalized settings
{
  set <- getAssessmentSettings()
  set$nbRuns             <- 10
  set$nextToProject      <- FALSE
  set$extendedEstimation <- TRUE
  set$useLin             <- FALSE
  set$initialParameters  <- data.frame(
    parameters = c("ka_pop", "V_pop", "Cl_pop"),
    min = c(0.1, 0.1, 0.01), max = c(10, 10, 1),
    fixed = FALSE)
}
runAssessment(settings = set)

AssessmentResults <- getAssessmentResults()

# -> Save the results to keep them with the model (they only exist in the R session otherwise)
saveRDS(AssessmentResults, "warfarin_final/AssessmentResults.rds")

# -> Plots: the lines of the runs should end on the same values
plot_convergence_assessment(AssessmentResults, which = "iterations") 

# Final estimate of each run (ci = TRUE adds the 95 % interval)
plot_convergence_assessment(AssessmentResults, which = "estimates", ci = FALSE) 


# The reshaped data, if you want to make your own plots/diagnostics
convergence_data <- tidy_convergence_assessment(AssessmentResults)
# Ex:
convergence_data$final |>
  group_by(parameter) |>
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    mean_rse_pct  = mean((se / estimate) * 100, na.rm = TRUE),
    cv_pct        = (sd(estimate, na.rm = TRUE) / mean(estimate, na.rm = TRUE)) * 100 
    )


# -> Reading the assessment back, without the R session
# The model folder (or its Assessment folder) can be given directly:
plot_convergence_assessment("warfarin_final", which = "estimates", ci = TRUE)
# convergence_data <- tidy_convergence_assessment("warfarin_final")
# AssessmentResults <- read_assessment("warfarin_final")

# NB: the run folders hold the final estimates only (the SAEM trajectories are not written there
# unless the charts of each run exist). "iterations" is drawn when the AssessmentResults.rds saved
# above is in the model folder or in its Assessment folder, which is why it is worth keeping.
plot_convergence_assessment("warfarin_final", which = "iterations")


# -> It can also be a section of the model report let's re-run diagnostics:
# assessment = NULL (default): uses <model_dir>/Assessment, or <model_dir>/AssessmentResults.rds
run_diagnostics("warfarin_final")
# or choose it: assessment = "warfarin_final/AssessmentResults.rds" (or a folder, or the object)



# 9. Testing variability        ================================================

# test_variability() adds a random effect to one parameter at a time and compares each model
# with the reference by a likelihood ratio test (LRT).
#   - level = an occasion level (e.g. "OCC") : inter-occasion variability (IOV)
#   - level = "id"                            : inter-individual variability (IIV), i.e. which
#                                               parameters need an omega
# Arguments and outputs are documented in the header of test_variability(). In short:
#   - forward = FALSE (default): one pass, every parameter against the base model.
#     forward = TRUE: stepwise selection like backward_elimination(), until nothing is
#     significant at `alpha`.
#   - Every model is saved as run<NNN>_I<k>_<parameter>_<level> with its own results (the
#     starting project is never overwritten); variability_log.csv holds every test.
#   - scenario: same as in backward_elimination() (fast or complete estimation of every model).

  ## a) IOV Example ------------------------------------------------------------
    ### 1/ Model         -------------------------------------------------------

# -> Load project
# We'll use a different example (not the warfarin one) to try IOV implementation
IOV_project <- paste0(getDemoPath(), "/5.models_for_individual_parameters/5.4.inter_occasion_variability/iov3_project.mlxtran")
loadProject(projectFile = IOV_project)

# Inter-Occasion Variability is already enabled on V:
getIndividualParameterModel()$variability

# Let's remove it for demonstration purposes:
setIndividualParameterVariability(OCC = list(V = FALSE))
getIndividualParameterModel()$variability

# Copy in working directory
saveProject(projectFile = "IOV_example")


# -> Check available parameters and variability levels
getIndividualParameterModel()$name
getIndividualParameterModel()$variability |> names()


# -> Personalized scenario
fast_scenario <- list(
  tasks = c(
    populationParameterEstimation       = TRUE,
    conditionalDistributionSampling     = TRUE,
    conditionalModeEstimation           = FALSE,
    standardErrorEstimation             = TRUE,
    logLikelihoodEstimation             = TRUE,
    plots                               = FALSE
  ),
  linearization = FALSE,
  plotList = c(
    "obspred", 
    "vpc", 
    "saemresults"
  )
)

    ### 2/ One pass test -------------------------------------------------------

# -> One pass: IOV on each parameter against the base model
IOV_one_pass <- test_variability(
  project_file = "IOV_example.mlxtran",
  output_dir   = "IOV_example/IOV_one_pass",
  # parameters   = c("ka", "V",  "Cl"), # Default = all
  level        = "OCC", # Default = first variability level that is not the individual one
  alpha        = 0.05,
  scenario     = fast_scenario
)

IOV_one_pass$runs            # one row per model: OFV, delta OFV, df, p-value, verdict
IOV_one_pass$final_project   # = the base model in a single pass: nothing is adopted
IOV_one_pass |> saveRDS("IOV_example/IOV_one_pass/IOV_one_pass.rds")
# IOV_one_pass <- read_rds("IOV_example/IOV_one_pass/IOV_one_pass.rds")



# -> Check results
# OFV
IOV_df <- read_models(file.path("IOV_example/IOV_one_pass"), full = FALSE)

IOV_df |> 
  mutate(
    delta_OFV  = OFV - OFV[model_id == 0],
    delta_BICc = BICc - BICc[model_id == 0],
  ) |> 
  select(model_name, OFV, delta_OFV, AIC, BIC, BICc, delta_BICc) |> 
  arrange(OFV)


# Parameters (differences to the base model, the first one)
IOV_comp <- compare_models(
  dir        = "IOV_example/IOV_one_pass", 
  values     = c("estimate", "rse_pct", "shrinkage"), 
  clean      = TRUE, 
  delta      = FALSE,
  bold_above = c(rse_pct = 30, shrinkage = 50)
)
IOV_comp


# Predictions
compare_obspred_metrics(
  first = "IOV_example/IOV_one_pass/run000_I0_base", 
  second = "IOV_example/IOV_one_pass/run003_I1_Cl_OCC",
  labels = c("Base", "IOV on Cl")
  )

op1 <- obspred_diagnostics("IOV_example/IOV_one_pass/run000_I0_base")
op2 <- obspred_diagnostics("IOV_example/IOV_one_pass/run003_I1_Cl_OCC")
op1$obspred | op2$obspred


# -> Conclusion : IOV doesn't improve the model



    ### 3/ Forward test  -------------------------------------------------------

# -> Forward selection: adopt the most significant IOV, then test the others on top of it
IOV_forward <- test_variability(
  project_file = "IOV_example.mlxtran",
  output_dir   = "IOV_example/IOV_forward",
  level        = "OCC",
  alpha        = 0.05,
  scenario     = fast_scenario,
  forward      = TRUE, 
  prune_alpha  = 0.1, 
  max_rse      = 100,
  # resume       = TRUE,   # safe to launch again after an interruption
  warm_start  =  TRUE
)

IOV_forward$history         # base model, then each IOV added
IOV_forward$stop_reason
IOV_forward$final_project

# -> Flowchart of the search 
# (same plot_flowchart() as for covariates; the data come from the result or
# from the output folder, also for a run that was interrupted)
plot_flowchart(
  build_variability_data(IOV_forward),     # or build_variability_data("IOV_example/IOV_forward")
  title = "IOV forward selection - OCC"
)




  ## b) IIV Example ------------------------------------------------------------

# Can also be used to test IIV : does a parameter need an omega at all?

# -> Load project
# We'll use a different example (not the warfarin one) to try IIV implementation
IOV_project <- paste0(getDemoPath(), "/5.models_for_individual_parameters/5.4.inter_occasion_variability/iov3_project.mlxtran")
loadProject(projectFile = IOV_project)

# Let's remove the variability (IIV/IOV) for demonstration purposes:
setIndividualParameterVariability(id = list(ka = FALSE, V = FALSE, Cl = FALSE), OCC = list(V = FALSE))
getIndividualParameterModel()$variability

# Copy in working directory
saveProject(projectFile = "IIV_example")


# -> Forward selection: adopt the most significant IIV, then test the others on top of it
IIV_test <- test_variability(
  project_file = "IIV_example.mlxtran",
  output_dir   = "IIV_example/IIV_test",
  level        = "id",
  alpha        = 0.05,
  scenario     = fast_scenario,
  forward      = TRUE, 
  prune_alpha  = 0.1, 
  max_rse      = 100,
  resume       = TRUE,   # safe to launch again after an interruption
  warm_start   = TRUE, 
  overwrite    = TRUE
)

IIV_test$history         # base model, then each IIV added
IIV_test$stop_reason
IIV_test$final_project


# -> Flowchart of the search 
# (same plot_flowchart() as for covariates; the data come from the result or
# from the output folder, also for a run that was interrupted)
plot_flowchart(
  build_variability_data(IIV_test),     # or build_variability_data("IIV_example/IIV_test")
  title = "IIV forward selection"
)

compare_two_models(
  first  = "IIV_example/IIV_test/run000_I0_base",  
  second = "IIV_example/IIV_test/run006_I3_ka_id", 
  values = c("estimate", "rse_pct"), 
  delta  = "abs", 
  bold_above = c(rse_pct = 30)
)


