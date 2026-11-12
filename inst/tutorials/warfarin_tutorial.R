
#' @----------------------------------------------------------------------------
#' Tutorial for model building using R only (no GUI)
#' Combining LixofConnector with personnal functions
#' Using a Monolix demo (warfarin_covariate1_project) as example
#' @----------------------------------------------------------------------------


# Set-up: packages & helpers ===================================================
# Packages
library(lixoftConnectors) 
initializeLixoftConnectors(software = "monolix")

library(tidyverse)
library(kableExtra)
library(rmarkdown)

# The package: readers, tables, diagnostics, comparisons, backward elimination, flowcharts
# install.packages("remotes")
# remotes::install_github("R-Quentin/mlxflow")
library(mlxflow)


# Directory
cur_dir <- this.path::here()
setwd(cur_dir)



# Load, save and run project ===================================================
# Load demo
project <- paste0(getDemoPath(), "/5.models_for_individual_parameters/5.2.covariate_model/warfarin_covariate1_project.mlxtran")
loadProject(projectFile = project)

# Copy in working directory
saveProject(projectFile = "warfarin_example")


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

computeChartsData()
saveProject()



# Explore model diagnostics  ===================================================

#' @If not done already, run: loadProject(projectFile = "warfarin_example.mlxtran")

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





  ## Nicer diagnostics    ------------------------------------------------------

# -> Summary based diagnostics
ex_sum <- read_summary("warfarin_example/summary.txt", full = TRUE) 

ex_sum_tbl <- build_param_table(ex_sum, shrinkage_method = "CD", add_criteria = TRUE) |>
  clean_param_labels(overwrite = FALSE)

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


plot_fisher_matrix(ex_sum, threshold = 0.6, ratio_warn = 20)


# -> Chats based diagnostics
obspred_res <- obspred_diagnostics("warfarin_example")
obspred_res$obspred # obspred_res$data / $metrics / $obspred / $aipe




  ## Generate full report ------------------------------------------------------

# Charts data need to be saved for diagnostic run !
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
run_diagnostics(model_dir = "warfarin_example")


# In "2.5 Covariate assessment" section, it looks like wt could be a covariate
# of V and CL, and sex a covariate of V as well.
# Let's build the covariate model.



# Covariate model building   ===================================================

#' @If not done already, run: loadProject(projectFile = "warfarin_example.mlxtran")

  ## Transform covs  -----------------------------------------------------------

plotCovariates()
getCovariateInformation()$covariate |> 
  summary()

addContinuousTransformedCovariate(logtWt = "log(wt/70)")
saveProject()

# getPopulationParameterInformation()


  ## SCM             -----------------------------------------------------------

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




  ## COSSAC          -----------------------------------------------------------

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




  ## Plot Flowchart  -----------------------------------------------------------
    ### ----- SCM -----
forward_df <- read_models(file.path("warfarin_example/SCM"), full = FALSE)

# Strip categorical suffixes from covariate labels before building graph
forward_df <- forward_df |>
  mutate(covariates = map(covariates, clean_cov_labels))

# Construct SCM data
# Alpha MUST be the level used by the Monolix run to reflect the Monolix's decisions. 
# Builder can choose a different value above (via set_SCM$threshold$lrt[1] = 0.01).
forward_data <- build_scm_data(forward_df, alpha = 0.01)

# Get full model
full_model <- forward_df %>%
  filter(model_name == forward_data$last_best_name) %>%
  pull(covariates)

message("Last selected model: ", forward_data$last_best_name, 
        "\nFull covariate model: ", paste(unlist(full_model), collapse = ", "))

# If alpha coincides, mlxflow "$last_best_name" is the same as Monolix "$bestModel"

# However you can alaways cross-check Monolix's own report 
# (modelBuilding.txt in the model building folder)
mb_report <- read_model_building_report("warfarin_example/SCM")
get_monolix_best_model(mb_report, forward_df)$model_name


# Visualisation
plot_flowchart(
  forward_data,
  title         = "Forward SCM - Warfarin PopPK",
  iter_label    = "Iteration",
  best_label    = "Best of iteration",
  node_size     = 30
)



    ### ----- COSSAC ----
# Read all models
cossac_df <- read_models(file.path("warfarin_example/COSSAC"), full = FALSE) |> 
  mutate(covariates = map(covariates, clean_cov_labels))

# Construct COSSAC data
cossac_data <- build_cossac_data(cossac_df, alpha = 0.01) # = set_COSSAC$threshold$lrt[1]


# Get full model
full_model <- cossac_df %>%
  filter(model_name == cossac_data$last_best_name) %>%
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





# Optional Backward deletion ===================================================
  ## Automated function --------------------------------------------------------

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
# Those are the best model number depending on model building strategy
# full_mod_path <- "warfarin_example/SCM/ModelBuilding_007/ModelBuilding_007.mlxtran"    # If SCM run
# full_mod_path <- "warfarin_example/COSSAC/ModelBuilding_002/ModelBuilding_002.mlxtran" # If COSSAC run

# But for demonstration purposes we will choose an over-parametrized model (with logtWt ~ Cl)
full_mod_path <- "warfarin_example/SCM/ModelBuilding_016/ModelBuilding_016.mlxtran"    # If SCM run
full_mod_path <- "warfarin_example/COSSAC/ModelBuilding_003/ModelBuilding_003.mlxtran" # If COSSAC run


# -> Output directory
backward_path <- "warfarin_example/Backward"


# -> Optional progress monitor
# backward_monitor() prints one line after the base model and after each adopted removal: number of
# relations, OFV, BICc (and its change), largest eta-shrinkage (flagged above shrinkage_warn).
# Any function(iter, project_file) can be used instead, e.g. to launch run_diagnostics() on each model.
my_diag <- backward_monitor(shrinkage_warn = 30)

# -> Run backward check
result <- backward_elimination(
  project_file   = full_mod_path,
  alpha          = 0.01,     #' = OFV +6.63 for 1 df (see below)
  resume         = TRUE,     #' safe to re-launch after an interruption
  fast_screening = FALSE,    #' explicit search (SCM), change to TRUE for wald pre-ordering
  save_all       = TRUE,     #' save all model run, change to FALSE to keep best only
  # resume         = TRUE,     #' resume from last finished run in case of failure
  overwrite      = TRUE,     #' overwrite existing output_dir (parameter changed?)
  diagnostic_fun = my_diag,
  output_dir     = backward_path
)

result$history       # step-by-step LRT summary (df, delta OFV, p-value) - step 0 = base model
result$stop_reason
result$final_project


# Just so you know:
# critical_ofv(alpha, df) = qchisq(1 - alpha, df)
# critical_ofv(0.01, 1) # 6.63
# critical_ofv(0.01, 2) # 9.21



  ## Plot Flowchart     --------------------------------------------------------

# ----- Backward Flowchart
backward_df <- read_models("warfarin_example/Backward", full = FALSE)

# Strip categorical suffixes from covariate labels before building graph
backward_df <- backward_df |>
  mutate(covariates = map(covariates, clean_cov_labels))

backward_data <- build_backward_data(backward_df, alpha = 0.01) # = alpha given to backward_elimination()

# head(backward_data$nodes |> select(model_id, iter, added_cov, delta_OFV, significant, is_best))

last_model <- backward_df %>%
  filter(model_name == backward_data$last_best_name) %>%
  pull(covariates)

message("Last selected model: ", backward_data$last_best_name,
        "\nLast covariate model: ", paste(unlist(last_model), collapse = ", "))


# Example of large colour management (semi-permanent)
{
  backward_colors <- eval(formals(plot_flowchart)$colors)
  backward_colors$base    <- list(fill = "green4", border = "darkgreen", font = "white") # Colours of forward last best (base model)
  backward_colors$best    <- list(fill = "brown3", border = "brown4",    font = "white",  arrow = "brown3") # Dropped (max)
  backward_colors$not_sig <- list(fill = "white",  border = "gray30",    font = "black",  arrow = "gray40") # Dropped
  backward_colors$sig     <- list(fill = "gray90", border = "gray30",    font = "grey30", arrow = "gray80") # Kept
}

plot_flowchart(
  backward_data,
  title         = "Backward Elimination - Warfarin PopPK",
  iter_label    = "Iteration",
  best_label    = "Dropped (best of iteration)",
  sig_label     = "Significant",
  base_label    = "Full model \n(all covariates)",
  colors        = backward_colors
)




# ----- Combined Flowchart
# Get a full forward + backward flowchart plot

# SCM
merged_data_SCM <- merge_scm_data(forward_data, backward_data, gap = 2)
# str(merged_data_SCM, max.level = 2) # Check 

plot_flowchart_combined(
  merged_data_SCM,
  title              = "Covariate Model Building - Warfarin PopPK",
  colors_backward    = backward_colors,
  node_size          = 30
)


# COSSAC + Colour management (optional)
merged_data_COSSAC <- merge_scm_data(cossac_data, backward_data, gap = 2)

plot_flowchart_combined(
  merged_data_COSSAC,
  title              = "Covariate Model Building - Warfarin PopPK",
  forward_title      = "Forward inclusion (COSSAC)",
  colors_forward     = modifyList(
    eval(formals(plot_flowchart_combined)$colors_forward),  
    list(best  = list(fill = "palegreen3", border = "palegreen4", font = "black",      arrow  = "palegreen4"))
  ), 
  colors_backward    = modifyList(
    eval(backward_colors),  
    list(base  = list(fill = "palegreen3", border = "palegreen4", font = "black",      arrow  = "palegreen4"))
  ), 
  node_size          = 30
)




# Final model                ===================================================

# Final model is logWT 
final_mod_path <- file.path(cur_dir, "warfarin_example/SCM/ModelBuilding_016/ModelBuilding_016.mlxtran")    # If SCM run
final_mod_path <- file.path(cur_dir, "warfarin_example/COSSAC/ModelBuilding_002/ModelBuilding_002.mlxtran") # If COSSAC run

# Let's run the full estimation and full report
loadProject(projectFile = final_mod_path)
saveProject(projectFile = "warfarin_final")

# Run missing estimations
getScenario() # SE and plots disabled
runStandardErrorEstimation()

computeChartsData()
saveProject()

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
run_comparison("warfarin_example", "warfarin_final", label_1 = "Base", label_2 = "Final")





# Model summary comparison   ===================================================

# Apart from the comparison template provided (see above), you might want to use
# the functions independently or compare more than two models at once.


  ## Two models comparison      ------------------------------------------------
    ### ----- Situation A -----
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

one_param_two_mod <- compare_two_models(
  first  = model_with, 
  second = model_without, 
  values = "estimate", 
  clean  = TRUE
  )

# A list: the data and a ready-made table
one_param_two_mod$data    # raw table
one_param_two_mod$kable   # styled table (the footnote gives the two models)

# The table is a kableExtra object, so it can be styled further
one_param_two_mod$kable |> 
  kable_styling(font_size = 9)


# You can also use
compare_criteria(first = model_with, second = model_without)
# (compare_obspred_metrics() compares the prediction errors of models whose ChartsData was exported:
#  see situation B)



    ### ----- Situation B -----
# What changed from the initial and final model ?

base  <- read_summary("warfarin_example/summary.txt", full = TRUE)
final <- read_summary("warfarin_final/summary.txt", full = TRUE)

# -> Now we look at more columns (and compute the differences)

#' ----- compare_two_models():
#' Fully join TWO models and compare one or several columns (delta columns computed)
#'     - first and second : models (read_summary(..., full = TRUE) objects, summary.txt paths or model folders)
#'     - values           : column(s) to compare (estimate, rse_pct etc.). Comparison following : second - first
#'     - clean            : if TRUE (default), clean covariates label (beta_CL_WT -> ~WT)

# bold_above: values above a threshold are printed in bold in the table (here RSE > 30 %)
all_param_two_mod <- compare_two_models(
  first      = base, 
  second     = final, 
  values     = c("estimate", "rse_pct"), 
  clean      = TRUE,
  bold_above = c(rse_pct = 30)
  )

all_param_two_mod$data
all_param_two_mod$kable

# Any other formatting goes through kableExtra, e.g. italic relative change of the estimates (column 5)
all_param_two_mod$kable |> 
  column_spec(5, italic = TRUE)

# Prediction errors (MPE, RMSE, ...) of the two models; read from their ChartsData
compare_obspred_metrics(first = base, second = final, labels = c("Base", "Final"))



  ## Multiple models comparison ------------------------------------------------

# Just for example : COSSAC workflow models comparison
# A real situation coul be : Is IOV useful? Try IOV implementation on each parameter and compare.

# One column
one_param_all_mod <- compare_models(dir = "warfarin_example/COSSAC", values = "estimate", clean = TRUE)
one_param_all_mod$data
one_param_all_mod$kable

# Multiple columns
all_param_all_mod <- compare_models(dir = "warfarin_example/COSSAC", values = c("estimate", "rse_pct"), clean = TRUE,
                                    bold_above = c(rse_pct = 30))
all_param_all_mod$data
all_param_all_mod$kable




# Convergence assessment     ===================================================

#' @If not done already, run: loadProject(projectFile = "warfarin_example.mlxtran")

# runAssessment()

# -> Personalized settings
set <- getAssessmentSettings()
set$nbRuns             <- 10
set$extendedEstimation <- TRUE
set$useLin             <- TRUE
set$initialParameters  <- data.frame(parameters = c("ka_pop", "V_pop", "Cl_pop"),
                                     min = c(0.1, 0.1, 0.01), max = c(10, 10, 1),
                                     fixed = FALSE)
runAssessment(settings = set)

AssessmentResults <- getAssessmentResults()

# -> Save the results to keep them with the model (they only exist in the R session otherwise)
saveRDS(AssessmentResults, "warfarin_example/AssessmentResults.rds")

# -> Plots: the lines of the runs should end on the same values
plot_convergence_assessment(AssessmentResults, which = "iterations")

# Final estimate of each run (ci = TRUE adds the 95 % interval)
plot_convergence_assessment(AssessmentResults, which = "estimates", ci = FALSE)

# The reshaped data, if you want to make your own plots
tidy_convergence_assessment(AssessmentResults)

# -> It can also be a section of the model report (see run_diagnostics()):
run_diagnostics("warfarin_example", assessment = "warfarin_example/AssessmentResults.rds")
