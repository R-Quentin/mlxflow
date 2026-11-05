# Functions: Structure & covariates  ===========================================
  ## 1. Read population parameters                              ----------------
#' Read the population parameters of a model
#'
#' Reads populationParameters.txt of a Monolix results folder.
#'
#' @param model_dir Monolix results folder of the model (the folder holding populationParameters.txt).
#'
#' @return tibble with one row per population parameter.
#'
#' @export
read_pop_params <- function(model_dir) {
  
  path  <- file.path(model_dir, "populationParameters.txt")
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  
  header <- strsplit(lines[1], ",")[[1]]
  n_cols <- length(header)
  
  rows <- map_dfr(lines[-1], function(ln) {
    parts        <- strsplit(ln, ",")[[1]]
    length(parts) <- n_cols          # pad with NA if fewer columns (fixed params)
    as.data.frame(t(parts), stringsAsFactors = FALSE)
  })
  
  colnames(rows) <- header
  
  rows |>
    mutate(across(-parameter, ~ suppressWarnings(as.numeric(.)))) |>
    as_tibble()
}





  ## 2. Parse .mlxtran                                          ----------------
    ### A) Resolve the .mlxtran path for a model folder                  -------
#' Standalone so other readers (e.g. read_obspred -> detect_dv_from_mlxtran)
#'
#' can resolve the file without going through the full [INDIVIDUAL]/
#' [COVARIATE] parsing done by parse_mlxtran(), which requires sections that
#' are irrelevant to DV detection and may be absent in some model exports.
#'
#' @noRd
resolve_mlxtran_path <- function(model_dir) {
  
  # Look for .mlxtran in the model folder, then one level up (project root).
  # If several .mlxtran files are found, prefer the one whose name matches
  # the model folder name (Monolix project folders can accumulate multiple
  # historical .mlxtran files, e.g. several model versions side by side).
  mlxtran_files <- list.files(model_dir, pattern = "\\.mlxtran$",
                              full.names = TRUE, recursive = FALSE)
  search_dir <- model_dir
  
  if (length(mlxtran_files) == 0) {
    search_dir    <- dirname(model_dir)
    mlxtran_files <- list.files(search_dir, pattern = "\\.mlxtran$",
                                full.names = TRUE, recursive = FALSE)
  }
  if (length(mlxtran_files) == 0) stop("No .mlxtran file found.")
  
  if (length(mlxtran_files) > 1) {
    
    model_folder_name <- basename(model_dir)
    name_match <- mlxtran_files[
      str_detect(basename(mlxtran_files), stringr::fixed(model_folder_name))
    ]
    
    if (length(name_match) >= 1) {
      mlxtran_path <- name_match[1]
    } else {
      mlxtran_path <- mlxtran_files[1]
    }
    
    warning(sprintf(
      "Several .mlxtran files found in '%s':\n  %s\nUsing: %s\nIf this is not the correct file, pass a model_dir whose .mlxtran is unambiguous.",
      search_dir,
      paste(basename(mlxtran_files), collapse = ", "),
      basename(mlxtran_path)
    ), call. = FALSE)
    
  } else {
    mlxtran_path <- mlxtran_files[1]
  }
  
  mlxtran_path
}




    ### B) Top-level parser                                              -------
#' parse_mlxtran
#'
#' @noRd
parse_mlxtran <- function(model_dir) {
  
  mlxtran_path <- resolve_mlxtran_path(model_dir)
  lines <- readLines(mlxtran_path, warn = FALSE)
  
  # B) [INDIVIDUAL] section -> DEFINITION lines --------------------------------
  indiv_start   <- which(str_detect(lines, "^\\[INDIVIDUAL\\]"))
  next_section  <- which(str_detect(lines, "^\\[") &
                           seq_along(lines) > indiv_start)[1]
  indiv_lines   <- lines[indiv_start:next_section]
  def_lines     <- indiv_lines[str_detect(indiv_lines, "distribution=")]
  
  param_defs <- map(def_lines, parse_individual_def) |>
    compact() |>
    bind_rows()
  
  # C) [COVARIATE] section -> continuous covariate transformations -------------
  #    and categorical transformations (DEFINITION: tXXX = {transform=XXX, ...})
  cov_start <- which(str_detect(lines, "^\\[COVARIATE\\]"))
  if (length(cov_start) > 0) {
    next_sec       <- which(str_detect(lines, "^\\[") &
                              seq_along(lines) > cov_start[1])[1]
    cov_lines      <- lines[cov_start[1]:next_sec]
    cov_transforms <- parse_covariate_transforms(cov_lines)
    cat_transforms <- parse_categorical_transforms(cov_lines)
  } else {
    cov_transforms <- tibble(cov_transformed     = character(),
                             cov_raw             = character(),
                             cov_center          = numeric(),
                             cov_scale           = numeric(),
                             cov_transform_type  = character())
    cat_transforms <- tibble(cov_transformed = character(),
                             cov_raw         = character(),
                             reference       = character())
  }
  
  list(param_defs     = param_defs,
       cov_transforms = cov_transforms,
       cat_transforms = cat_transforms,
       mlxtran_path   = mlxtran_path)
}




    ### C) Parse one INDIVIDUAL DEFINITION line                          -------
#' Handles the three coefficient syntaxes produced by Monolix:
#'
#'   1. scalar      : coefficient=beta_Cl_logtAGE
#'   2. plain vector: coefficient=\{beta_ka2_logtWT, beta_ka2_logtBMI_J0\}
#'   3. nested      : coefficient=\{\{0, beta_F1_tSEX_G_Female\}, beta_F1_logtWT, beta_F1_logtAGE\}
#'      (the inner \{\} holds the categorical dummy; remaining are continuous)
#'
#' The covariate vector always matches the coefficient vector positionally:
#'   covariate=\{tSEX, logtWT, logtAGE\}  <->  coefficient=\{\{0, beta_...\}, beta_..., beta_...\}
#' For the categorical slot, the inner \{ref_value, non_ref_beta\} is expanded
#' into two rows: one with coefficient "0" and one with the actual beta name.
#'
#' @noRd
parse_individual_def <- function(line) {
  
  line <- str_trim(line)
  if (!str_detect(line, "=\\s*\\{")) return(NULL)
  
  param_name <- str_extract(line, "^[A-Za-z0-9_]+")
  
  # Extract content between the outer braces
  content <- str_match(line, "\\{(.*)\\}")[2]
  if (is.na(content)) return(NULL)
  
  extract_field <- function(field) {
    m <- str_match(content, paste0(field, "\\s*=\\s*([A-Za-z0-9_]+)"))
    if (!is.na(m[1])) return(str_trim(m[2]))
    NA_character_
  }
  
  distribution <- extract_field("distribution")
  typical      <- extract_field("typical")
  no_var       <- str_detect(content, "no-variability")
  
  # --- Parse covariate vector (may be scalar or {a, b, c}) -------------------
  cov_match <- str_match(content, "covariate\\s*=\\s*\\{([^}]*)\\}")
  if (!is.na(cov_match[1])) {
    covariates <- str_split(cov_match[2], ",")[[1]] |> str_trim()
  } else {
    cov_scalar <- str_match(content, "covariate\\s*=\\s*([A-Za-z0-9_]+)")
    covariates <- if (!is.na(cov_scalar[1])) str_trim(cov_scalar[2]) else character(0)
  }
  
  # --- Parse coefficient vector (scalar / plain vector / nested vector) -------
  #
  # The coefficient vector can mix plain scalars and one categorical group
  # {ref_val, non_ref_beta} at ANY position -- not necessarily first. E.g.:
  #   coefficient={{0, beta_F1_tSEX_G_Female}, beta_F1_logtWT, beta_F1_logtAGE}
  #   coefficient={beta_F1_logtWT, {0, beta_F1_tSEX_G_Female}}
  # Both must resolve to the same covariate <-> coefficient pairing, driven
  # by position, not by assuming the categorical group comes first.
  #
  # tokenize_top_level() splits a comma-separated string while respecting one
  # level of {} nesting, so "{0, beta_X}, beta_Y" -> c("{0, beta_X}", "beta_Y").
  #
  # extract_braced_field() finds the value of a named field inside `content`
  # by walking character-by-character, tracking brace depth. This handles
  # fields like  coefficient={{...}, beta_X}, varlevel={id, id*occ}, sd={...}
  # correctly -- it stops at the matching closing brace, not at the greedy
  # last brace in the string. This is critical when varlevel= or other
  # brace-containing fields follow coefficient= in the same definition.
  tokenize_top_level <- function(s) {
    depth <- 0
    tokens <- character(0)
    current <- ""
    chars <- strsplit(s, "")[[1]]
    for (ch in chars) {
      if (ch == "{") depth <- depth + 1
      if (ch == "}") depth <- depth - 1
      if (ch == "," && depth == 0) {
        tokens <- c(tokens, current)
        current <- ""
      } else {
        current <- paste0(current, ch)
      }
    }
    tokens <- c(tokens, current)
    str_trim(tokens)
  }
  
  extract_braced_field <- function(content, field) {
    # Find "field = {" in content, then extract up to the matching "}"
    pat   <- paste0(field, "\\s*=\\s*\\{")
    start <- regexpr(pat, content, perl = TRUE)
    if (start[1] == -1) return(NA_character_)
    
    # Position of the opening "{" of the field value
    brace_open <- start[1] + attr(start, "match.length")[1] - 1
    chars <- strsplit(content, "")[[1]]
    depth <- 0
    for (i in brace_open:length(chars)) {
      if (chars[i] == "{") depth <- depth + 1
      if (chars[i] == "}") {
        depth <- depth - 1
        if (depth == 0) {
          # Return the content BETWEEN the outer braces
          return(paste(chars[(brace_open + 1):(i - 1)], collapse = ""))
        }
      }
    }
    NA_character_  # unbalanced braces -- should not happen in valid .mlxtran
  }
  
  # Try to extract coefficient as a braced field first; fall back to scalar
  coef_inner <- extract_braced_field(content, "coefficient")
  
  if (!is.na(coef_inner) && str_detect(content, "coefficient\\s*=\\s*\\{")) {
    
    coef_tokens <- tokenize_top_level(coef_inner)
    
    # A token like "{0, beta_F1_tSEX_G_Female}" is the categorical group;
    # a token like "beta_F1_logtWT" is a plain continuous coefficient.
    is_cat_token <- str_detect(coef_tokens, "^\\{")
    
    if (any(is_cat_token)) {
      
      cat_pos <- which(is_cat_token)[1]
      inner    <- str_match(coef_tokens[cat_pos], "^\\{(.*)\\}$")[2]
      inner_parts <- str_split(inner, ",")[[1]] |> str_trim()
      ref_val     <- inner_parts[1]
      nonref_beta <- inner_parts[2]
      
      # covariates and coef_tokens are positionally aligned (1:1, same order
      # as written in the .mlxtran), so the covariate at cat_pos is the
      # categorical one, regardless of where it falls in the sequence.
      rows <- map_dfr(seq_along(coef_tokens), function(j) {
        if (j == cat_pos) {
          tibble(parameter = param_name, distribution = distribution,
                 typical = typical, no_variability = no_var,
                 covariate = covariates[j],
                 coefficient = c(ref_val, nonref_beta))
        } else {
          tibble(parameter = param_name, distribution = distribution,
                 typical = typical, no_variability = no_var,
                 covariate = covariates[j],
                 coefficient = coef_tokens[j])
        }
      })
      
      return(rows)
    }
    
    # No nested {} categorical group found above. But there is a second,
    # flatter syntax Monolix uses when the categorical covariate is the
    # ONLY covariate on this parameter: covariate is then a plain scalar
    # (e.g. "tSEX", not "{tSEX, ...}") while coefficient is a flat
    # two-element vector "{0, beta_X}" -- no inner braces at all, because
    # there's nothing else to wrap it against. Without this check, the
    # code below would zip a length-1 covariate vector against a length-2
    # coefficient vector and silently drop the real beta behind an NA
    # covariate. Signature: exactly 2 coefficient tokens, the first being
    # the literal reference value "0", and covariate already scalar.
    if (length(coef_tokens) == 2 && coef_tokens[1] == "0" && length(covariates) == 1) {
      return(tibble(
        parameter    = param_name,
        distribution = distribution,
        typical      = typical,
        no_variability = no_var,
        covariate    = covariates[1],
        coefficient  = coef_tokens   # c("0", "beta_X")
      ))
    }
    
    # No categorical group: plain vector, positionally aligned
    coefficients <- coef_tokens
    n <- max(length(covariates), length(coefficients))
    length(covariates)   <- n
    length(coefficients) <- n
    
    return(tibble(
      parameter    = param_name,
      distribution = distribution,
      typical      = typical,
      no_variability = no_var,
      covariate    = covariates,
      coefficient  = coefficients
    ))
  }
  
  # Case 1: scalar coefficient (single covariate, no braces)
  {
    coef_scalar  <- str_match(content, "coefficient\\s*=\\s*([A-Za-z0-9_]+)")
    coefficients <- if (!is.na(coef_scalar[1])) str_trim(coef_scalar[2]) else character(0)
  }
  
  # Ensure at least one row per parameter even when it has no covariate at
  # all (e.g. "no-variability" parameters like D3, ka1, V3...), so that
  # build_formula_table() can still list every PK parameter, not just the
  # ones with an active covariate effect.
  n <- max(length(covariates), length(coefficients), 1L)
  length(covariates)   <- n
  length(coefficients) <- n
  
  tibble(
    parameter    = param_name,
    distribution = distribution,
    typical      = typical,
    no_variability = no_var,
    covariate    = covariates,
    coefficient  = coefficients
  )
}




    ### D) Parse centred covariate transformations from [COVARIATE]      -------
#' Two forms produced by Monolix, both centred continuous covariates:
#'
#'   log-centred:    logtAGE = log(AGE / 44)
#'   linear-centred: FBMI    = (BMI_J0 - 24.5) / 24.5
#'
#' @return one combined tibble with a "cov_transform_type" column ("log" or "linear") so downstream code (compute_typical(), compute_covariate_effects()) can apply the correct formula instead of assuming log-centring for every continuous covariate.
#'
#' @noRd
parse_covariate_transforms <- function(cov_lines) {
  
  log_lines <- cov_lines[str_detect(cov_lines, "=\\s*log\\(")]
  log_hits <- map_dfr(log_lines, function(ln) {
    m <- str_match(str_trim(ln),
                   "^(\\w+)\\s*=\\s*log\\(([A-Za-z0-9_]+)\\s*/\\s*([0-9.]+)\\)")
    if (is.na(m[1])) return(NULL)
    tibble(cov_transformed     = m[2],
           cov_raw             = m[3],
           cov_center          = as.numeric(m[4]),
           cov_scale           = NA_real_,
           cov_transform_type  = "log")
  })
  
  lin_lines <- cov_lines[str_detect(
    cov_lines, "=\\s*\\(\\s*[A-Za-z0-9_]+\\s*-\\s*[0-9.]+\\s*\\)\\s*/\\s*[0-9.]+"
  )]
  lin_hits <- map_dfr(lin_lines, function(ln) {
    m <- str_match(str_trim(ln),
                   "^(\\w+)\\s*=\\s*\\(\\s*([A-Za-z0-9_]+)\\s*-\\s*([0-9.]+)\\s*\\)\\s*/\\s*([0-9.]+)")
    if (is.na(m[1])) return(NULL)
    tibble(cov_transformed     = m[2],
           cov_raw             = m[3],
           cov_center          = as.numeric(m[4]),
           cov_scale           = as.numeric(m[5]),
           cov_transform_type  = "linear")
  })
  
  bind_rows(log_hits, lin_hits)
}




    ### E) Parse categorical transformations from [COVARIATE] DEFINITION -------
#' e.g
#'
#'   tSEX = \{
#'     transform = SEX,
#'     categories = \{ 'G_Male' = \{'Male'\}, 'G_Female' = \{'Female'\} \},
#'     reference = 'G_Male'
#'   \}
#'
#' @return tibble mapping the transformed name (tSEX) back to the raw covariate name (SEX) and its reference category (G_Male). This lets downstream code strip the "t" prefix correctly -- only for names that are genuinely categorical transforms, never a naive substring removal (so a covariate literally named "Tabac" is untouched unless it has its own "tTabac = \{transform = Tabac, ...\}" block).
#'
#' @noRd
parse_categorical_transforms <- function(cov_lines) {
  
  def_start <- which(str_detect(cov_lines, "^DEFINITION:"))
  if (length(def_start) == 0) {
    return(tibble(cov_transformed = character(),
                  cov_raw         = character(),
                  reference       = character()))
  }
  
  block <- paste(cov_lines[def_start:length(cov_lines)], collapse = " ")
  
  # Each transform definition: name = { transform = RAW, ... reference = 'REF' }
  matches <- str_match_all(
    block,
    "(\\w+)\\s*=\\s*\\{\\s*transform\\s*=\\s*(\\w+).*?reference\\s*=\\s*'([^']+)'\\s*\\}"
  )[[1]]
  
  if (nrow(matches) == 0) {
    return(tibble(cov_transformed = character(),
                  cov_raw         = character(),
                  reference       = character()))
  }
  
  tibble(cov_transformed = matches[, 2],
         cov_raw         = matches[, 3],
         reference       = matches[, 4])
}
#' Strips a Monolix-generated "t" prefix from a covariate name, but only when
#'
#' it matches a declared categorical transform (cat_transforms) -- never a
#' blind substring removal. Falls back to the original name otherwise.
#'   strip_categorical_prefix("tSEX",   cat_transforms) -> "SEX"
#'   strip_categorical_prefix("Tabac",  cat_transforms) -> "Tabac"  (untouched)
#'
#' @noRd
strip_categorical_prefix <- function(cov_name, cat_transforms) {
  hit <- cat_transforms |> filter(cov_transformed == cov_name)
  if (nrow(hit) > 0) hit$cov_raw[1] else cov_name
}





  ## 3. Read individual parameters (for observed percentiles)   ----------------
#' Read the estimated individual parameters of a model
#'
#' Reads IndividualParameters/estimatedIndividualParameters.txt and keeps one row per individual.
#'
#' @param model_dir Monolix results folder of the model.
#'
#' @return tibble with one row per individual.
#'
#' @export
read_indiv_params <- function(model_dir) {
  path <- file.path(model_dir,
                    "IndividualParameters/estimatedIndividualParameters.txt")
  read_csv(path, show_col_types = FALSE) |>
    group_by(id) |>
    slice(1) |>
    ungroup()
}





  ## 4. Compute typical parameter value along a covariate range ----------------
#' All other covariates are held at their reference value
#'
#' For centred covariates (log- or linear-centred), the reference is the
#' centre (transformed value = 0). For categorical indicators the reference
#' is coded 0.
#'
#' @noRd
compute_typical <- function(param_name, cov_name, cov_values,
                            param_defs, pop_params, cov_transforms) {
  
  def <- param_defs |> filter(parameter == param_name)
  if (nrow(def) == 0) return(NULL)
  
  distribution <- def$distribution[1]
  typical_val  <- pop_params |>
    filter(parameter == def$typical[1]) |>
    pull(value)
  if (length(typical_val) == 0 || is.na(typical_val)) return(NULL)
  
  # Coefficient for the covariate of interest
  beta_row <- def |> filter(covariate == cov_name, coefficient != "0")
  if (nrow(beta_row) == 0) return(NULL)
  
  beta_val <- pop_params |>
    filter(parameter == beta_row$coefficient[1]) |>
    pull(value)
  if (length(beta_val) == 0 || is.na(beta_val)) return(NULL)
  
  # Transform raw covariate values to model space
  tf <- cov_transforms |> filter(cov_transformed == cov_name)
  if (nrow(tf) > 0 && tf$cov_transform_type[1] == "log") {
    cov_model <- log(cov_values / tf$cov_center[1])              # log-centred
  } else if (nrow(tf) > 0 && tf$cov_transform_type[1] == "linear") {
    cov_model <- (cov_values - tf$cov_center[1]) / tf$cov_scale[1]  # linear-centred
  } else {
    cov_model <- as.numeric(cov_values)               # categorical 0 / 1
  }
  
  # Apply PK formula (other covariates contribute 0 at their reference)
  if (distribution == "logNormal") {
    param_values <- typical_val * exp(beta_val * cov_model)
  } else if (distribution == "logitNormal") {
    logit_ref    <- log(typical_val / (1 - typical_val))
    logit_p      <- logit_ref + beta_val * cov_model
    param_values <- exp(logit_p) / (1 + exp(logit_p))
  } else {
    param_values <- typical_val + beta_val * cov_model
  }
  
  tibble(parameter     = param_name,
         covariate     = cov_name,
         cov_value     = cov_values,
         param_typical = param_values)
}





  ## 5. Build formula table                                     ----------------
#' Reconstructs, for every PK parameter defined in [INDIVIDUAL], the formula
#'
#' shown by Monolix (e.g. "log(Cl) = log(Cl_pop) + beta_Cl_logtAGE*logtAGE"),
#' using the link function implied by the distribution and resolving raw
#' covariate / parameter names where possible.
#'
#' @return tibble: [parameter, distribution, formula] - one row per PK parameter (covariates of a given parameter are combined into a single formula string, matching the Monolix interface display)
#'
#' @noRd
build_formula_table <- function(param_defs, cov_transforms, cat_transforms) {
  
  link_fun <- function(distribution, param_name) {
    switch(distribution,
           logNormal   = sprintf("log(%s)", param_name),
           logitNormal = sprintf("log( %s / ( 1 - %s ) )", param_name, param_name),
           param_name)
  }
  
  params <- param_defs |> distinct(parameter, distribution, typical, no_variability)
  
  map_dfr(seq_len(nrow(params)), function(i) {
    
    pk_param     <- params$parameter[i]
    distribution <- params$distribution[i]
    typical_name <- params$typical[i]
    no_var       <- params$no_variability[i]
    
    lhs <- link_fun(distribution, pk_param)
    rhs_typical <- link_fun(distribution, typical_name)
    
    terms <- param_defs |>
      filter(parameter == pk_param, !is.na(covariate), !is.na(coefficient))
    
    cov_terms <- if (nrow(terms) > 0) {
      terms_active <- terms |> filter(coefficient != "0")
      if (nrow(terms_active) == 0) {
        ""
      } else {
        paste0(" + ", paste(
          sprintf("%s*%s", terms_active$coefficient, terms_active$covariate),
          collapse = " + "
        ))
      }
    } else {
      ""
    }
    
    eta_term <- if (isTRUE(no_var)) "" else sprintf(" + eta_ID_%s", pk_param)
    
    formula <- sprintf("%s = %s%s%s", lhs, rhs_typical, cov_terms, eta_term)
    
    tibble(parameter    = pk_param,
           distribution = distribution,
           formula      = formula)
  })
}





  ## 6. Main function: compute_covariate_effects()                ----------------
#' Covariate effects analysis for a PopPK Monolix model
#'
#' @param model_dir   Path to the Monolix model output folder.
#' @param n_points    Number of points for continuous covariate curves (default 200).
#' @param ref_color   Colour used for the reference lines / bars.
#' @param curve_color Colour used for the effect curve / bars.
#'
#' @return A named list:
#'   $table  tibble with relative effects at P5 / P50 / P95 (or category levels)
#'   $plots  named list of patchwork panels, one per PK parameter
#'   $data   named list of raw curve tibbles (one per parameter x covariate pair)
#' @export
compute_covariate_effects <- function(model_dir,
                                    n_points    = 200,
                                    ref_color   = "grey60",
                                    curve_color = "dodgerblue3") {
  
  ## 5A. Read files ------------------------------------------------------------
  message("Reading Monolix files...")
  
  # Z) summary.txt
  summary <- read_summary(file.path(model_dir, "summary.txt"), full = TRUE)
  # model_name <- summary$model_name
  message("  [ok] summary.txt  (Model: ", summary$model_name, ")")
  
  
  # A) Population parameters
  pop_params <- read_pop_params(model_dir)
  message("  [ok] populationParameters.txt  (", nrow(pop_params), " parameters)")
  
  # B) .mlxtran
  mlx            <- parse_mlxtran(model_dir)
  param_defs     <- mlx$param_defs
  cov_transforms <- mlx$cov_transforms
  cat_transforms <- mlx$cat_transforms
  mlxtran_path   <- mlx$mlxtran_path
  message("  [ok] .mlxtran parsed  (",
          n_distinct(param_defs$parameter), " PK parameters defined)")
  message("       file used: ", basename(mlxtran_path))
  
  formula_table <- build_formula_table(param_defs, cov_transforms, cat_transforms)
  message("  [ok] Formulas reconstructed for ", nrow(formula_table), " PK parameters")
  
  # C) Individual parameters 
  indiv_data <- read_indiv_params(model_dir)
  message("  [ok] Individual parameters  (", nrow(indiv_data), " subjects)")
  
  # D) Active covariate-parameter pairs (exclude reference category rows)
  active_pairs <- param_defs |>
    filter(!is.na(covariate), !is.na(coefficient), coefficient != "0") |>
    select(parameter, covariate, coefficient) |>
    distinct()
  
  # Reorder parameters and covariates to match the model's own order, as
  # printed in summary.txt / populationParameters.txt, rather than the
  # .mlxtran DEFINITION order or an alphabetical sort. This makes the table
  # and combined plots easier to scan against the summary side by side.
  #
  # Parameter order: taken from the "*_pop" rows of pop_params (one row per
  # PK parameter, e.g. "D3_pop", "F3_pop", "Cl_pop", ...), in file order.
  # Covariate order, within a parameter: each active pair already carries
  # its exact coefficient name (e.g. "beta_F1_tSEX_G_Female"), so its rank
  # is just its position among the "beta_" rows of pop_params -- no need to
  # reconstruct the name from parameter/covariate.
  param_order <- pop_params$parameter |>
    keep(~ str_detect(.x, "_pop$")) |>
    str_remove("_pop$")
  
  if (length(param_order) > 0) {
    
    cov_order_lookup <- pop_params$parameter |>
      keep(~ str_detect(.x, "^beta_"))
    
    active_pairs <- active_pairs |>
      mutate(
        .param_rank = match(parameter, param_order),
        .param_rank = if_else(is.na(.param_rank), Inf, as.numeric(.param_rank)),
        .cov_rank   = {
          r <- match(coefficient, cov_order_lookup)
          if_else(is.na(r), Inf, as.numeric(r))
        }
      ) |>
      arrange(.param_rank, .cov_rank) |>
      select(-.param_rank, -.cov_rank)
  }
  
  message("\n", nrow(active_pairs), " active covariate effect(s) detected:")
  walk(seq_len(nrow(active_pairs)), function(i)
    message("   * ", active_pairs$covariate[i],
            " -> ", active_pairs$parameter[i],
            "  (", active_pairs$coefficient[i], ")")
  )
  
  
  ## 5B. Compute curves --------------------------------------------------------
  all_curves <- list()
  
  for (i in seq_len(nrow(active_pairs))) {
    
    pk_param  <- active_pairs$parameter[i]
    cov_name  <- active_pairs$covariate[i]
    tf        <- cov_transforms |> filter(cov_transformed == cov_name)
    
    if (nrow(tf) > 0) {
      
      # ---- Continuous covariate (log- or linear-centred)
      raw_name <- tf$cov_raw[1]
      center   <- tf$cov_center[1]
      
      if (raw_name %in% colnames(indiv_data)) {
        raw_vals <- indiv_data[[raw_name]]
        p5  <- quantile(raw_vals, 0.05, na.rm = TRUE)
        p50 <- quantile(raw_vals, 0.50, na.rm = TRUE)
        p95 <- quantile(raw_vals, 0.95, na.rm = TRUE)
      } else {
        p5  <- center * 0.6; p50 <- center; p95 <- center * 1.6
      }
      
      cov_seq <- seq(p5, p95, length.out = n_points)
      curve   <- compute_typical(pk_param, cov_name, cov_seq,
                                 param_defs, pop_params, cov_transforms)
      if (is.null(curve)) next
      
      curve$cov_raw  <- raw_name
      curve$cov_center <- center
      curve$p5  <- p5; curve$p50 <- p50; curve$p95 <- p95
      curve$cov_type <- "continuous"
      
    } else {
      
      # ---- Categorical covariate
      beta_name   <- active_pairs$coefficient[i]
      beta_val    <- pop_params |> filter(parameter == beta_name) |> pull(value)
      def         <- param_defs |> filter(parameter == pk_param)
      distribution <- def$distribution[1]
      typical_val <- pop_params |>
        filter(parameter == def$typical[1]) |> pull(value)
      
      if (length(typical_val) == 0 || is.na(typical_val) ||
          length(beta_val)    == 0 || is.na(beta_val))     next
      
      cat_values <- c(0, 1)
      cat_label  <- str_remove(beta_name,
                               paste0("beta_", pk_param, "_", cov_name, "_"))
      raw_name   <- strip_categorical_prefix(cov_name, cat_transforms)
      
      if (distribution == "logNormal") {
        param_vals <- typical_val * exp(beta_val * cat_values)
      } else if (distribution == "logitNormal") {
        logit_ref  <- log(typical_val / (1 - typical_val))
        param_vals <- exp(logit_ref + beta_val * cat_values) /
          (1 + exp(logit_ref + beta_val * cat_values))
      } else {
        param_vals <- typical_val + beta_val * cat_values
      }
      
      curve <- tibble(
        parameter     = pk_param,
        covariate     = cov_name,
        cov_raw       = raw_name,
        cov_value     = cat_values,
        cov_label     = c("Reference", cat_label),
        param_typical = param_vals,
        cov_type      = "categorical"
      )
    }
    
    all_curves[[paste(pk_param, cov_name, sep = "_")]] <- curve
  }
  
  
  ## 5C. Summary table ---------------------------------------------------------
  message("\nBuilding summary table...")
  
  summary_table <- map_dfr(all_curves, function(curve) {
    
    if (curve$cov_type[1] == "continuous") {
      
      ref_y  <- curve$param_typical[which.min(abs(curve$cov_value - curve$p50[1]))]
      p5_y   <- curve$param_typical[which.min(abs(curve$cov_value - curve$p5[1]))]
      p95_y  <- curve$param_typical[which.min(abs(curve$cov_value - curve$p95[1]))]
      
      tibble(
        Parameter           = curve$parameter[1],
        Covariate           = curve$cov_raw[1],
        Relationship         = sprintf("%s ~ %s", curve$parameter[1], curve$cov_raw[1]),
        Type                = "Continuous",
        Level               = c(sprintf("P5  (%.1f)",      curve$p5[1]),
                                sprintf("P50 (%.1f) [ref]", curve$p50[1]),
                                sprintf("P95 (%.1f)",      curve$p95[1])),
        Param_value         = c(p5_y, ref_y, p95_y),
        Relative_change_pct = round((c(p5_y, ref_y, p95_y) / ref_y - 1) * 100, 1)
      )
      
    } else {
      ref_y <- curve$param_typical[curve$cov_value == 0]
      tibble(
        Parameter           = curve$parameter[1],
        Covariate           = curve$cov_raw[1],
        Relationship         = sprintf("%s ~ %s", curve$parameter[1], curve$cov_raw[1]),
        Type                = "Categorical",
        Level               = curve$cov_label,
        Param_value         = curve$param_typical,
        Relative_change_pct = round((curve$param_typical / ref_y - 1) * 100, 1)
      )
    }
  }) |>
    mutate(Relationship = factor(Relationship, levels = unique(Relationship))) |>
    arrange(Relationship) |>
    mutate(Relationship = as.character(Relationship))
  
  
  ## 5C. Plots -----------------------------------------------------------------
  message("Generating plots...")
  
  pk_params <- unique(map_chr(all_curves, ~ .$parameter[1]))
  
  plots           <- list()
  sub_plots_flat  <- list()   # every individual P~COV plot, flat (for combine_all_plots)
  
  for (pk_param in pk_params) {
    
    curves_param  <- keep(all_curves, ~ .$parameter[1] == pk_param)
    n_cols_param  <- min(3L, length(curves_param))
    
    # The physical width allotted to a panel shrinks as more panels share
    # the same row. Annotation text size stays fixed (in points), so a
    # margin calibrated for a single, full-width panel would clip once
    # placed in a 2- or 3-column grid. Each individual plot is therefore
    # built with a margin generous enough for the worst case (3 columns,
    # the max used anywhere in this module -- see n_cols_param above and
    # combine_all_plots' default ncol_panels). This keeps every individual
    # plot visually consistent and safe to reuse in any grid width up to 3
    # columns, whether shown alone, grouped by parameter, or in the
    # flat combined figure.
    margin_scale <- 1 + 0.45 * (3L - 1)
    
    sub_plots <- map(curves_param, function(curve) {
      
      if (curve$cov_type[1] == "continuous") {
        
        # ---- Continuous: line plot
        ref_x  <- curve$p50[1]
        ref_y  <- curve$param_typical[which.min(abs(curve$cov_value - ref_x))]
        p5_x   <- curve$p5[1]
        p5_y   <- curve$param_typical[which.min(abs(curve$cov_value - p5_x))]
        p95_x  <- curve$p95[1]
        p95_y  <- curve$param_typical[which.min(abs(curve$cov_value - p95_x))]
        rel_p5  <- round((p5_y  / ref_y - 1) * 100, 1)
        rel_p95 <- round((p95_y / ref_y - 1) * 100, 1)
        
        # Expand both axes so endpoint labels are never clipped, without
        # relying on coord_cartesian(clip = "off"). The margin is computed
        # from the actual data range of the curve (not just the P5-P95
        # window), since the curve itself defines panel limits in ggplot,
        # and scaled up when several panels share a row (margin_scale).
        x_span   <- diff(range(curve$cov_value))
        y_span   <- diff(range(curve$param_typical))
        x_margin <- x_span * 0.15 * margin_scale
        y_margin <- y_span * 0.18 * margin_scale
        
        # y position for the covariate reference label: just above the
        # bottom of the (expanded) panel, so it reads like an axis tick
        y_floor <- min(curve$param_typical) - y_margin * 0.55
        
        # x position for the parameter reference label: near the left edge
        # of the (expanded) panel, so it reads like a y-axis tick rather
        # than floating mid-plot
        x_left <- min(curve$cov_value) - x_margin * 0.85
        
        p <- ggplot(curve, aes(x = cov_value, y = param_typical)) +
          # Shaded P5-P95 band
          annotate("rect",
                   xmin = p5_x, xmax = p95_x, ymin = -Inf, ymax = Inf,
                   fill = "#DDEEFF", alpha = 0.6) +
          # Reference lines
          geom_vline(xintercept = ref_x, color = ref_color,
                     linetype = "dashed", linewidth = 0.55) +
          geom_hline(yintercept = ref_y, color = ref_color,
                     linetype = "dashed", linewidth = 0.55) +
          # Effect curve
          geom_line(color = curve_color, linewidth = 1.1) +
          # P5 / P95 point annotations
          # annotate() rather than geom_point(aes(x = <scalar>, ...)): the
          # latter inherits the 200-row curve data and draws 200 identical
          # points (and warns "All aesthetics have length 1, but the data
          # has 200 rows" at print time).
          annotate("point", x = p5_x,  y = p5_y,  color = curve_color, size = 2.8) +
          annotate("point", x = p95_x, y = p95_y, color = curve_color, size = 2.8) +
          annotate("text", x = p5_x,  y = p5_y,
                   label  = sprintf("%+.0f%%", round(rel_p5) + 0),
                   hjust  = 1.15, vjust = -0.5, size = 3.2,
                   color  = curve_color, fontface = "bold") +
          annotate("text", x = p95_x, y = p95_y,
                   label  = sprintf("%+.0f%%", round(rel_p95) + 0),
                   hjust  = -0.15, vjust = -0.5, size = 3.2,
                   color  = curve_color, fontface = "bold") +
          # Reference label for the parameter (y axis side) -- placed near
          # the left edge of the panel rather than at ref_x, so it reads
          # like an axis annotation instead of floating mid-curve
          annotate("text", x = x_left, y = ref_y,
                   label  = sprintf("ref: %.3g", ref_y),
                   hjust  = 0, vjust = -0.9, size = 2.8,
                   color  = ref_color) +
          # Reference label for the covariate (x axis side)
          annotate("text", x = ref_x, y = y_floor,
                   label  = sprintf("ref: %.3g", ref_x),
                   hjust  = -0.08, vjust = 0, size = 2.8,
                   color  = ref_color) +
          # Expanded limits on both axes to prevent label clipping
          scale_x_continuous(expand = expansion(add = x_margin)) +
          scale_y_continuous(expand = expansion(add = y_margin)) +
          labs(title = sprintf("%s ~ %s", pk_param, curve$cov_raw[1]),
               x     = curve$cov_raw[1],
               y     = pk_param) +
          theme_bw(base_size = 11) +
          theme(plot.title      = element_text(face = "bold", size = 11),
                panel.grid.minor = element_blank(),
                plot.margin     = margin(8, 12, 8, 8))
        
      } else {
        
        # ---- Categorical: bar plot
        ref_y      <- curve$param_typical[curve$cov_value == 0]
        # Bars are always plotted from a 0 baseline (geom_col), so the
        # headroom must be sized relative to the full 0-to-max axis span,
        # not just the gap between the reference and effect values -- using
        # diff(range(param_typical)) alone underestimates the margin whenever
        # that gap is small relative to the absolute values (e.g. Cl ~ Tabac:
        # 3.17 -> 3.88 is only a modest jump far from 0, so its own range
        # barely reflects the true 0-3.88 panel scale), which clipped the
        # top label.
        y_span     <- diff(range(c(0, curve$param_typical)))
        y_margin   <- y_span * 0.25 * margin_scale   # headroom for value labels
        
        curve_plot <- curve |>
          mutate(rel_change = (param_typical / ref_y - 1) * 100,
                 cov_label  = factor(cov_label, levels = cov_label),
                 fill_col   = ifelse(cov_value == 0, ref_color, curve_color))
        
        p <- ggplot(curve_plot, aes(x = cov_label, y = param_typical,
                                    fill = fill_col)) +
          geom_col(width = 0.45) +
          geom_hline(yintercept = ref_y, color = ref_color,
                     linetype = "dashed", linewidth = 0.55) +
          geom_text(aes(label = ifelse(
            cov_value == 0,
            sprintf("ref\n%.3g", param_typical),
            sprintf("%+.0f%%\n(%.3g)", round(rel_change) + 0, param_typical))),
            vjust = -0.3, size = 3.2,
            color = "grey20", fontface = "bold") +
          scale_fill_identity() +
          scale_y_continuous(expand = expansion(mult = c(0.05, 0), add = c(0, y_margin))) +
          labs(title    = sprintf("%s ~ %s", pk_param, curve$cov_raw[1]),
               x        = curve$cov_raw[1],
               y        = pk_param) +
          theme_bw(base_size = 11) +
          theme(plot.title       = element_text(face = "bold", size = 11),
                panel.grid.minor = element_blank(),
                axis.text.x      = element_text(angle = 20, hjust = 1),
                plot.margin      = margin(8, 12, 8, 8))
      }
      
      sub_plots_flat[[paste(pk_param, curve$cov_raw[1], sep = "_")]] <<- p
      p
    })
    
    # Assemble all sub-plots for this PK parameter with patchwork
    plots[[pk_param]] <- wrap_plots(sub_plots, ncol = n_cols_param) +
      plot_annotation(
        title = sprintf("Covariate effects on %s  \u2014  %s", pk_param, summary$model_name),
        theme = theme(plot.title = element_text(face = "bold", size = 13))
      )
  }
  
  plots <- compact(plots)
  
  message("\nDone.  Access results via $model_name, $formulas, $table, $plots, $plots_flat, $data")
  
  list(model_name   = summary$model_name,
       mlxtran_path = mlxtran_path,
       formulas     = formula_table,
       table        = summary_table,
       plots        = plots,
       plots_flat   = sub_plots_flat,
       data         = all_curves,
       summary      = summary)
  
  
  
}





  ## 7. Combine all plots into a single patchwork figure        ----------------
#' By default, combines every individual P~COV plot flat (results$plots_flat),
#'
#' each built with margins calibrated for the panel width it will actually
#' get in this layout. This avoids the cramped annotations seen when
#' combining whole per-parameter blocks (which mix panels with 1, 2, or 3
#' columns into the same grid).
#'
#' Set by_parameter = TRUE to instead stack the existing per-parameter blocks
#' (results$plots) -- faster, but inherits each block's own internal layout.
#'
#' Usage:
#'   combined <- combine_all_plots(results)
#'   print(combined)
#'   ggsave("all_covariate_effects.png", combined, width = 14, height = 10)
#'
#' @param results Output of compute_covariate_effects().
#' @param ncol_panels Number of panel columns of the combined figure (default 3).
#' @param by_parameter If TRUE, one panel per PK parameter instead of one panel per covariate effect.
#'
#' @export
combine_all_plots <- function(results, ncol_panels = 3, by_parameter = FALSE) {
  
  if (by_parameter) {
    wrap_plots(results$plots, ncol = min(2L, length(results$plots))) +
      plot_annotation(
        title = sprintf("Covariate effects \u2014 all PK parameters  (%s)",
                        results$model_name),
        theme = theme(plot.title = element_text(face = "bold", size = 14))
      )
  } else {
    wrap_plots(results$plots_flat, ncol = ncol_panels) +
      plot_annotation(
        title = sprintf("Covariate effects \u2014 all PK parameters  (%s)",
                        results$model_name),
        theme = theme(plot.title = element_text(face = "bold", size = 14))
      )
  }
}
