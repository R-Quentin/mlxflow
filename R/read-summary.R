  ## read_summary() ------------------------------------------------------------
    ### ---- MODULE 1a: read_summary() ----
#' ---- INTERNAL HELPERS
#'
#' @return index of the first line matching a pattern, NA if absent
#'
#' @noRd
find_line <- function(lines, pat) {
  idx <- which(str_detect(lines, pat))
  if (length(idx) == 0) NA_integer_ else idx[1]
}
#' Removes lines that are blank or separator-only (---, ___, ***)
#'
#' @noRd
clean_block <- function(blk) {
  blk[!str_detect(blk, "^\\s*$|^[-_*]{3,}")]
}
#' Strips the _pop suffix from population parameter names for join with
#'
#' individual parameter names (e.g. "Cl_pop" -> "Cl")
#'
#' @noRd
strip_pop <- function(x) str_remove(x, "_pop$")
#' Extracts numeric values from the portion of a line that follows the parameter
#'
#' name, using a greedy float tokenizer.
#'
#' This avoids two pitfalls:
#'   - digits embedded in parameter names (e.g. "ka1_pop", "beta_F3_logtBMI_J0")
#'     being counted as values
#'   - pairs of numbers printed without a separating space by Monolix when values
#'     are small (e.g. "0.0989280.00020671" or "-0.0936450.00056686"), which a
#'     split-based approach mis-parses
#'
#' @noRd
nums_after_name <- function(line) {
  line <- str_trim(line)
  rest <- str_remove(line, "^[A-Za-z_][A-Za-z0-9_]*\\s*")
  as.numeric(str_extract_all(
    rest,
    "-?(?:[0-9]+\\.?[0-9]*|[0-9]*\\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
  )[[1]])
}
#' Parses a Monolix summary.txt and extracts:
#'
#'   - Model name
#'   - Information criteria (OFV, AIC, BICc, BIC)
#'   - Active covariate relations from beta_ fixed effects
#'   - Dataset information
#'
#' When full = TRUE, also extracts:
#'   - Population parameters with SE, RSE\%, and 90\% CI when available
#'   - Individual parameters (Conditional Distribution and/or Conditional Mode)
#'     with shrinkage percentages
#'   - Fisher correlation matrix (full symmetric) and eigen values
#'
#' @param file path to the summary.txt file
#' @param full if TRUE, parse all estimation sections (slower for batch use)
#'
#' @return named list: $model_name    character(1) $criteria      tibble  [criterion, value] $ll_method     "is" (importance sampling), "lin" (linearization), ... or NA $covariates    character vector of active covariate relations $dataset_info  tibble  [info, value] $pop_params    tibble or NULL $indiv_params  tibble or NULL $fisher_matrix data.frame (symmetric) or NULL $eigen_values  tibble or NULL
#'
#' @export
read_summary <- function(file, full = FALSE) {
  
  lines <- readLines(file, warn = FALSE)
  
  # --- Model name ---
  model_name <- lines[str_detect(lines, "\\.mlxtran")] |>
    str_trim() |>
    str_extract("[A-Za-z0-9_]+(?=\\.mlxtran)") |>
    na.omit() |>
    first()
  
  
  # --- Information criteria ---
  criteria_pats <- c(
    OFV  = "^-2 x log-likelihood\\s+\\(OFV\\)",
    AIC  = "^Akaike Information Criteria\\s+\\(AIC\\)",
    BICc = "^Corrected Bayesian Information Criteria\\s+\\(BICc\\)",
    BIC  = "^Bayesian Information Criteria\\s+\\(BIC\\)"
  )
  
  criteria <- imap_dfr(criteria_pats, function(pat, nm) {
    val <- lines[str_detect(lines, pat)] |>
      str_extract("[0-9]+\\.[0-9]+") |>
      as.numeric()
    tibble(criterion = nm, value = val)
  })
  
  
  # --- Log-likelihood method: the "(is)" / "(lin)" line printed just above the OFV ---
  # OFV, AIC and BIC are only comparable between models computed the same way
  ofv_line  <- which(str_detect(lines, criteria_pats[["OFV"]]))[1]
  ll_method <- NA_character_
  if (!is.na(ofv_line) && ofv_line > 1) {
    above <- str_trim(lines[seq_len(ofv_line - 1)])
    above <- above[nzchar(above)]
    if (length(above) > 0) ll_method <- str_match(utils::tail(above, 1), "^\\((\\w+)\\)$")[, 2]
  }
  
  
  # --- Active covariates ---
  # Only the population-parameter block is read ("beta_xxx   :   value" lines).
  # Searching the whole file also matched the first column of the Fisher
  # correlation matrix whenever a beta has the widest label (no leading
  # spaces): duplicated or glued names such as "V_logtWt-0.0089253".
  pop_first <- find_line(lines, "ESTIMATION OF THE POPULATION PARAMETERS")
  pop_last  <- coalesce(find_line(lines, "ESTIMATION OF THE INDIVIDUAL PARAMETERS"),
                        find_line(lines, "ESTIMATION OF THE FISHER INFORMATION MATRIX"),
                        find_line(lines, "ESTIMATION OF THE LOG-LIKELIHOOD"),
                        length(lines) + 1L)
  pop_lines <- if (is.na(pop_first)) lines else lines[seq(pop_first + 1L, max(pop_first + 1L, pop_last - 1L))]
  
  covariates <- pop_lines |>
    str_trim() |>
    str_subset("^beta_[^\\s:]+\\s*:") |>
    str_extract("^beta_[^\\s:]+") |>
    str_remove("^beta_")
  
  
  # --- Dataset information ---
  ds_start <- find_line(lines, "^DATASET INFORMATION")
  
  dataset_info <- if (!is.na(ds_start)) {
    ds_lines <- lines[(ds_start + 1):length(lines)]
    ds_lines <- ds_lines[str_detect(ds_lines, ":") & str_detect(ds_lines, "\\d")]
    ds_lines <- ds_lines[!str_detect(ds_lines, "^\\s*$")]
    map_dfr(ds_lines, function(ln) {
      parts <- str_split_fixed(str_trim(ln), "\\s*:\\s*", n = 2)
      tibble(info  = str_trim(parts[1]),
             value = suppressWarnings(as.integer(str_trim(parts[2]))))
    })
  } else {
    tibble(info = character(), value = integer())
  }
  
  
  # --- Early return when full = FALSE ---
  if (!full) {
    return(list(
      model_name    = model_name,
      criteria      = criteria,
      ll_method     = ll_method,
      covariates    = covariates,
      dataset_info  = dataset_info,
      pop_params    = NULL,
      indiv_params  = NULL,
      fisher_matrix = NULL,
      eigen_values  = NULL
    ))
  }
  
  
  # --- Section boundary indices ---
  pop_start    <- find_line(lines, "ESTIMATION OF THE POPULATION PARAMETERS")
  indiv_start  <- find_line(lines, "ESTIMATION OF THE INDIVIDUAL PARAMETERS")
  fisher_start <- find_line(lines, "ESTIMATION OF THE FISHER INFORMATION MATRIX")
  ll_start     <- find_line(lines, "ESTIMATION OF THE LOG-LIKELIHOOD")
  
  
  # --- Population parameters ---
  section_types <- c(
    "Fixed Effects",
    "Standard Deviation of the Random Effects",
    "Error Model Parameters"
  )
  
  pop_end <- coalesce(indiv_start, fisher_start, ll_start)
  
  type_idx <- which(
    str_detect(lines, paste(section_types, collapse = "|")) &
      seq_along(lines) > pop_start &
      seq_along(lines) < pop_end
  )
  type_names <- str_extract(lines[type_idx],
                            paste(section_types, collapse = "|"))
  
  elapsed_idx <- which(
    str_detect(lines, "^Elapsed time") &
      seq_along(lines) > pop_start &
      seq_along(lines) < pop_end
  )[1]
  
  type_idx <- c(type_idx, elapsed_idx)
  
  pop_list <- list()
  
  for (i in seq_along(type_names)) {
    blk <- lines[(type_idx[i] + 1):(type_idx[i + 1] - 1)]
    blk <- clean_block(blk)
    blk <- str_trim(blk)
    blk <- blk[!str_detect(blk, "Elapsed|Exploratory|Smoothing")]
    data_lines <- blk[str_detect(blk, ":") & str_detect(blk, "\\d")]
    if (length(data_lines) == 0) next
    
    rows <- map_dfr(data_lines, function(ln) {
      param <- str_trim(str_extract(ln, "^[^:]+"))
      nums  <- nums_after_name(str_remove(ln, "^[^:]+:\\s*"))
      tibble(
        parameter = param,
        estimate  = if (length(nums) >= 1) nums[1] else NA_real_,
        se_sa     = if (length(nums) >= 2) nums[2] else NA_real_,
        rse_pct   = if (length(nums) >= 3) nums[3] else NA_real_,
        ci_low    = if (length(nums) >= 4) nums[4] else NA_real_,
        ci_high   = if (length(nums) >= 5) nums[5] else NA_real_
      )
    })
    
    rows$type <- type_names[i]
    pop_list[[type_names[i]]] <- rows
  }
  
  pop_params <- bind_rows(pop_list) |> 
    select(where(function(x) any(!is.na(x)))) |> # Remove absent (empty) columns
    relocate(type)
  
  
  # --- Individual parameters ---
  indiv_params <- NULL
  
  if (!is.na(indiv_start)) {
    
    cd_start <- find_line(lines, "Conditional Distribution")
    cm_start <- find_line(lines, "Conditional Mode")
    indiv_end <- coalesce(fisher_start, ll_start)
    
    parse_indiv_block <- function(block_start, block_end) {
      if (is.na(block_start)) return(NULL)
      blk <- lines[(block_start + 1):(block_end - 1)]
      blk <- clean_block(blk)
      blk <- str_trim(blk)
      data_lines <- blk[
        str_detect(blk, ":") &
          str_detect(blk, "\\d") &
          !str_detect(blk, "Elapsed|Iterations|min|Q1|shrinkage")
      ]
      if (length(data_lines) == 0) return(NULL)
      
      map_dfr(data_lines, function(ln) {
        param <- str_trim(str_extract(ln, "^[^:]+"))
        nums  <- nums_after_name(str_remove(ln, "^[^:]+:\\s*"))
        tibble(
          parameter = param,
          indiv_min = if (length(nums) >= 1) nums[1] else NA_real_,
          q1        = if (length(nums) >= 2) nums[2] else NA_real_,
          median    = if (length(nums) >= 3) nums[3] else NA_real_,
          q3        = if (length(nums) >= 4) nums[4] else NA_real_,
          indiv_max = if (length(nums) >= 5) nums[5] else NA_real_,
          shrinkage = if (length(nums) >= 6) nums[6] else NA_real_
        )
      })
    }
    
    # Conditional Distribution ends where Conditional Mode begins (or at
    # the section end if only one method was run)
    cd_end <- if (!is.na(cm_start) && !is.na(cd_start) && cm_start > cd_start) {
      cm_start
    } else {
      indiv_end
    }
    
    indiv_cd <- parse_indiv_block(cd_start, cd_end)
    indiv_cm <- parse_indiv_block(cm_start, indiv_end)
    
    if (!is.null(indiv_cd)) indiv_cd$method <- "Conditional Distribution"
    if (!is.null(indiv_cm)) indiv_cm$method <- "Conditional Mode"
    
    indiv_params <- bind_rows(indiv_cd, indiv_cm)
  }
  
  
  # --- Fisher information matrix ---
  fisher_matrix <- NULL
  eigen_values  <- NULL
  
  if (!is.na(fisher_start)) {
    
    corr_start  <- find_line(lines, "Correlation Matrix")
    eigen_idx   <- find_line(lines, "Eigen values")
    
    if (!is.na(corr_start) && !is.na(eigen_idx)) {
      
      # Eigen values ---
      eigen_nums <- nums_after_name(lines[eigen_idx])
      eigen_values <- tibble(
        eigen_min = eigen_nums[1],
        eigen_max = eigen_nums[2],
        ratio     = eigen_nums[3]
      )
      
      # Correlation matrix ---
      # Each row: "   param_name   val1  val2  ...  1"
      # Numbers must be extracted from the portion after the param name to
      # avoid picking up digits embedded in names (e.g. "ka1", "beta_F3").
      corr_lines <- lines[(corr_start + 1):(eigen_idx - 1)]
      corr_lines <- clean_block(corr_lines)
      corr_lines <- corr_lines[str_detect(corr_lines, "\\d")]
      corr_lines <- str_trim(corr_lines)
      
      if (length(corr_lines) > 0) {
        
        param_names <- str_extract(corr_lines, "^[A-Za-z_][A-Za-z0-9_]*")
        n <- length(param_names)
        
        mat <- matrix(NA_real_, nrow = n, ncol = n,
                      dimnames = list(param_names, param_names))
        
        for (i in seq_len(n)) {
          # Parse only what comes after the parameter name
          vals <- nums_after_name(corr_lines[i])
          # Row i of a lower-triangular matrix has exactly i values (cols 1..i)
          # Guard against malformed lines
          k <- min(length(vals), i)
          for (j in seq_len(k)) {
            mat[i, j] <- vals[j]
            mat[j, i] <- vals[j]
          }
        }
        
        fisher_matrix <- as.data.frame(mat)
      }
    }
  }
  
  
  # --- Return ---
  list(
    model_name    = model_name,
    criteria      = criteria,
    ll_method     = ll_method,
    covariates    = covariates,
    dataset_info  = dataset_info,
    pop_params    = pop_params,
    indiv_params  = indiv_params,
    fisher_matrix = fisher_matrix,
    eigen_values  = eigen_values
  )
}
