  ## MODULE 3: build_x_data()       ============================================
    ### ---- HELPERS: critical_ofv() / .lrt_p() ---
#' p-value of an OFV improvement `delta` (> 0 = significant direction) on `df`
#'
#' degrees of freedom; a negative delta gives p = 1
#'
#' @noRd
.lrt_p <- function(delta, df) pchisq(pmax(delta, 0), df = df, lower.tail = FALSE)


    ### ---- MODULE 3a: build_scm_data() ---
#' From the master models table, builds nodes and edges for the forward SCM
#'
#' @param models_df tibble from read_models()
#' @param alpha LRT significance level, REQUIRED. Must equal the level used by the run being drawn (e.g. set_SCM$threshold$lrt[1]). The critical OFV drop is critical_ofv(alpha, df), df = number of beta coefficients added vs the true parent (3.84 for df = 1)
#'
#' @return list: $nodes          : tibble with added_cov, delta_OFV, df, p_value, significant, is_best, iter $edges          : tibble with from / to model_id pairs $last_best_name : model_name of the best model in the final iteration
#'
#' @export
build_scm_data <- function(models_df, alpha) {
  if (missing(alpha)) stop("`alpha` is required: use the SAME LRT level as the run being drawn ",
                           "(e.g. set_SCM$threshold$lrt[1] for a Monolix forward step, or the alpha given to ",
                           "backward_elimination()).", call. = FALSE)
  
  iters    <- sort(unique(models_df$iter))
  iter_min <- min(iters)
  
  best_per_iter <- models_df |>
    group_by(iter) |>
    slice_min(OFV, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(iter, best_model_id = model_id, best_OFV = OFV)
  
  # Find minimum covariate-distance parent among all previous models
  .find_alt_parent <- function(mid, covs) {
    candidates <- models_df |> filter(model_id < mid)
    if (nrow(candidates) == 0) return(NA_integer_)
    distances <- map_int(candidates$covariates, function(cc) {
      length(setdiff(covs, cc)) + length(setdiff(cc, covs))
    })
    min_dist <- min(distances)
    best_idx <- which(distances == min_dist)
    # Among equidistant candidates, prefer the most recent (highest model_id)
    candidates$model_id[best_idx[which.max(candidates$model_id[best_idx])]]
  }
  
  # Format covariate change label relative to a parent covariate set
  .cov_label <- function(covs, parent_covs, current_iter) {
    if (current_iter == iter_min) return(paste(covs, collapse = ", "))
    added   <- setdiff(covs, parent_covs)
    removed <- setdiff(parent_covs, covs)
    if (length(added) == 1 && length(removed) == 0) return(paste(added))
    if (length(added) == 0 && length(removed) == 1) return(paste0("-", removed))
    parts <- c(
      if (length(added)   > 0) paste0("+", paste(added,   collapse = ", ")),
      if (length(removed) > 0) paste0("-", paste(removed, collapse = ", "))
    )
    paste(parts, collapse = " / ")
  }
  
  # Compute true parent for each model
  # Standard (distance==1): best of previous iter
  # Alternative (distance>1): closest prior model by covariate overlap
  parent_lookup <- map_dfr(seq_len(nrow(models_df)), function(i) {
    mid          <- models_df$model_id[i]
    covs         <- models_df$covariates[[i]]
    current_iter <- models_df$iter[i]
    
    if (current_iter == iter_min)
      return(tibble(model_id = mid, true_parent_id = NA_integer_, is_alt = FALSE))
    
    prev_iter <- max(best_per_iter$iter[best_per_iter$iter < current_iter])
    std_pid   <- best_per_iter$best_model_id[best_per_iter$iter == prev_iter]
    std_covs  <- models_df$covariates[[which(models_df$model_id == std_pid)]]
    std_dist  <- length(setdiff(covs, std_covs)) + length(setdiff(std_covs, covs))
    
    if (std_dist == 1) {
      tibble(model_id = mid, true_parent_id = std_pid, is_alt = FALSE)
    } else {
      tibble(model_id = mid, true_parent_id = .find_alt_parent(mid, covs), is_alt = TRUE)
    }
  })
  
  nodes <- models_df |>
    left_join(parent_lookup, by = "model_id") |>
    mutate(
      parent_OFV_correct = map_dbl(true_parent_id, function(pid) {
        if (is.na(pid)) return(NA_real_)
        models_df$OFV[models_df$model_id == pid]
      }),
      
      delta_OFV = OFV - parent_OFV_correct,
      
      # df = number of beta coefficients added/removed vs the true parent
      parent_n_cov = map_int(true_parent_id, function(pid) {
        if (is.na(pid)) return(NA_integer_)
        models_df$n_cov[models_df$model_id == pid]
      }),
      df      = pmax(abs(n_cov - parent_n_cov), 1L),
      p_value = .lrt_p(-delta_OFV, df),
      
      added_cov = map2_chr(model_id, covariates, function(mid, covs) {
        current_iter <- models_df$iter[models_df$model_id == mid]
        if (length(covs) == 0) return("Base model")
        if (current_iter == iter_min) return(paste(covs, collapse = ", "))
        pid <- parent_lookup$true_parent_id[parent_lookup$model_id == mid]
        if (is.na(pid)) return(paste(covs, collapse = ", "))
        parent_covs <- models_df$covariates[[which(models_df$model_id == pid)]]
        .cov_label(covs, parent_covs, current_iter)
      }),
      
      significant       = !is.na(delta_OFV) & (delta_OFV < -critical_ofv(alpha, df)),
      is_best           = model_id %in% best_per_iter$best_model_id & significant,
      parent_id_correct = true_parent_id
    ) |>
    select(model_id, model_name, iter, OFV, AIC, BICc, BIC,
           covariates, n_cov, added_cov, delta_OFV,
           significant, is_best, parent_id_correct, df, p_value)
  
  # Standard edges: best(iter k) -> standard children(iter k+1)
  std_edges <- map_dfr(iters[-length(iters)], function(w) {
    from_id <- best_per_iter$best_model_id[best_per_iter$iter == w]
    to_ids  <- parent_lookup |>
      filter(model_id %in% models_df$model_id[models_df$iter == w + 1],
             !is_alt) |>
      pull(model_id)
    if (length(to_ids) == 0) return(tibble(from = integer(), to = integer()))
    tibble(from = from_id, to = to_ids)
  })
  
  # Alternative edges: true parent -> alternative branch models
  alt_edges <- parent_lookup |>
    filter(is_alt, !is.na(true_parent_id)) |>
    transmute(from = true_parent_id, to = model_id)
  
  edges <- bind_rows(std_edges, alt_edges)
  
  last_best <- nodes |> filter(is_best) |>
    slice_max(iter, n = 1, with_ties = FALSE)
  last_best_name <- if (nrow(last_best) > 0) last_best$model_name else NA_character_
  
  list(nodes = nodes, edges = edges,
       alpha = alpha, type = "forward",
       last_best_name = last_best_name)
}



    ### ---- MODULE 3e: read_model_building_report() / get_monolix_best_model() ---
#' Read the model building report written by Monolix
#'
#' Monolix writes a plain-text report of a model building run (modelBuilding.txt): for every model, its -2LL, its BICc and the parameter x covariate matrix of the relations it contains (an X marks an active relation), then the line "Best Model: n". This function parses it, so the model retained by Monolix can be known without loading Monolix or the results object.
#'
#' @param path The modelBuilding.txt file, or the folder that contains it.
#'
#' @return list with `models` (tibble: model_id, ofv, bicc, covariates, n_rel; covariates are "<parameter>_<covariate>" labels, like in the summaries) and `best_model_id` (the number after "Best Model:", NA if absent).
#'
#' @export
read_model_building_report <- function(path) {
  
  if (dir.exists(path)) {
    found <- list.files(path, pattern = "^modelBuilding.*\\.txt$", ignore.case = TRUE, full.names = TRUE)
    if (length(found) == 0) stop("No modelBuilding.txt in ", path, call. = FALSE)
    path <- found[1]
  }
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  
  lines  <- str_remove(readLines(path, warn = FALSE), "\r$")
  starts <- which(str_detect(lines, "^Model\\s+\\d+\\s*$"))
  if (length(starts) == 0) stop("No 'Model k' block found in ", path, call. = FALSE)
  ends   <- c(starts[-1] - 1L, length(lines))
  
  number_after <- function(block, label) {
    hit <- str_subset(block, paste0("^\\s*", label, "\\s*:"))
    if (length(hit) == 0) return(NA_real_)
    as.numeric(str_extract(hit[1], "-?[0-9.]+([eE][-+]?[0-9]+)?\\s*$"))
  }
  
  models <- map2(starts, ends, function(a, b) {
    block <- lines[a:b]
    
    # Matrix: header line lists the covariates, each row "| param | cell | cell |"
    cells <- function(x) str_trim(str_split(x, fixed("|"))[[1]])
    header_i <- which(str_detect(block, "^\\s+\\|"))[1]
    covs     <- cells(block[header_i]); covs <- covs[nzchar(covs)]
    rows     <- block[str_detect(block, "^\\|\\s*\\S+\\s*\\|")]
    active <- map(rows, function(r) {
      parts <- cells(r)[-1]                                    # first piece is empty
      param <- parts[1]
      marks <- c(parts[-1], rep("", length(covs)))[seq_along(covs)]
      if (any(nzchar(marks))) paste(param, covs[nzchar(marks)], sep = "_") else character(0)
    }) |> unlist()
    
    tibble(model_id   = as.integer(str_extract(lines[a], "\\d+")),
           ofv        = number_after(block, "-2 \\* LL"),
           bicc       = number_after(block, "BICc"),
           covariates = list(as.character(active)),
           n_rel      = length(active))
  }) |> list_rbind()
  
  best <- str_subset(lines, "^\\s*Best Model\\s*:")
  list(models        = models,
       best_model_id = if (length(best)) as.integer(str_extract(best[1], "\\d+")) else NA_integer_)
}
#' Model retained by Monolix at the end of a model building run
#'
#' Returns the model flagged as best in the report. If a table of models is given, the model name is looked up in it; when the best model has no folder of its own (it can happen), the run with the same set of relations is used instead.
#'
#' @param report Output of read_model_building_report().
#' @param models_df Optional table from read_models(), with covariate labels cleaned by clean_cov_labels(). Used to return the model name.
#'
#' @return list with `model_id`, `covariates`, `model_name` and `matched_by`.
#'
#' @export
get_monolix_best_model <- function(report, models_df = NULL) {
  
  best <- report$models |> filter(model_id == report$best_model_id)
  if (nrow(best) != 1) stop("Best model not found in the report", call. = FALSE)
  out <- list(model_id = best$model_id, covariates = best$covariates[[1]],
              model_name = NA_character_, matched_by = NA_character_)
  
  if (!is.null(models_df)) {
    by_id <- models_df |> filter(model_id == best$model_id)
    if (nrow(by_id) == 1) {
      out$model_name <- by_id$model_name; out$matched_by <- "model number"
    } else {
      key  <- function(x) paste(sort(as.character(x)), collapse = "|")
      same <- models_df |> filter(map_chr(covariates, key) == key(best$covariates[[1]])) |> arrange(model_id)
      if (nrow(same) == 0) stop("No run with the covariate set of the best model (",
                                paste(best$covariates[[1]], collapse = ", "), ")", call. = FALSE)
      out$model_id <- same$model_id[1]; out$model_name <- same$model_name[1]
      out$matched_by <- "covariate set (the best model has no folder of its own)"
    }
  }
  out
}



    ### ---- MODULE 3b: build_cossac_data() ---
#' Flowchart data for a sequential covariate search (chain version)
#'
#' Reads sequential model comparisons where the number of covariates is not strictly monotone. The parent of model N is model N-1 (testing order); added or removed covariates are detected by comparing with the parent covariate set.
#'
#' @param models_df Table from read_models(), covariate labels cleaned by clean_cov_labels().
#' @param alpha LRT significance level, REQUIRED. Must equal the level used by the run being drawn.
#'
#' @return list with `nodes`, `edges`, `alpha`, `type` and `last_best_name`.
#'
#' @export
build_cossac_data_linear <- function(models_df, alpha) {
  if (missing(alpha)) stop("`alpha` is required: use the SAME LRT level as the run being drawn ",
                           "(e.g. set_SCM$threshold$lrt[1] for a Monolix forward step, or the alpha given to ",
                           "backward_elimination()).", call. = FALSE)
  
  models_df <- models_df |> arrange(model_id)
  
  nodes <- models_df |>
    mutate(
      parent_id = lag(model_id),  # parent = previous model in sequence
      
      parent_OFV_correct = map_dbl(parent_id, function(pid) {
        if (is.na(pid)) return(NA_real_)
        models_df$OFV[models_df$model_id == pid]
      }),
      
      parent_id_correct = parent_id,
      
      delta_OFV = OFV - parent_OFV_correct,
      df        = pmax(abs(n_cov - lag(n_cov)), 1L),
      p_value   = .lrt_p(-delta_OFV, df),
      
      # Detect what changed vs parent: addition or removal
      added_cov = map2_chr(model_id, covariates, function(mid, covs) {
        pid <- models_df$model_id[models_df$model_id == mid - 1]
        if (length(pid) == 0) return(paste(covs, collapse = ", "))
        parent_covs <- models_df$covariates[models_df$model_id == pid][[1]]
        added   <- setdiff(covs, parent_covs)
        removed <- setdiff(parent_covs, covs)
        label <- c(
          if (length(added)   > 0) paste0("+", paste(added,   collapse = ", ")),
          if (length(removed) > 0) paste0("-", paste(removed, collapse = ", "))
        )
        if (length(label) == 0) paste(covs, collapse = ", ")
        else paste(label, collapse = " / ")
      }),
      
      significant = !is.na(delta_OFV) & (delta_OFV < -critical_ofv(alpha, df)),
      
      # iter = position in sequence (for layout purposes only)
      iter    = row_number() - 1L,
      is_best = FALSE  # no best-of-wave concept in COSSAC
    ) |>
    select(model_id, model_name, iter, OFV, AIC, BICc, BIC,
           covariates, n_cov, added_cov, delta_OFV,
           significant, is_best, parent_id_correct, df, p_value)
  
  # Edges: each model -> next model (linear chain)
  edges <- tibble(
    from = models_df$model_id[-nrow(models_df)],
    to   = models_df$model_id[-1]
  )
  
  list(nodes = nodes, edges = edges,
       alpha = alpha, type = "forward",
       last_best_name = nodes |> filter(significant) |>
         slice_max(iter, n = 1, with_ties = FALSE) |> pull(model_name))
}
#' Flowchart data for a COSSAC covariate search
#'
#' Sequential covariate model building: the parent of each model is the last accepted (significant) model before it, and non-significant models branch off their true parent instead of chaining. `is_best` is TRUE for the accepted models that form the backbone of the final model. The degrees of freedom of each test are the number of coefficients added or removed.
#'
#' @param models_df Table from read_models(), covariate labels cleaned by clean_cov_labels().
#' @param alpha LRT significance level, REQUIRED. Must equal the level of the Monolix run (e.g. the first value of set_COSSAC$threshold$lrt): the builder replays the decisions, it does not take them.
#'
#' @return list with `nodes`, `edges`, `alpha`, `type` and `last_best_name`.
#'
#' @export
build_cossac_data <- function(models_df, alpha) {
  if (missing(alpha)) stop("`alpha` is required: use the SAME LRT level as the run being drawn ",
                           "(e.g. set_SCM$threshold$lrt[1] for a Monolix forward step, or the alpha given to ",
                           "backward_elimination()).", call. = FALSE)
  
  models_df <- models_df |> arrange(model_id)
  n         <- nrow(models_df)
  
  # Sequential scan: track the last accepted model at each step
  true_parent_idx <- rep(NA_integer_, n)
  delta_ofv       <- rep(NA_real_,    n)
  significant     <- rep(FALSE,       n)
  df_vec          <- rep(NA_integer_, n)
  p_vec           <- rep(NA_real_,    n)
  last_acc_idx    <- 1L  # starting model is always the first
  
  for (i in 2:n) {
    d  <- models_df$OFV[i] - models_df$OFV[last_acc_idx]
    df <- max(abs(models_df$n_cov[i] - models_df$n_cov[last_acc_idx]), 1L)
    true_parent_idx[i] <- last_acc_idx
    delta_ofv[i]       <- d
    df_vec[i]          <- df
    p_vec[i]           <- .lrt_p(-d, df)
    if (d < -critical_ofv(alpha, df)) {
      significant[i] <- TRUE
      last_acc_idx   <- i
    }
  }
  
  true_parent_id <- models_df$model_id[true_parent_idx]
  
  # Assign iter:
  # - Starting model        : iter 0
  # - Accepted models       : sequential position in acceptance chain (1, 2, 3...)
  # - Non-accepted models   : same iter as their true parent (branch off it)
  acc_counter <- 0L
  iter_vec    <- integer(n)
  
  for (i in seq_len(n)) {
    if (i == 1) {
      iter_vec[i] <- 0L
    } else if (significant[i]) {
      acc_counter  <- acc_counter + 1L
      iter_vec[i]  <- acc_counter
    } else {
      # Same iter as true parent
      iter_vec[i] <- iter_vec[true_parent_idx[i]]
    }
  }
  
  # Covariate label: what changed vs true parent
  added_cov <- map2_chr(seq_len(n), models_df$covariates, function(i, covs) {
    pidx <- true_parent_idx[i]
    if (is.na(pidx)) return(paste(covs, collapse = ", "))
    parent_covs <- models_df$covariates[[pidx]]
    added   <- setdiff(covs, parent_covs)
    removed <- setdiff(parent_covs, covs)
    parts <- c(
      if (length(added)   > 0) paste0("+", paste(added,   collapse = ", ")),
      if (length(removed) > 0) paste0("-", paste(removed, collapse = ", "))
    )
    if (length(parts) == 0) paste(covs, collapse = ", ") else paste(parts, collapse = " / ")
  })
  
  nodes <- models_df |>
    mutate(
      iter               = iter_vec,
      parent_id_correct  = true_parent_id,
      delta_OFV          = delta_ofv,
      added_cov          = added_cov,
      significant        = significant,
      df                 = df_vec,
      p_value            = p_vec,
      # Accepted models form the backbone — mark them as is_best for green display
      is_best            = significant
    ) |>
    select(model_id, model_name, iter, OFV, AIC, BICc, BIC,
           covariates, n_cov, added_cov, delta_OFV,
           significant, is_best, parent_id_correct, df, p_value)
  
  # Edges: true parent -> each model (skip starting model)
  edges <- nodes |>
    filter(!is.na(parent_id_correct)) |>
    transmute(from = parent_id_correct, to = model_id)
  
  last_best_name <- nodes |>
    filter(is_best) |>
    slice_max(iter, n = 1, with_ties = FALSE) |>
    pull(model_name)
  
  if (length(last_best_name) == 0) last_best_name <- nodes$model_name[1]
  
  list(nodes = nodes, edges = edges,
       alpha = alpha, type = "forward",
       last_best_name = last_best_name)
}




    ### ---- MODULE 3c: build_backward_data() ---
#' Builds nodes and edges for the backward elimination flowchart
#'
#' Logic (mirror of backward_elimination()):
#'   - Iteration 0  : full model
#'   - Iteration k  : candidates with one fewer relation than the model ADOPTED
#'                    at iteration k-1
#'   - delta_OFV    : positive = removal worsened fit (relation is useful)
#'   - df           : number of beta coefficients removed vs the parent, so a
#'                    multi-level categorical covariate is tested on its own df
#'   - p_value      : LRT p-value of the removal, 1 - pchisq(delta_OFV, df)
#'   - significant  : p_value < alpha -> relation should be KEPT
#'   - is_best      : adopted candidate = highest p-value of the iteration (ties:
#'                    lowest OFV), only flagged if NOT significant (i.e. the
#'                    relation can be dropped). This is the very rule applied by
#'                    backward_elimination(), so the chart shows the path the
#'                    engine actually followed.
#'
#' @param models_df tibble from read_models()
#' @param alpha LRT significance level to KEEP a relation, REQUIRED (the alpha given to backward_elimination(), or the threshold you used by hand)
#'
#' @return list: $nodes, $edges, $alpha, $type, $last_best_name
#'
#' @export
build_backward_data <- function(models_df, alpha) {
  if (missing(alpha)) stop("`alpha` is required: use the SAME LRT level as the run being drawn ",
                           "(e.g. set_SCM$threshold$lrt[1] for a Monolix forward step, or the alpha given to ",
                           "backward_elimination()).", call. = FALSE)
  
  models_df <- models_df |>
    mutate(n_rel = map_int(covariates, ~ length(unique(clean_cov_labels(.x)))),
           iter  = max(n_rel) - n_rel)
  
  iters  <- sort(unique(models_df$iter))
  parent <- NULL
  nodes  <- vector("list", length(iters))
  
  # Sequential scan: the parent of iteration k is the model adopted at k-1
  for (k in seq_along(iters)) {
    it <- models_df |> filter(iter == iters[k])
    
    if (is.null(parent)) {
      # iteration 0 = full model (lowest OFV if it was run more than once)
      it <- it |>
        mutate(parent_id_correct = NA_integer_, delta_OFV = NA_real_,
               df = NA_integer_, p_value = NA_real_,
               significant = FALSE, is_best = FALSE)
      parent <- it |> slice_min(OFV, n = 1, with_ties = FALSE)
    } else {
      it <- it |>
        mutate(parent_id_correct = parent$model_id,
               delta_OFV   = OFV - parent$OFV,
               df          = pmax(parent$n_cov - n_cov, 1L),
               p_value     = .lrt_p(delta_OFV, df),
               significant = p_value < alpha)
      best_id <- it |> arrange(desc(p_value), OFV) |> slice(1) |> pull(model_id)
      it      <- it |> mutate(is_best = model_id == best_id & !significant)
      parent  <- it |> filter(model_id == best_id)
    }
    nodes[[k]] <- it
  }
  
  nodes <- bind_rows(nodes) |>
    mutate(added_cov = map2_chr(covariates, parent_id_correct, function(covs, pid) {
      if (is.na(pid)) return("Full model")
      parent_covs <- models_df$covariates[[which(models_df$model_id == pid)]]
      dropped <- unique(setdiff(parent_covs, covs))
      if (length(dropped) == 0) paste(covs, collapse = ", ") else paste(dropped, collapse = ", ")
    })) |>
    select(model_id, model_name, iter, OFV, AIC, BICc, BIC,
           covariates, n_cov, added_cov, delta_OFV,
           significant, is_best, parent_id_correct, df, p_value)
  
  edges <- nodes |>
    filter(!is.na(parent_id_correct)) |>
    transmute(from = parent_id_correct, to = model_id)
  
  # If nothing was dropped, the full model (iter 0) is the final model
  last_best_name <- if (any(nodes$is_best)) {
    nodes |> filter(is_best) |> slice_max(iter, n = 1, with_ties = FALSE) |> pull(model_name)
  } else {
    nodes |> filter(iter == 0) |> slice_min(OFV, n = 1, with_ties = FALSE) |> pull(model_name)
  }
  
  list(nodes = nodes, edges = edges,
       alpha = alpha, type = "backward",
       last_best_name = last_best_name)
}



    ### ---- MODULE 3d: merge_scm_data() ---
#' Merges forward and backward scm_data lists into a single structure for the
#'
#' combined flowchart (option 4).
#'
#' Layout logic:
#'   - Forward nodes:  x = iter (0 .. n)
#'   - Visual gap:     x = n+1 .. n+gap-1  (empty columns, creates white space)
#'   - Backward nodes: x = n+gap+iter  (iter 0 = Full model, iter 1.. = removals)
#'   - NO arrow between last forward best and backward iter=0 node — the gap
#'     is the visual separator, no connecting edge is added.
#'
#' Unique vertex ids (uid): "F_<model_id>" / "B_<model_id>"
#' This avoids collisions when forward and backward share model_id numbering.
#'
#' @param forward_data list from build_scm_data()
#' @param backward_data list from build_backward_data()
#' @param gap white-space columns between phases (default 2)
#'
#' @return list: $nodes                  : combined tibble (uid, phase, x_pos added) $edges                  : combined edges using uid $forward_data           : original forward list (for individual plot) $backward_data          : original backward list (for individual plot) $alpha_forward          : forward LRT level $alpha_backward         : backward LRT level $last_best_name         : final selected model (from backward)
#'
#' @export
merge_scm_data <- function(forward_data, backward_data, gap = 2) {
  
  max_fwd_x <- max(forward_data$nodes$iter)
  
  # Add phase, uid, x_pos to each phase
  fwd_nodes <- forward_data$nodes |>
    mutate(
      phase = "forward",
      uid   = paste0("F_", model_id),
      x_pos = iter
    )
  
  # Backward iter=0 (Full model) starts at max_fwd_x + gap
  bwd_nodes <- backward_data$nodes |>
    mutate(
      phase = "backward",
      uid   = paste0("B_", model_id),
      x_pos = max_fwd_x + gap + iter
    )
  
  all_nodes <- bind_rows(fwd_nodes, bwd_nodes)
  
  # Prefix edges with phase to match uid references
  fwd_edges <- forward_data$edges |>
    mutate(from = paste0("F_", from), to = paste0("F_", to))
  
  bwd_edges <- backward_data$edges |>
    mutate(from = paste0("B_", from), to = paste0("B_", to))
  
  # No edge between forward last best and backward iter=0 — gap is intentional
  all_edges <- bind_rows(fwd_edges, bwd_edges)
  
  list(
    nodes                  = all_nodes,
    edges                  = all_edges,
    forward_data           = forward_data,
    backward_data          = backward_data,
    alpha_forward          = forward_data$alpha,
    alpha_backward         = backward_data$alpha,
    last_best_name         = backward_data$last_best_name
  )
}
