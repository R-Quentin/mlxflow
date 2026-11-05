  ## MODULE 0: renumber_model_folder()    ============================================
#' Renames all folders and files in a ModelBuilding_B directory so their
#'
#' numbering continues from ModelBuilding_A.
#'
#' Renaming order (bottom-up to avoid path conflicts):
#'   1. Update exportpath in .mlxtran file content
#'   2. Rename .mlxtran file
#'   3. Rename inner homonymous subfolder
#'   4. Rename outer run folder
#'
#' @param b_path path to ModelBuilding_B folder
#' @param offset integer added to each B run number
#' @param dry_run if TRUE (default), only print planned operations
#'
#' @export
renumber_model_folder <- function(b_path, offset, dry_run = TRUE) {
  
  if (dry_run) message("--- DRY RUN: no files will be modified ---\n")
  
  pad3 <- function(n) sprintf("%03d", n)
  
  b_dirs  <- dir_ls(b_path, type = "directory") |>
    keep(~ str_detect(basename(.x), "^[A-Za-z]+_\\d+$"))
  
  old_ids <- as.integer(str_extract(basename(b_dirs), "\\d+"))
  new_ids <- old_ids + offset
  
  name_map <- tibble(
    old_num  = old_ids,
    new_num  = new_ids,
    old_name = paste0(str_extract(basename(b_dirs), "^[A-Za-z]+_"), pad3(old_ids)),
    new_name = paste0(str_extract(basename(b_dirs), "^[A-Za-z]+_"), pad3(new_ids))
  )
  
  message(glue("  Offset: +{offset}"))
  message(glue("  B runs: {min(old_ids)} to {max(old_ids)}"))
  message(glue("  Will become: {min(new_ids)} to {max(new_ids)}\n"))
  
  log <- list()
  
  for (i in rev(seq_len(nrow(name_map)))) {
    
    old_nm  <- name_map$old_name[i]
    new_nm  <- name_map$new_name[i]
    run_dir <- file.path(b_path, old_nm)
    
    # Step 1: update .mlxtran content
    mlxtran_old <- file.path(run_dir, paste0(old_nm, ".mlxtran"))
    mlxtran_new <- file.path(run_dir, paste0(new_nm, ".mlxtran"))
    
    if (file.exists(mlxtran_old)) {
      content  <- readLines(mlxtran_old, warn = FALSE)
      updated  <- str_replace_all(content, stringr::fixed(old_nm), new_nm)
      n_rep    <- sum(content != updated)
      
      if (!dry_run) writeLines(updated, mlxtran_old)
      
      log[[length(log) + 1]] <- tibble(
        step = "mlxtran_content", old = mlxtran_old,
        new = mlxtran_new, details = glue("{n_rep} line(s) updated")
      )
      message(glue("  [content] {old_nm}.mlxtran -> {n_rep} line(s) updated"))
    }
    
    # Step 2: rename .mlxtran file
    if (file.exists(mlxtran_old)) {
      if (!dry_run) file.rename(mlxtran_old, mlxtran_new)
      log[[length(log) + 1]] <- tibble(
        step = "rename_mlxtran", old = mlxtran_old, new = mlxtran_new, details = ""
      )
      message(glue("  [rename]  {old_nm}.mlxtran  ->  {new_nm}.mlxtran"))
    }
    
    # Step 3: rename inner subfolder
    inner_old <- file.path(run_dir, old_nm)
    inner_new <- file.path(run_dir, new_nm)
    
    if (dir.exists(inner_old)) {
      if (!dry_run) file.rename(inner_old, inner_new)
      log[[length(log) + 1]] <- tibble(
        step = "rename_inner_dir", old = inner_old, new = inner_new, details = ""
      )
      message(glue("  [rename]  {old_nm}/{old_nm}/  ->  {old_nm}/{new_nm}/"))
    }
    
    # Step 4: rename outer run folder
    outer_new <- file.path(b_path, new_nm)
    if (!dry_run) file.rename(run_dir, outer_new)
    log[[length(log) + 1]] <- tibble(
      step = "rename_outer_dir", old = run_dir, new = outer_new, details = ""
    )
    message(glue("  [rename]  {old_nm}/  ->  {new_nm}/"))
  }
  
  log_df <- bind_rows(log)
  
  if (dry_run) {
    message("\n--- DRY RUN complete: re-run with dry_run = FALSE to apply ---")
  } else {
    message(glue("\n  Done. {nrow(log_df)} operations completed."))
  }
  
  invisible(log_df)
}
