  ## MODULE 4: plot_flowchart() ============================================

    ### ---- MODULE 4a: plot_flowchart() ---
#' Renders forward SCM or backward elimination flowchart using ggraph
#'
#' @param scm_data list from build_scm_data() or build_backward_data()
#' @param title plot title (the LRT level shown in the legend is taken from scm_data$alpha, set in build_*_data())
#' @param iter_label header prefix, e.g. "Iteration"
#' @param best_label legend label for best-of-iteration nodes
#' @param sig_label legend label for significant nodes; "Not <sig_label>" auto-derived
#' @param base_label label displayed inside the base/full model node
#' @param iter0_label header of the first column (default: "Full model" for a backward chart, "Base" otherwise)
#' @param node_size size passed to geom_node_point (default 26)
#' @param node_spacing vertical spacing multiplier between nodes (default 1.2)
#' @param colors named list with one sub-list per category (base, best, sig, not_sig) and an arrow sub-list. Each category accepts: fill, border, font. Arrow accepts: best, sig, not_sig.
#'
#' @export
plot_flowchart <- function(
    scm_data,
    title         = "SCM Flowchart",
    iter_label    = "Iteration",
    best_label    = "Best of iteration",
    sig_label     = "Significant",
    base_label    = "Base model\n(no covariates)",
    iter0_label   = NULL,
    node_size     = 25,
    node_spacing  = 1.2,
    colors        = list(
      base    = list(fill = "lightblue", border = "steelblue4", font = "black"),
      sig     = list(fill = "white",     border = "gray30",     font = "black",  arrow = "black"),
      not_sig = list(fill = "gray90",    border = "gray30",     font = "grey30", arrow = "gray80"),
      best    = list(fill = "green4",    border = "darkgreen",  font = "white",  arrow = "darkgreen")
    )) {
  
  nodes <- scm_data$nodes
  edges <- scm_data$edges
  alpha <- scm_data$alpha
  if (is.null(iter0_label))
    iter0_label <- if (identical(scm_data$type, "backward")) "Full model" else "Base"

  g <- igraph::graph_from_data_frame(
    d        = edges,
    vertices = nodes |> select(model_id, everything()),
    directed = TRUE
  )
  
  # Manual layout: x = iter, y = vertically centered within iter
  layout_df <- nodes |>
    group_by(iter) |>
    mutate(
      n_in_iter    = n(),
      rank_in_iter = row_number(),
      x = iter,
      y = -(rank_in_iter - (n_in_iter + 1) / 2) * node_spacing
    ) |>
    ungroup() |>
    arrange(model_id)
  
  layout_mat <- as.matrix(layout_df[, c("x", "y")])
  rownames(layout_mat) <- as.character(layout_df$model_id)
  
  vertex_order <- as.integer(igraph::V(g)$name)
  layout_mat   <- layout_mat[as.character(vertex_order), , drop = FALSE]
  
  idx     <- match(vertex_order, nodes$model_id)
  # is_base <- nodes$iter[idx] == 0
  is_base <- nodes$iter[idx] == min(nodes$iter)
  is_best <- nodes$is_best[idx]
  is_sig  <- nodes$significant[idx]
  
  # Node aesthetics — base takes priority over best to avoid iter-0 conflicts
  node_fill <- case_when(
    is_base ~ colors$base$fill,
    is_best ~ colors$best$fill,
    is_sig  ~ colors$sig$fill,
    TRUE    ~ colors$not_sig$fill
  )
  
  node_border <- case_when(
    is_base ~ colors$base$border,
    is_best ~ colors$best$border,
    is_sig  ~ colors$sig$border,
    TRUE    ~ colors$not_sig$border
  )
  
  font_color <- case_when(
    is_base ~ colors$base$font,
    is_best ~ colors$best$font,
    is_sig  ~ colors$sig$font,
    TRUE    ~ colors$not_sig$font
  )
  
  # Arrow color based on destination node status
  edge_colors <- map_chr(seq_len(igraph::ecount(g)), function(i) {
    to_id  <- as.integer(ends(g, i)[2])
    to_idx <- match(to_id, nodes$model_id)
    if      (nodes$is_best[to_idx])      colors$best$arrow
    else if (nodes$significant[to_idx])  colors$sig$arrow
    else                                 colors$not_sig$arrow
  })
  
  # Node labels
  node_labels  <- nodes$added_cov[idx]
  delta_labels <- ifelse(
    is.na(nodes$delta_OFV[idx]), "",
    sprintf("\u0394OFV = %.2f", nodes$delta_OFV[idx])
  )
  full_labels          <- paste0(node_labels, "\n", delta_labels)
  full_labels[is_base] <- base_label
  
  iter_headers <- layout_df |>
    group_by(iter) |>
    summarise(x = first(x), .groups = "drop") |>
    mutate(label = ifelse(iter == 0, iter0_label, paste(iter_label, iter)))
  
  # Cap radius: in ggplot2 size is in mm (diameter), so radius = node_size / 2
  cap_mm <- node_size / 2
  
  y_header    <- max(layout_df$y) + 1.2
  y_separator <- max(layout_df$y) + 0.85
  y_legend    <- min(layout_df$y) - 1.5
  x_legend    <- min(iter_headers$x)
  
  p <- ggraph(g, layout = layout_mat) +
    
    geom_edge_link(
      aes(color = edge_colors),
      arrow     = arrow(length = unit(2.5, "mm"), type = "closed"),
      end_cap   = circle(cap_mm, "mm"),
      start_cap = circle(cap_mm, "mm"),
      width     = 0.6,
      alpha     = 0.85
    ) +
    scale_edge_color_identity() +
    
    geom_node_point(
      shape  = 22,
      size   = node_size,
      fill   = node_fill,
      color  = node_border,
      stroke = 1.2
    ) +
    
    geom_node_text(
      aes(label = full_labels),
      size       = 2.5,
      lineheight = 1.15,
      color      = font_color
    ) +
    
    annotate("text",
             x = iter_headers$x, y = y_header,
             label = iter_headers$label,
             fontface = "bold", size = 3.5, color = "gray20") +
    
    annotate("segment",
             x    = iter_headers$x - 0.45, xend = iter_headers$x + 0.45,
             y    = y_separator,           yend = y_separator,
             color = "gray60", linewidth = 0.4) +
    
    annotate("point", shape = 22, size = 5,
             x = x_legend, y = y_legend,
             fill = colors$best$fill, color = colors$best$border) +
    annotate("text", hjust = 0, size = 2.8, color = "gray20",
             x = x_legend + 0.15, y = y_legend,
             label = best_label) +
    
    annotate("point", shape = 22, size = 5,
             x = x_legend, y = y_legend - 0.45,
             fill = colors$sig$fill, color = colors$sig$border) +
    annotate("text", hjust = 0, size = 2.8, color = "gray20",
             x = x_legend + 0.15, y = y_legend - 0.45,
             label = paste0(sig_label, " (LRT p < ", alpha, ")")) +
    
    annotate("point", shape = 22, size = 5,
             x = x_legend, y = y_legend - 0.9,
             fill = colors$not_sig$fill, color = colors$not_sig$border) +
    annotate("text", hjust = 0, size = 2.8, color = "gray20",
             x = x_legend + 0.15, y = y_legend - 0.9,
             label = paste0("Not ", tolower(sig_label))) +
    
    labs(title = title) +
    
    theme_void() +
    theme(
      plot.title      = element_text(size = 14, face = "bold",
                                     color = "gray10", margin = margin(b = 15)),
      plot.margin     = margin(25, 25, 25, 25),
      plot.background = element_rect(fill = "gray98", color = NA)
    )
  
  p
}



  ### ---- MODULE 4b: plot_flowchart_combined() ---
#' Renders forward SCM and backward elimination as a single plot
#'
#' Visual conventions:
#'   - Left section  : forward SCM, iterations labeled "Iteration 1..n"
#'   - Dashed separator with phase labels ("Forward SCM" / "Backward Elimination")
#'   - Right section : backward, iterations labeled independently "Iteration 1..m"
#'     The backward iter=0 (Full model) is labeled as such, not "Iteration N+1"
#'   - Each phase has its own color palette and legend
#'   - No arrow between the two sections
#'
#' @param merged_data list from merge_scm_data()
#' @param title overall plot title
#' @param forward_title Title of the forward section.
#' @param forward_iter_label header prefix for forward columns (default "Iteration")
#' @param forward_base_label text inside the forward base model node
#' @param forward_sig_label significance label for forward legend
#' @param forward_best_label best-of-iter label for forward legend
#' @param backward_title Title of the backward section.
#' @param backward_iter_label header prefix for backward columns (default "Iteration")
#' @param backward_base_label text inside the backward full model node
#' @param backward_sig_label significance label for backward legend
#' @param backward_best_label best-of-iter label for backward legend
#' @param node_size ggraph point size (default 25)
#' @param node_spacing vertical spacing multiplier (default 1.2)
#' @param colors_forward colors list (same structure as plot_flowchart)
#' @param colors_backward colors list for backward phase
#' @param separator_color dashed separator line color
#'
#' @export
plot_flowchart_combined <- function(
    merged_data,
    title                = "SCM Flowchart",
    forward_title        = "Forward inclusion (SCM)",
    forward_iter_label   = "Iteration",
    forward_base_label   = "Base model\n(no covariates)",
    forward_sig_label    = "Significant",
    forward_best_label   = "Best of iteration",
    backward_title       = "Backward deletion",
    backward_iter_label  = "Iteration",
    backward_base_label  = "Full model\n(all covariates)",
    backward_sig_label   = "Kept",
    backward_best_label  = "Dropped (worst)",
    node_size            = 25,
    node_spacing         = 1.2,
    colors_forward       = list(
      base    = list(fill = "lightblue",   border = "steelblue4", font = "black"),
      sig     = list(fill = "white",       border = "gray30",     font = "black",  arrow = "black"),
      not_sig = list(fill = "gray90",      border = "gray30",     font = "grey30", arrow = "gray80"),
      best    = list(fill = "green4",      border = "darkgreen",  font = "white",  arrow = "darkgreen")
    ),
    colors_backward      = list(
      base    = list(fill = "lightgreen",  border = "darkgreen",  font = "black"),
      sig     = list(fill = "lightyellow", border = "gray30",     font = "black",  arrow = "black"),
      not_sig = list(fill = "brown3",      border = "gray30",     font = "black",  arrow = "gray40"),
      best    = list(fill = "steelblue",   border = "steelblue4", font = "white",  arrow = "steelblue")
    ),
    separator_color      = "gray50") {
  
  nodes <- merged_data$nodes
  edges <- merged_data$edges
  alpha_fwd <- merged_data$alpha_forward
  alpha_bwd <- merged_data$alpha_backward
  
  # Build igraph using uid as vertex identifier
  g <- igraph::graph_from_data_frame(
    d        = edges,
    vertices = nodes |> select(uid, everything()),
    directed = TRUE
  )
  
  # Manual layout: x = x_pos, y = centered within x_pos column
  layout_df <- nodes |>
    group_by(x_pos) |>
    mutate(
      n_in_col    = n(),
      rank_in_col = row_number(),
      x = x_pos,
      y = -(rank_in_col - (n_in_col + 1) / 2) * node_spacing
    ) |>
    ungroup() |>
    arrange(uid)
  
  layout_mat <- as.matrix(layout_df[, c("x", "y")])
  rownames(layout_mat) <- layout_df$uid
  
  vertex_order <- igraph::V(g)$name
  layout_mat   <- layout_mat[vertex_order, , drop = FALSE]
  
  idx        <- match(vertex_order, nodes$uid)
  node_phase <- nodes$phase[idx]
  node_iter  <- nodes$iter[idx]
  is_best    <- nodes$is_best[idx]
  is_sig     <- nodes$significant[idx]
  
  min_fwd_iter <- min(nodes$iter[nodes$phase == "forward"])
  is_base_fwd  <- node_phase == "forward"  & node_iter == min_fwd_iter
  is_base_bwd  <- node_phase == "backward" & node_iter == 0
  
  # Node aesthetics: base takes priority over best to avoid iter-0 conflicts
  node_fill <- case_when(
    is_base_fwd                              ~ colors_forward$base$fill,
    is_base_bwd                              ~ colors_backward$base$fill,
    node_phase == "forward"  & is_best       ~ colors_forward$best$fill,
    node_phase == "backward" & is_best       ~ colors_backward$best$fill,
    node_phase == "forward"  & is_sig        ~ colors_forward$sig$fill,
    node_phase == "backward" & is_sig        ~ colors_backward$sig$fill,
    node_phase == "forward"                  ~ colors_forward$not_sig$fill,
    TRUE                                     ~ colors_backward$not_sig$fill
  )
  
  node_border <- case_when(
    is_base_fwd                              ~ colors_forward$base$border,
    is_base_bwd                              ~ colors_backward$base$border,
    node_phase == "forward"  & is_best       ~ colors_forward$best$border,
    node_phase == "backward" & is_best       ~ colors_backward$best$border,
    node_phase == "forward"                  ~ colors_forward$not_sig$border,
    TRUE                                     ~ colors_backward$not_sig$border
  )
  
  font_color <- case_when(
    is_base_fwd                              ~ colors_forward$base$font,
    is_base_bwd                              ~ colors_backward$base$font,
    node_phase == "forward"  & is_best       ~ colors_forward$best$font,
    node_phase == "backward" & is_best       ~ colors_backward$best$font,
    node_phase == "forward"  & is_sig        ~ colors_forward$sig$font,
    node_phase == "backward" & is_sig        ~ colors_backward$sig$font,
    node_phase == "forward"                  ~ colors_forward$not_sig$font,
    TRUE                                     ~ colors_backward$not_sig$font
  )
  
  # Arrow color based on destination node phase and status
  edge_colors <- map_chr(seq_len(igraph::ecount(g)), function(i) {
    to_uid  <- ends(g, i)[2]
    to_idx  <- match(to_uid, nodes$uid)
    ph      <- nodes$phase[to_idx]
    col_set <- if (ph == "forward") colors_forward else colors_backward
    if      (nodes$is_best[to_idx])      col_set$best$arrow
    else if (nodes$significant[to_idx])  col_set$sig$arrow
    else                                 col_set$not_sig$arrow
  })
  
  # Node labels
  node_labels  <- nodes$added_cov[idx]
  delta_labels <- ifelse(
    is.na(nodes$delta_OFV[idx]), "",
    sprintf("\u0394OFV = %.2f", nodes$delta_OFV[idx])
  )
  full_labels             <- paste0(node_labels, "\n", delta_labels)
  full_labels[is_base_fwd] <- forward_base_label
  full_labels[is_base_bwd] <- backward_base_label
  
  # Column headers — forward and backward each restart from "Iteration 1"
  col_headers <- nodes |>
    group_by(x_pos, phase, iter) |>
    summarise(.groups = "drop") |>
    mutate(label = case_when(
      phase == "forward"  & iter == min_fwd_iter ~ "Base",
      phase == "backward" & iter == 0            ~ "Full\nmodel",
      phase == "forward"                         ~ paste(forward_iter_label,  iter - min_fwd_iter),
      phase == "backward"                        ~ paste(backward_iter_label, iter)
    ))
  
  # Layout boundaries
  max_fwd_x   <- max(nodes$x_pos[nodes$phase == "forward"])
  min_bwd_x   <- min(nodes$x_pos[nodes$phase == "backward"])
  max_bwd_x   <- max(nodes$x_pos[nodes$phase == "backward"])
  sep_x       <- (max_fwd_x + min_bwd_x) / 2
  x_leg_fwd   <- min(nodes$x_pos[nodes$phase == "forward"])
  x_leg_bwd   <- min_bwd_x
  
  y_max       <- max(layout_df$y)
  y_min       <- min(layout_df$y)
  y_header    <- y_max + 1.2
  y_separator <- y_max + 0.85
  y_phase_lbl <- y_max + 1.7
  y_legend    <- y_min - 1.5
  cap_mm      <- node_size / 2
  
  p <- ggraph(g, layout = layout_mat) +
    
    geom_edge_link(
      aes(color = edge_colors),
      arrow      = arrow(length = unit(2.5, "mm"), type = "closed"),
      end_cap    = circle(cap_mm, "mm"),
      start_cap  = circle(cap_mm, "mm"),
      width      = 0.6,
      alpha      = 0.85
    ) +
    scale_edge_color_identity() +
    
    geom_node_point(
      shape  = 22,
      size   = node_size,
      fill   = node_fill,
      color  = node_border,
      stroke = 1.2
    ) +
    
    geom_node_text(
      aes(label = full_labels),
      size       = 2.5,
      lineheight = 1.15,
      color      = font_color
    ) +
    
    # Column headers
    annotate("text",
             x = col_headers$x_pos, y = y_header,
             label    = col_headers$label,
             fontface = "bold", size = 3.2, color = "gray20") +
    
    annotate("segment",
             x    = col_headers$x_pos - 0.45, xend = col_headers$x_pos + 0.45,
             y    = y_separator,              yend = y_separator,
             color = "gray60", linewidth = 0.4) +
    
    # Vertical dashed separator between phases
    annotate("segment",
             x = sep_x, xend = sep_x,
             y = y_min - 0.8, yend = y_phase_lbl - 0.1,
             color = separator_color, linewidth = 0.5, linetype = "dashed") +
    
    # Phase labels
    annotate("text",
             x        = (x_leg_fwd + max_fwd_x) / 2,
             y        = y_phase_lbl,
             label    = forward_title,
             fontface = "bold.italic", size = 3.5, color = "gray30") +
    
    annotate("text",
             x        = (min_bwd_x + max_bwd_x) / 2,
             y        = y_phase_lbl,
             label    = backward_title,
             fontface = "bold.italic", size = 3.5, color = "gray30") +
    
    # Forward legend
    annotate("point", shape = 22, size = 4,
             x = x_leg_fwd, y = y_legend,
             fill = colors_forward$best$fill, color = colors_forward$best$border) +
    annotate("text", hjust = 0, size = 2.6, color = "gray20",
             x = x_leg_fwd + 0.15, y = y_legend,
             label = forward_best_label) +
    
    annotate("point", shape = 22, size = 4,
             x = x_leg_fwd, y = y_legend - 0.4,
             fill = colors_forward$sig$fill, color = colors_forward$sig$border) +
    annotate("text", hjust = 0, size = 2.6, color = "gray20",
             x = x_leg_fwd + 0.15, y = y_legend - 0.4,
             label = paste0(forward_sig_label, " (LRT p < ", alpha_fwd, ")")) +
    
    annotate("point", shape = 22, size = 4,
             x = x_leg_fwd, y = y_legend - 0.8,
             fill = colors_forward$not_sig$fill, color = colors_forward$not_sig$border) +
    annotate("text", hjust = 0, size = 2.6, color = "gray20",
             x = x_leg_fwd + 0.15, y = y_legend - 0.8,
             label = paste0("Not ", tolower(forward_sig_label))) +
    
    # Backward legend
    annotate("point", shape = 22, size = 4,
             x = x_leg_bwd, y = y_legend,
             fill = colors_backward$best$fill, color = colors_backward$best$border) +
    annotate("text", hjust = 0, size = 2.6, color = "gray20",
             x = x_leg_bwd + 0.15, y = y_legend,
             label = backward_best_label) +
    
    annotate("point", shape = 22, size = 4,
             x = x_leg_bwd, y = y_legend - 0.4,
             fill = colors_backward$sig$fill, color = colors_backward$sig$border) +
    annotate("text", hjust = 0, size = 2.6, color = "gray20",
             x = x_leg_bwd + 0.15, y = y_legend - 0.4,
             label = paste0(backward_sig_label, " (LRT p < ", alpha_bwd, ")")) +
    
    annotate("point", shape = 22, size = 4,
             x = x_leg_bwd, y = y_legend - 0.8,
             fill = colors_backward$not_sig$fill, color = colors_backward$not_sig$border) +
    annotate("text", hjust = 0, size = 2.6, color = "gray20",
             x = x_leg_bwd + 0.15, y = y_legend - 0.8,
             label = paste0("Not ", tolower(backward_sig_label))) +
    
    labs(title = title) +
    
    theme_void() +
    theme(
      plot.title      = element_text(size = 14, face = "bold",
                                     color = "gray10", margin = margin(b = 15)),
      plot.margin     = margin(25, 25, 25, 25),
      plot.background = element_rect(fill = "gray98", color = NA)
    )
  
  p
}
