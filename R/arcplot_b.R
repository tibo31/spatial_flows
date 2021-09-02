# modification of the function arcdiagram (Sanchez, 2018)
# to be adapted for a matrix of flows 
arcplot_b <- function (edgelist, vertices, sorted = FALSE, decreasing = FALSE, 
          ordering = NULL, labels = NULL, horizontal = TRUE, above = NULL, 
          col.arcs = "#5998ff77", lwd.arcs = 1.8, lty.arcs = 1, lend = 1, 
          ljoin = 2, lmitre = 1, show.nodes = TRUE, pch.nodes = 19, 
          cex.nodes = 1, col.nodes = "gray80", bg.nodes = "gray80", 
          lwd.nodes = 1, show.labels = TRUE, col.labels = "gray55", 
          cex.labels = 0.9, las = 2, font = 1, line = 0, outer = FALSE, 
          adj = NA, padj = NA, axes = FALSE, xlim = NULL, ylim = NULL, 
          ...) {
  if (hasArg(vertices)) {
    nodes_edges = graph_info(edgelist, vertices = vertices, 
                             sorted = sorted, decreasing = decreasing, ordering = ordering, 
                             labels = labels)
  }
  else {
    nodes_edges = graph_info(edgelist, sorted = sorted, decreasing = decreasing, 
                             ordering = ordering, labels = labels)
  }
  nodes = nodes_edges$nodes
  num_nodes = nodes_edges$num_nodes
  num_edges = nodes_edges$num_edges
  aux_ord = nodes_edges$aux_ord
  labels = nodes_edges$labels
  aux_ord <- 1:num_nodes
  centers = xynodes(num_nodes, aux_ord, labels)
  above = above_below(edgelist, above)
  radios_locs = arc_radius_locs(edgelist, nodes, centers)
  radios = radios_locs$radios
  locs = radios_locs$locs
  if (length(col.arcs) != num_edges) 
    col.arcs = rep(col.arcs, length = num_edges)
  if (length(lwd.arcs) != num_edges) 
    lwd.arcs = rep(lwd.arcs, length = num_edges)
  if (length(lty.arcs) != num_edges) 
    lty.arcs = rep(lty.arcs, length = num_edges)
  if (length(pch.nodes) != num_nodes) {
    pch.nodes = rep(pch.nodes, length = num_nodes)
  }
  pch.nodes = pch.nodes[aux_ord]
  if (length(cex.nodes) != num_nodes) {
    cex.nodes = rep(cex.nodes, length = num_nodes)
  }
  cex.nodes = cex.nodes[aux_ord]
  if (length(col.nodes) != num_nodes) {
    col.nodes = rep(col.nodes, length = num_nodes)
  }
  col.nodes = col.nodes[aux_ord]
  if (length(bg.nodes) != num_nodes) {
    bg.nodes = rep(bg.nodes, length = num_nodes)
  }
  bg.nodes = bg.nodes[aux_ord]
  if (length(lwd.nodes) != num_nodes) {
    lwd.nodes = rep(lwd.nodes, length = num_nodes)
  }
  lwd.nodes = lwd.nodes[aux_ord]
  if (length(col.labels) != num_nodes) {
    col.labels = rep(col.labels, length = num_nodes)
  }
  col.labels = col.labels[aux_ord]
  if (length(cex.labels) != num_nodes) {
    cex.labels = rep(cex.labels, length = num_nodes)
  }
  cex.labels = cex.labels[aux_ord]
  z = seq(0, pi, length.out = 100)
  if (horizontal) {
    side = 1
  }
  else {
    side = 2
  }
  if (is.null(xlim)) {
    if (horizontal) {
      xlim = c(-0.015, 1.015)
      x_nodes = centers
    }
    else {
      xlims = min_max_margin(radios, above)
      xlim = c(xlims$min, xlims$max)
      x_nodes = rep(0, num_nodes)
    }
  }
  else {
    if (horizontal) {
      x_nodes = centers
    }
    else {
      x_nodes = rep(0, num_nodes)
    }
  }
  if (is.null(ylim)) {
    if (horizontal) {
      ylims = min_max_margin(radios, above)
      ylim = c(ylims$min, ylims$max)
      y_nodes = rep(0, num_nodes)
    }
    else {
      ylim = c(-0.015, 1.015)
      y_nodes = centers
    }
  }
  else {
    if (horizontal) {
      y_nodes = rep(0, num_nodes)
    }
    else {
      y_nodes = centers
    }
  }
  plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n", xlab = "", 
       ylab = "", axes = axes, ...)
  for (i in 1L:num_edges) {
    radio = radios[i]
    if (horizontal) {
      x_arc = locs[i] + radio * cos(z)
      if (above[i]) {
        y_arc = radio * sin(z)
      }
      else {
        y_arc = radio * sin(-z)
      }
    }
    else {
      y_arc = locs[i] + radio * cos(z)
      if (above[i]) {
        x_arc = radio * sin(z)
      }
      else {
        x_arc = radio * sin(-z)
      }
    }
    lines(x_arc, y_arc, col = col.arcs[i], lwd = lwd.arcs[i], 
          lty = lty.arcs[i], lend = lend, ljoin = ljoin, lmitre = lmitre)
    if (show.nodes) {
      points(x = x_nodes, y = y_nodes, pch = pch.nodes, 
             col = col.nodes, bg = bg.nodes, cex = cex.nodes, 
             lwd = lwd.nodes)
    }
    if (show.labels) {
      #mtext(labels, side = side, line = line, at = centers, 
      #      cex = cex.labels, outer = outer, col = col.labels, 
      #      las = las, font = font, adj = adj, padj = padj, 
      #      ...)
      text(centers, -0.04, labels, pos = 2, cex = cex.labels, srt = 45,       # rotation des étiquettes
           xpd = T)
    }
    # plot the intra flow
    if (edgelist[i, 1] == edgelist[i, 2]) {
      R <- diff(centers)[1] / 5 
      x_0 <- c(0, -0.01, -0.015, -0.017, -0.018, -0.0185, -0.0185, -0.018, -0.017, -0.015, -0.01, -0.005, 
             0.005, .01, 0.015, 0.017, 0.018, 0.0185, 0.0185, 0.018, 0.017, 0.015, 0.01, 0)
      x_0 <- centers[edgelist[i, 1]] + x_0 * R / max(x_0)
      y_0 <- c(seq(0, -0.04, length.out = 12), seq(-0.04, 0, length.out = 12))
      lines(x_0, y_0, xpd = T, col = col.arcs[i], lwd = lwd.arcs[i], 
            lty = lty.arcs[i])
    }
  }
}