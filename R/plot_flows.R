plot_flows <- function(y, index_o, index_d, type_plot = "flow_map", 
                       xy_sf = NULL, contours_map = xy_sf, ordering = "distance",
                       sankey.options = list(
                         labels_od = c("origin", "destination"),
                         ylab = "",
                         title = "",
                         nrows = 3),
                       arc.options = list(
                         maxcex = 1,
                         maxlwd = 1,
                         alpha.q = 0.75
                       ),
                       flow_map.options = list(
                         maxlwd = 1,
                         alpha.q = 0.75,
                         max_bar = 1,
                         round.values = 0
                       ),
                       griffith.options = list(
                         round.label = 0,
                         scale_xy = c(-1/20, - 1),
                         maxsize = 1,
                         maxlwd = 1,
                         alpha.q = 0.75,
                         max_bar = 1.5,
                         title.legend = "Outflows/Inflows",
                         quantiles.legend = c(0.2, 0.4, 0.6, 0.8, 0.95)
                       ),
                       heatmap.options = list(
                         style_class = "kmeans",
                         n.class = 7,
                         title = "")) {
  
  # verification
  stopifnot(ordering %in% c("none", 
                            "longitude", 
                            "latitude", 
                            "descending",
                            "clustering",
                            "distance"))
  
  # Install packages if not already installed
  if (!require(colorspace)) {
    install.packages("colorspace")
  }
  if (!require(sf)) {
    install.packages("sf")
  }
  if (!require(tidyverse)) {
    install.packages("tidyverse")
  }
  
  ############### Initialisation   
  # number of flows 
  N <- length(y)
  # index of the origin
  O <- levels(unique(factor(index_o)))
  n_o <- length(O)
  D <- levels(unique(factor(index_d)))
  n_d <- length(D)
  # number of unique site in S
  S <- union(O, D)
  n <- length(S)
  # matricial form 
  Y <- matrix(0, n_o, n_d)
  Y_s <- matrix(0, n, n)
  dimnames(Y) <- list(O, D)
  dimnames(Y_s) <- list(S, S)
  for (k in 1:N) {
    Y[as.character(index_o)[k], as.character(index_d)[k]] <- y[k]
    Y_s[as.character(index_o)[k], as.character(index_d)[k]] <- y[k]
  }
  # outflows/inflows
  outflows <- apply(Y_s, 1, function(x) sum(x, na.rm = T))
  inflows <- apply(Y_s, 2, function(x) sum(x, na.rm = T))

  
  ############  Check on the spatial data 
  if (!is.null(xy_sf)) {
    stopifnot("S" %in% names(xy_sf), all(S %in% xy_sf$S))
  }
  
  # Order the sites S : Y or N
  if (ordering == "distance") {
    stopifnot(!is.null(xy_sf))
    mat_distance <- st_distance(xy_sf)
    my_index <- 1
    for(k in 1:n) {
      my_id <- 2
      my_order <- order(mat_distance[my_index[k], ])
      my_choice <- my_order[my_id]
      while(my_choice %in% my_index) {
        my_id <- my_id + 1
        my_choice <- my_order[my_id]
      }
      my_index <- c(my_index, my_choice)
    }
    order_S <- xy_sf$S[na.omit(my_index)] 
  } else {
    if (ordering == "longitude") {
      order_S <- xy_sf$S[order(st_coordinates(xy_sf)[, 1])]
    } else {
      if (ordering == "latitude") {
        order_S <- xy_sf$S[order(st_coordinates(xy_sf)[, 2])]
      } else {
        if (ordering == "descending") {
          order_O <- O[order(outflows[O])]
          order_D <- D[order(inflows[D])]
          order_S <- S[order(outflows + inflows)]
        } else {
          if (ordering == "clustering") {
            clust_origin <- hclust(dist(Y))
            clust_dest <- hclust(dist(t(Y)))
            clust_both <- hclust(dist(cbind(Y_s, t(Y_s))))
            order_O <- O[clust_origin$order]
            order_D <- D[clust_dest$order]
            order_S <- S[clust_both$order]
          } else {
            order_S <- S  
          }
        }
      }
    } 
  }
 
  # colors 
  q4 <- qualitative_hcl(n, palette = "Dark 3")
  names(q4) <- order_S
  
  ######  Sankey Diagram 
  if (type_plot == "Sankey") {
    # Install packages if not already installed
    if (!require(ggalluvial)) {
      install.packages("ggalluvial")
      library(ggalluvial)
    }
    
    # the optional argument
    if(is.null(sankey.options$labels_od)) 
      labels_od <- c("origin", "destination")
    else
      labels_od <- sankey.options$labels_od
    
    if(is.null(sankey.options$ylab)) 
      ylab <- ""
    else
      ylab <- sankey.options$ylab
    
    if(is.null(sankey.options$title)) 
      title <- ""
    else
      title <- sankey.options$title
  
    if(is.null(sankey.options$nrows)) 
      nrows <- 3
    else
      nrows <- sankey.options$nrows
    
    # preparation of the data
    data_long_2 <- data.frame(
      index_o = index_o,
      index_d = index_d,
      y = y)
    
    # create the name of the flow ("origin_destination")
    data_long_2$names_flow <- paste(data_long_2$index_o, data_long_2$index_d, sep = "_")
    # tidy the data
    data_long_2 <- pivot_longer(data_long_2, cols = c(1, 2),
                                names_to = "survey", values_to = "response")
    # we use factor
    data_long_2$survey <- factor(data_long_2$survey, levels = c("index_o", "index_d"),
                                 labels = labels_od)
    
    if (ordering %in% c("descending", "clustering")) {
      special_levels_1 <- order_O 
      data_long_2$names_flow <- paste0(data_long_2$names_flow, " ")
      boolScale <- scale_fill_manual(name = "Zone", values = rev(q4[order_O]))
    } else {
      special_levels_1 <- (order_S)
      boolScale <- scale_fill_manual(name = "Zone", values = rev(q4[order_S]))
      
    }
    
    data_long_2$zone <- factor(data_long_2$response, levels = special_levels_1)
    
    if (ordering %in% c("descending", "clustering")) {
      special_levels_2 <- c((order_O), (paste0(order_D, " ")))
      data_long_2$zone2 <- ifelse(data_long_2$survey == labels_od[1], 
                                  as.character(data_long_2$zone), paste0(data_long_2$zone, " "))
    } else {
      special_levels_2 <- order_S
      data_long_2$zone2 <- data_long_2$zone
    }
    
    data_long_2$zone2 <- factor(data_long_2$zone2, 
                                levels = special_levels_2)
    
    p <- ggplot(data_long_2,
                aes(x = survey, stratum = zone2, alluvium = names_flow,
                    y = y, fill = zone, label = zone)) +
      boolScale +
      scale_x_discrete(expand = c(.1, .1)) +
      geom_flow() +
      geom_stratum(alpha = .5) + 
      guides(fill = guide_legend(nrow = nrows, byrow = T)) +
      theme(legend.position = "bottom") +
      coord_flip() +
      xlab("") +
      ylab(ylab) +
      ggtitle(title)
    print(p)
  }
  
  if (type_plot == "circular") {
    if (!require(circlize)) {
      install.packages("circlize")
      library(circlize)
    }
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 60 / n)
    chordDiagram(x = Y, directional = 1, order = order_S, grid.col = q4, 
                 annotationTrack = "grid",
                 transparency = 0.25, annotationTrackHeight = c(0.1, 0.1),
                 preAllocateTracks = list(track.height = 0.1), diffHeight = -0.04)
    #add in labels and axis
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      sector.index = get.cell.meta.data("sector.index")
      # text direction (dd) and adjusmtents (aa)
      theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
      dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
      aa = c(1, 0.5)
      if(theta < 90 || theta > 270) aa =c(0, 0.5)
      circos.text(x = mean(xlim), y = 0.1, labels = sector.index, 
                  facing = dd, adj = aa)
    }, bg.border = NA)
  }
  
  if (type_plot == "arc") {
    if (!require(igraph)) {
      install.packages("igraph")
      library(igraph)
    }
    if (!require(arcdiagram)) {
      if (!require(devtools)) {
        install.packages("devtools")
      }
      devtools::install_github('gastonstat/arcdiagram')   # Arc diagram plot
      library(arcdiagram)
    }
    
    # modification of the arcdiagram function
    source("https://raw.githubusercontent.com/tibo31/spatial_flows/master/R/arcplot_b.R")
    
    # preparation of the options
    if(is.null(arc.options$maxcex)) 
      maxcex <- 15 
    else
      maxcex <- 15 * arc.options$maxcex
    
    if(is.null(arc.options$maxlwd))
      maxlwd <- 15 
    else
      maxlwd <- 15 * arc.options$maxlwd
    
    if(is.null(arc.options$alpha.q))
      alpha.q <- 0.75
    else
      alpha.q <- arc.options$alpha.q
    
    # create the adjacency matrix
    # create the graph
    nodes <- data.frame(S = order_S, outflows = outflows[order_S], inflows = inflows[order_S])
    relations <- data.frame(from = index_o, to = index_d, y = y)
    g <- graph_from_data_frame(relations, directed = TRUE, vertices = nodes)
    # create the adjacency matrix
    star_edges <- get.edgelist(g)
    
    # compute the size of each nodes
    lwd.nodes <- outflows[order_S] + inflows[order_S]
    lwd.nodes <- maxcex * lwd.nodes / max(lwd.nodes)
    # compute the width of the arcs
    lwd.arcs <- as.vector(maxlwd * y / max(y, na.rm = T))
    # compute the colors of the arcs
    col.arcs <- q4[star_edges[, 1]] 

    arcplot_b(star_edges[y > quantile(y, alpha.q),], 
            vertices = order_S,
            lwd.arcs = lwd.arcs[y > quantile(y, alpha.q)],
            col.arcs = col.arcs[y > quantile(y, alpha.q)], ordering = order_S,
            show.nodes = TRUE, pch.nodes = 21,
            col.nodes = q4[order_S], bg.nodes = q4[order_S], lwd.nodes = lwd.nodes,
            las = 2, font = 1, cex.labels = 0.9, col.labels = rgb(0.3, 0.3, 0.3))
  }
  
  # fonction qui créé un arc entre deux points 
  my_arc <- function(xA, yA, xB, yB) {
    xC <- (xA + xB)/2
    yC <- (yA + yB)/2
    R <- sqrt((yB - yA) ^ 2 + (xB - xA)^2)
    slope <- (yB - yA)/(xB - xA)
    s <- -1 / slope
    k <- sqrt(3) / 2 * R
    if (xB > xA) {
      xO <- xC - sqrt(k^2 / (s^2 + 1))
      yO <- yC - k * s / sqrt(s^2 + 1)
    } else {
      xO <- xC + sqrt(k^2 / (s^2 + 1))
      yO <- yC + k * s / sqrt(s^2 + 1)     
    }
    my_x <- seq(xA, xB, length.out = 100)
    if (yO < min(yA,yB)) {
      my_y <- yO + sqrt(R^2 - (my_x - xO)^2)
    } else {
      if (yO > max(yA, yB)) {
        my_y <- yO - sqrt(R^2 - (my_x - xO)^2)
      } else {
        if(xA > xB) {
          my_x_1 <- seq(xA, xO - R, length.out = 100)
          my_x_2 <- seq(xO - R, xB, length.out = 100)
        } else {
          my_x_1 <- seq(xA, xO + R, length.out = 100)
          my_x_2 <- seq(xO + R, xB, length.out = 100) 
        }
        
        cond_positiv <- (R^2 - (my_x_1 - xO)^2 > 0) & (R^2 - (my_x_2 - xO)^2 > 0)
        if(yA > yO) {
          my_y_1 <- yO + sqrt(R^2 - (my_x_1[cond_positiv] - xO)^2)  
          my_y_2 <- yO - sqrt(R^2 - (my_x_2[cond_positiv] - xO)^2)
        } else { 
          my_y_1 <- yO - sqrt(R^2 - (my_x_1[cond_positiv] - xO)^2)
          my_y_2 <- yO + sqrt(R^2 - (my_x_2[cond_positiv] - xO)^2)
        }
        my_y<-c(my_y_1, my_y_2)
        my_x<-c(my_x_1[cond_positiv], my_x_2[cond_positiv])
      }
    }
    my_na <- is.na(my_y)
    return(cbind(my_x[!my_na], my_y[!my_na]))
  }
  
  if (type_plot == "flow_map") {
    
    # initialisation 
    if (is.null(flow_map.options$maxlwd))
      maxlwd <- 1 * 7
    else
      maxlwd <- flow_map.options$maxlwd * 7
      
    if (is.null(flow_map.options$alpha.q))
      alpha.q <- 0.85
    else
      alpha.q <- flow_map.options$alpha.q
    
    if (is.null(flow_map.options$max_bar))
      max_bar <- 1 * 1.8
    else
      max_bar <- flow_map.options$max_bar * 1.8
    
    if (is.null(flow_map.options$round.values))
      round.values <- 0
    else
      round.values <- flow_map.options$round.values
    
    # width of the flows
    maxlwd <- maxlwd * y / max(y, na.rm = T)
    # vector of colors for the flows
    my_col_flow <- q4[as.numeric(factor(index_o, levels = order_S))] 
    my_col_bar <- q4[as.numeric(factor(O, levels = order_S))]
    
    # parameter to shift the destination barplot
    shift <- 1 / 50
    
    # define the coordinates of origin site: it corresponds to the coordinates
    # of s in S, slightly shifted
    shift_coords_x <- diff(range(st_coordinates(xy_sf)[, 1])) * shift / 2
    shift_coords_y <- diff(range(st_coordinates(xy_sf)[, 2])) * shift / 12
    
    xy_coord <- st_coordinates(xy_sf)
    rownames(xy_coord) <- xy_sf$S
    
    xy_origin <- xy_coord[O, ]
    xy_origin[, 1] <- xy_origin[, 1] - shift_coords_x
    xy_origin[, 2] <- xy_origin[, 2] - shift_coords_y
    
    xy_dest <- xy_coord[D, ]
    xy_dest[, 1] <- xy_dest[, 1] + shift_coords_x
    xy_dest[, 2] <- xy_dest[, 2] + shift_coords_y
    
    # shift of the two bar
    shift_bar <- diff(range(xy_origin[, 1])) * shift / 3
    
    # maximum height for the bars 
    max_bar <- shift_bar * 6 * max_bar
    bar_out <-  max_bar * outflows[O] / max(c(inflows, outflows), na.rm = T) 
    bar_in <- max_bar * inflows[D] / max(c(inflows, outflows), na.rm = T) 
    
    # pdf("figures/map.pdf", width = 8, height = 8)
    # par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    plot(st_geometry(contours_map), border = "lightgrey", lwd = 0.5)
    
    # vectorial version 
    # plot the highest flows
    ind_biggest <- which(y > quantile(y, alpha.q))
    for(i in ind_biggest) {
      A <- xy_origin[as.character(index_o[i]), ]
      B <- xy_dest[as.character(index_d[i]), ]
      xA <- A[1]
      yA <- A[2]
      xB <- B[1]
      yB <- B[2]
      my_arc_don <- my_arc(xA, yA, xB, yB)
      lines(my_arc_don[, 1], my_arc_don[, 2], 
            lwd = maxlwd[i], col = my_col_flow[i])
    }
    
    bar_1x <- cbind(xy_origin[,1] - shift_bar, xy_origin[,1] + shift_bar, 
                    xy_origin[,1] + shift_bar, xy_origin[,1] - shift_bar,
                    xy_origin[,1] - shift_bar)
    bar_2x <- cbind(xy_dest[,1] - shift_bar, xy_dest[,1] + shift_bar,
                    xy_dest[,1] + shift_bar, xy_dest[,1] - shift_bar,
                    xy_dest[,1] - shift_bar)
    
    # outflows bar
    for(k in 1:n_o) {
      bar_1y <- cbind(xy_origin[k, 2], xy_origin[k, 2], 
                      xy_origin[k, 2] + bar_out[k], 
                      xy_origin[k, 2] + bar_out[k],
                      xy_origin[k, 2])
      polygon(bar_1x[k, ], bar_1y, col = my_col_bar[k])
    }
    
    # redefine the y_coordinates of the barplot of the inflows origin = destination
    xy_dest[, 2] <- xy_dest[, 2] - 2 * shift_coords_y
    
    # inflows bar
    for(k in 1:n_d) {
      local_inflow <- numeric(n_o)
      for (i in 1:n_o) {
        does_ind_exist <- which(index_o == O[i] & index_d == D[k])
        if (length(does_ind_exist) == 1)
          local_inflow[i] <- y[does_ind_exist]
      }
      my_cum_sum <- xy_dest[k, 2] + c(0, bar_in[k] * cumsum(local_inflow) / sum(local_inflow)) 
      
      for (j in 1:n_o) {
        bar_2y <- cbind(my_cum_sum[j], my_cum_sum[j], 
                        my_cum_sum[j+1], my_cum_sum[j+1],
                        my_cum_sum[j])
        
        polygon(bar_2x[k,], bar_2y, col = my_col_bar[j], border = my_col_bar[j])
      }
      polygon(bar_2x[k, ], cbind(my_cum_sum[1], my_cum_sum[1], 
                                 my_cum_sum[n_o + 1], my_cum_sum[n_o + 1],
                                 my_cum_sum[1]))
    }
    
    # plot the legend
    # max of the outflows
    polygon(cbind(par()$xaxp[1] - shift_bar, 
                  par()$xaxp[1] + shift_bar, 
                  par()$xaxp[1] + shift_bar, 
                  par()$xaxp[1] - shift_bar,
                  par()$xaxp[1] - shift_bar),
            cbind(par()$yaxp[1], par()$yaxp[1], 
                  par()$yaxp[1] + max(bar_out), 
                  par()$yaxp[1] + max(bar_out),
                  par()$yaxp[1]))
    
    # max of the inflows
    polygon(cbind(par()$xaxp[1] - shift_bar + diff(range(xy_origin[, 1])) * shift, 
                  par()$xaxp[1] + shift_bar + diff(range(xy_origin[, 1])) * shift,
                  par()$xaxp[1] + shift_bar + diff(range(xy_origin[, 1])) * shift, 
                  par()$xaxp[1] - shift_bar + diff(range(xy_origin[, 1])) * shift,
                  par()$xaxp[1] - shift_bar + diff(range(xy_origin[, 1])) * shift),
            cbind(par()$yaxp[1], par()$yaxp[1], 
                  par()$yaxp[1] + max(bar_in), 
                  par()$yaxp[1] + max(bar_in),
                  par()$yaxp[1]))
    
    # Print out and In
    text(par()$xaxp[1], par()$yaxp[1], "Out", cex = 0.6, pos = 1)
    text(par()$xaxp[1] + diff(range(xy_origin[, 1])) * shift, par()$yaxp[1], "In", cex = 0.6, pos = 1)
    
    # Print the arrows 
    arrows(par()$xaxp[1] + 2 * diff(range(xy_origin[, 1])) * shift,   
           par()$yaxp[1], par()$xaxp[1] + 2 * diff(range(xy_origin[, 1])) * shift, 
           par()$yaxp[1] + max(bar_out, bar_in),
           length = 0.1)
    
    # Print the values  
    text(par()$xaxp[1] + 2.1 * diff(range(xy_origin[, 1])) * shift,
         par()$yaxp[1], round.values, cex = 0.5, pos = 4)
    text(par()$xaxp[1] + 2.1 * diff(range(xy_origin[, 1])) * shift, 
         par()$yaxp[1] + max(bar_out, bar_in), 
         round(max(outflows, inflows), round.values), 
         cex = 0.5, pos = 4)
    text(par()$xaxp[1] + 2.1 * diff(range(xy_origin[, 1])) * shift,
         par()$yaxp[1] + max(bar_out, bar_in) / 2, 
         round(max(outflows, inflows) / 2, round.values), 
         cex = 0.5, pos = 4)
    legend(par()$xaxp[2], par()$yaxp[2],
           legend = round(quantile(y, c(0.5, 0.95, 0.99, 0.999)), round.values)[-1], 
           lty = 1,
           lwd = quantile(maxlwd, c(0.5, 0.95, 0.99, 0.999))[-1],
           cex = 0.6)
    # dev.off()
    
  }
  
  if (type_plot == "griffith") {
 
     
    # round the values in the legend
    if(is.null(griffith.options$round.label))
      round.label <- 0
    else
      round.label <- griffith.options$round.label
    
    # shift origin
    if(is.null(griffith.options$scale_xy))
      scale_xy <- c(-1/20, - 1)
    else
      scale_xy <- griffith.options$scale_xy
    
    # maxsize of the bubbles
    if(is.null(griffith.options$maxsize))
      maxsize <- 5
    else
      maxsize <- griffith.options$maxsize * 5
    maxsize <- maxsize / max(sqrt(c(inflows, outflows)), na.rm = T)
    
    # maximum size of the width
    if(is.null(griffith.options$maxlwd))
      maxlwd <- 10
    else
      maxlwd <- griffith.options$maxlwd * 10

    # values of griffith.options quantile for representing the biggest flows
    if(is.null(griffith.options$alpha.q))
      alpha.q <- 0.75
    else
      alpha.q <- griffith.options$alpha.q
    
    # height of the bar
    if(is.null(griffith.options$max_bar))
      max_bar <- 1.5
    else
      max_bar <- griffith.options$max_bar
    
    # legend 
    if(is.null(griffith.options$title.legend))
      title.legend <- "Outflows/Inflows"
    else
      title.legend <- griffith.options$title.legend
    
    if(is.null(griffith.options$quantiles.legend))
      quantiles.legend <- c(0.2, 0.4, 0.6, 0.8, 0.95)
    else
      quantiles.legend <- griffith.options$quantiles.legend
    
    # width of the flows
    maxlwd <- maxlwd * y / max(y, na.rm = T)
    # vector of colors for the flows
    my_col_flow <- q4[as.numeric(as.factor(index_o))] 
    my_col_bar <- q4[as.numeric(as.factor(O))]
    # create two contours : one for origin, another for destination
    poly_sf_o <- contours_map
    xy_sf_o <- xy_sf
    bbox <- st_bbox(contours_map)
    poly_sf_d <- st_geometry(contours_map) + c((bbox[3] - bbox[1]) * scale_xy[1], 
                                               (bbox[4] - bbox[2]) * scale_xy[2])
    xy_sf_d <- st_geometry(xy_sf) + c((bbox[3] - bbox[1]) * scale_xy[1], 
                                      (bbox[4] - bbox[2]) * scale_xy[2])
    
    plot(st_geometry(poly_sf_o),  
         xlim = range(c(st_bbox(poly_sf_o)[c(1, 3)], st_bbox(poly_sf_d)[c(1, 3)])),
         ylim = range(c(st_bbox(poly_sf_o)[c(2, 4)], st_bbox(poly_sf_d)[c(2, 4)])))
    plot(poly_sf_d, add = T, border = "lightgrey")
    
    # new coordinates of the origin sites
    xy_coord <- st_coordinates(xy_sf_o)
    rownames(xy_coord) <- xy_sf$S
    xy_origin <- xy_coord[O, ]
    # new coordinates of the destination sites
    xy_coord <- st_coordinates(xy_sf_d)
    rownames(xy_coord) <- xy_sf$S
    xy_dest <- xy_coord[D, ]
    # print the bubbles
    points(xy_origin, 
           cex = maxsize * sqrt(abs(outflows[O])), 
           pch = 16, 
           col = q4[O])
    
    points(xy_dest, 
           cex = maxsize * sqrt(abs(inflows[D])),  
           pch = 16, 
           col = q4[D])
    
    # vectorial version 
    # plot the highest flows
    ind_biggest <- which(y > quantile(y, alpha.q))
    
    for(i in ind_biggest) {
      A <- xy_origin[as.character(index_o[i]), ]
      B <- xy_dest[as.character(index_d[i]), ]
      xA <- A[1]
      yA <- A[2]
      xB <- B[1]
      yB <- B[2]
      my_arc_don <- my_arc(xA, yA, xB, yB)
      lines(my_arc_don[, 1], my_arc_don[, 2], 
            lwd = maxlwd[i], col = my_col_flow[i])
    }
    
    
    legend("topleft", legend = round(quantile(c(outflows, inflows), quantiles.legend), 
                                     round.label), 
           pch = 16,
           pt.cex = maxsize * sqrt(quantile(c(outflows, inflows), quantiles.legend)),
           title = title.legend)
    legend("bottomleft",
           legend = round(quantile(y, c(0.5, 0.95, 0.99, 0.999)), round.label)[-1], 
           lty = 1,
           lwd = quantile(maxlwd, c(0.5, 0.95, 0.99, 0.999))[-1],
           cex = 0.7, title = "Flows")
  }
  
  
  if (type_plot == "heatmap") {
    
    if (!require(classInt)) {
      install.packages("classInt")
      library(classInt)
    }
    if (!require(gplots)) {
      install.packages("gplots")
      library(gplots)
    } 
    if (!require(viridis)) {
      install.packages("viridis")
      library(viridis)
    }
    
    gplots::heatmap.2(Y, col = rev(magma(7)), 
                      breaks = classInt::classIntervals(as.numeric(Y), 
                                                        style  = heatmap.options$style_class, 
                                                        n = heatmap.options$n.class)$brks,
                      key = FALSE,
                      srtRow = -50,     # angle
                      srtCol = 30,      # angle
                      trace= "none",
                      main = heatmap.options$title)
  }
  
}

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