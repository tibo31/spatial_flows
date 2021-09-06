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
                         round.values = 0,
                         col_geometry = "lightgrey" , 
                         col_border = "white"
                       ),
                       griffith.options = list(
                         round.label = 0,
                         scale_xy = c(-1/20, - 1),
                         maxsize = 1,
                         maxlwd = 1,
                         alpha.q = 0.75,
                         max_bar = 1.5,
                         title.legend = "Outflows/Inflows",
                         quantiles.legend = c(0.2, 0.4, 0.6, 0.8, 0.95),
                         col_geometry_o = "lightgrey" , 
                         col_border_o = "white",
                         col_geometry_d = "lightgrey" , 
                         col_border_d = "white"
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
    order_O <- order_S[order_S %in% O]
    order_D <- order_S[order_S %in% D]
  } else {
    if (ordering == "longitude") {
      order_S <- xy_sf$S[order(st_coordinates(xy_sf)[, 1])]
      order_O <- order_S[order_S %in% O]
      order_D <- order_S[order_S %in% D]
    } else {
      if (ordering == "latitude") {
        order_S <- xy_sf$S[order(st_coordinates(xy_sf)[, 2])]
        order_O <- order_S[order_S %in% O]
        order_D <- order_S[order_S %in% D]
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
            order_O <- O
            order_D <- D
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
    
    data_long_2$zone <- factor(data_long_2$response, levels = rev(special_levels_1))
    
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
    
    if (is.null(flow_map.options$col_geometry))
      col_geometry <- "lightgrey" 
    else
      col_geometry <- flow_map.options$col_geometry
    
    if (is.null(flow_map.options$col_border))
      col_border <- "white" 
    else
      col_border <- flow_map.options$col_border    
    
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
    plot(st_geometry(contours_map), 
         col = col_geometry,
         border = col_border, lwd = 0.2)
    
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
    
    if (is.null(griffith.options$col_geometry_o))
      col_geometry_o <- "lightgrey" 
    else
      col_geometry_o <- griffith.options$col_geometry_o
    
    if (is.null(griffith.options$col_border_o))
      col_border_o <- "white" 
    else
      col_border_o <- griffith.options$col_border_o   
    
    if (is.null(griffith.options$col_geometry_d))
      col_geometry_d  <- "lightgrey" 
    else
      col_geometry_d  <- griffith.options$col_geometry_d 
    
    if (is.null(griffith.options$col_border_d))
      col_border_d  <- "white" 
    else
      col_border_d  <- griffith.options$col_border_d  
    
    # width of the flows
    maxlwd <- maxlwd * y / max(y, na.rm = T)
    
    # vector of colors for the flows
    my_col_flow <- q4[as.character(index_o)] 

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
         ylim = range(c(st_bbox(poly_sf_o)[c(2, 4)], st_bbox(poly_sf_d)[c(2, 4)])),
         col = col_geometry_o, border = col_border_o)
    plot(poly_sf_d, add = T, col = col_geometry_d, border = col_border_d)
    
    # new coordinates of the origin sites
    xy_coord <- st_coordinates(xy_sf_o)
    rownames(xy_coord) <- xy_sf$S
    xy_origin <- xy_coord[O, ]
    # new coordinates of the destination sites
    xy_coord <- st_coordinates(xy_sf_d)
    rownames(xy_coord) <- xy_sf$S
    xy_dest <- xy_coord[D, ]
    
    # print the bubbles
    points(xy_origin[O, ], 
           cex = maxsize * sqrt(abs(outflows[O])), 
           pch = 16, 
           col = q4[O])
    
    points(xy_dest[D, ], 
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
    if (!require(plot.matrix)) {
      install.packages("plot.matrix")
      library(plot.matrix)
    }
    
    
    if(is.null(heatmap.options$title))
      title <- ""
    else
      title <- heatmap.options$title
    
    if (is.null(heatmap.options$n.class))
      n.class <- 7
    else
      n.class <- heatmap.options$n.class
    
    if (is.null(heatmap.options$style_class))
      style_class <- "kmeans" 
    else
      style_class <- heatmap.options$style_class   
    

    
    if (ordering == "clustering") {
      gplots::heatmap.2(Y, col = rev(magma(7)), 
                      breaks = classInt::classIntervals(as.numeric(Y), 
                                                        style  = style_class, 
                                                        n = n.class)$brks,
                      key = FALSE,
                      srtRow = -50,     # angle
                      srtCol = 30,      # angle
                      trace= "none",
                      main = title)
    } else {
      Y_permut_2 <- Y[order_O, order_D]
      plot(Y_permut_2,
           na.cell = F,
           axis.col = list(cex.axis = 0.01), axis.row = list(cex.axis = 0.85),
           breaks = classInt::classIntervals(as.numeric(Y_permut_2), 
                                             style  = style_class, n = n.class)$brks, 
           col = rev(magma(7)), main = title, 
           xlab = "",
           ylab = "")
      
      text(seq(1, length(D), by = 1), par("usr")[3] - 1, 
           srt = 60, adj = 1, xpd = TRUE,
           labels = colnames(Y_permut_2), cex = 0.85)
    }
  }
  
}
