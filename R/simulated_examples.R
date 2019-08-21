##################################
# Australia
n_au <- 8
id_region_au <- c("NT", "QLD", "WA", "SA", 
                  "NSW", "ACT", "VIC", "TAS")
# source : https://itt.abs.gov.au/itt/r.jsp?databyregion
au_df <- data.frame(id = id_region_au,
                    x = c(20, 40, 7, 10, 30, 25, 15, 10),
                    wage = c(56783, 47177, 52691, 46619, 49256, 64901, 48878, 46800),
                    pop = c(247327, 5011216, 2595192, 1736422, 4988241, 420960, 4963349, 232606),
                    age = c(32.6, 37.1, 36.6, 40, 37.5, 35, 35.6, 39.8),
                    row.names = "id")

au_df$ln_wage <- log(au_df$wage) 
au_df$ln_pop <- log(au_df$pop) 
au_df$ln_age <- log(au_df$age) 

sp_au <- SpatialPointsDataFrame(cbind(c(1, 2, 0, 1, 2, 3, 2, 3),
                                      c(3, 3, 2, 2, 2, 2, 1, 0)),
                                au_df)
au_flows <- matrix(c(771843, 20075, 35107, 4578, 8058, 2117, 2822, 9667,
                     17523, 595782, 16172, 5023, 7485, 2573, 2370, 2440,
                     28953, 16657, 616281, 4363, 9348, 2963, 3834, 2727,
                     4292, 5907, 5363, 174941, 2983, 734, 1567, 764, 
                     6273, 7339, 7589, 2258, 286941, 1351, 1572, 699, 
                     1859, 2842, 2809, 628, 1590, 53427, 268, 295,
                     2832, 2619, 5253, 2664, 2469, 444, 23045, 478, 
                     8389, 2471, 2852, 548, 700, 221, 376, 36373), 8, 8,
                   dimnames = list(c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT"), 
                                   c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT"))
)

##################################
# Germany
id_region_ge <- c("SH", "HH", "MV", "NW", "HB", "BB", "BE", "RP", 
                  "NI", "ST", "SN", "SL", "HE", "TH", "BW", "BY")
ge_df <- data.frame(id = id_region_ge,
                    x = c(10, 15, 20, 7, 20, 25, 15, 10, 
                          30, 20, 15, 10, 15, 10, 7, 7), 
                    row.names = "id")
x_ge <- c(1, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 1, 2) 
y_ge <- c(5, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 0, 0)
sp_ge <- SpatialPointsDataFrame(cbind(x_ge, y_ge), ge_df)

##################################
# USA
id_region_usa <- c("AK", "ME", "WI", "VT", "NH", "WA", "ID", "MT", "ND", 
                   "MN", "IL", "MI", "NY", "MA", "OR", "NV", "WY", "SD", 
                   "IA", "IN", "OH", "PA", "NJ", "CT", "RI", "CA", "UT", 
                   "CO", "NE", "MO", "KY", "WV", "VA", "MD", "DE", "AZ", 
                   "NM", "KS", "AR", "TN", "NC", "SC", "DC", "OK", "LA",
                   "MS", "AL", "GA", "HI", "TX", "FL")
usa_df <- data.frame(id = id_region_usa,
                     x = c(35, 29, 30, 29, 26, 35, 31, 28, 32, 40, 32, 27, 33, 32,
                           31, 32, 25, 35, 32, 31, 35, 32, 38, 29, 35, 31, 27, 29,
                           30, 33, 38, 30, 30, 30, 40, 35, 34, 34, 33, 37, 31, 31, 
                           31, 27, 32, 31, 30, 28, 29, 29, 34), 
                     row.names = "id")
x_usa <- c(0, 10, 5, 9, 10, 0, 1, 2, 3, 4, 5, 6, 8, 9, 0:10, 0:9, 
           1:8, 3:7, 0, 3, 8)
y_usa <- c(7, 7, rep(6, 3), rep(5, 9), rep(4, 11), rep(3, 10),
           rep(2, 8), rep(1, 5), rep(0, 3))
sp_usa <- SpatialPointsDataFrame(cbind(x_usa, y_usa), usa_df)

##################################
# grid
# origin
x_grid_o <- rep(0:5, 5)
y_grid_o <- rep(4:0, each = 6)
id_region_grid_o <- paste0("O_", 1:(6 * 5))
grid_df_o <- data.frame(id = id_region_grid_o,
                     x = round(15 + sqrt(1 + 5 * x_grid_o^2 + 4 * y_grid_o^2), 0), 
                     row.names = "id")
sp_grid_o <- SpatialPointsDataFrame(cbind(x_grid_o, y_grid_o), grid_df_o)

# destination
x_grid_d <- rep(8:10, 4)
y_grid_d <- rep(4:1, each = 3)
id_region_grid_d <- paste0("D_", 1:(4 * 3))
grid_df_d <- data.frame(id = id_region_grid_d,
                        z = round(30 + sqrt(1 + 5 * x_grid_d^2 + 4 * y_grid_d^2), 0), 
                        row.names = "id")
sp_grid_d <- SpatialPointsDataFrame(cbind(x_grid_d, y_grid_d), grid_df_d)


create_grid <- function(sp_centroid) {
  # initialization
  # define the grid
  coords <- coordinates(sp_centroid)
  id <- row.names(sp_centroid)
  nx <- length(unique(coords[, 1])) 
  ny <- length(unique(coords[, 2]))  
  gt <- GridTopology(c(min(coords[, 1]), min(coords[, 2])), c(1, 1), c(nx, ny))
  grd <- SpatialGrid(gt)
  spix <- as(grd, "SpatialPixels")
  spol <- as(spix, "SpatialPolygons")
  o <- over(spol, sp_centroid)[, 1]
  simu_spdf <- spol[!is.na(o), ]
  row.names(simu_spdf) <- id
  simu_spdf <- SpatialPolygonsDataFrame(simu_spdf, sp_centroid@data) 
  return(simu_spdf)
}

spdf_au <- create_grid(sp_centroid = sp_au)
spdf_ge <- create_grid(sp_centroid = sp_ge)
spdf_usa <- create_grid(sp_centroid = sp_usa)
spdf_grid_o <- create_grid(sp_centroid = sp_grid_o)
spdf_grid_d <- create_grid(sp_centroid = sp_grid_d)

spdf_au$NOM <- row.names(spdf_au)
spdf_ge$NOM <- row.names(spdf_ge)
spdf_usa$NOM <- row.names(spdf_usa)
spdf_grid_o$NOM <- row.names(spdf_grid_o)
spdf_grid_d$NOM <- row.names(spdf_grid_d)


# # Big data 
# # Australia
# n_big <- 5000
# id_region_big <- as.character(1:n_big)
# # source : https://itt.abs.gov.au/itt/r.jsp?databyregion
# au_big <- data.frame(id = id_region_big,
#                     x1 = rnorm(n_big),
#                     x2 = rbinom(n_big, size = 1000, prob = 0.5),
#                     row.names = "id")
# flows_big_d <- sapply(au_big[, c("x1", "x2")], 
#                      function (x) kronecker(rep(1, n_big), x))
# colnames(flows_big_d) <- paste0(colnames(flows_big_d), "_d")
# 
# flows_big_o <- sapply(au_big[, c("x1", "x2")], 
#                      function (x) kronecker(x, rep(1, n_big)))
# colnames(flows_big_o) <- paste0(colnames(flows_big_o), "_o")
# 
# flows_big <- as.data.frame(cbind(flows_big_d, flows_big_o))
# y_big <- as.matrix(cbind(1, flows_big)) %*% c(1, 2, 3, 4, 5)
# flows_big$g <- log(1:n_big ^ 2)
# flows_big[, "y"] <- as.vector(y_big) - 10 * flows_big$g + rnorm(n_big ^ 2)
# 
# system.time(lm(y ~ x1_d + x2_d + x1_o + x2_o + g, data = flows_big))
# 
# Y_big <- matrix(flows_big$y, n_big, n_big)
# G_big <- matrix(flows_big$g, n_big, n_big)
# 
# system.time(
#   gravity_model(au_big, Y_big, G_big, centered = F)
# )
