# Australia
n_au <- 8
id_region_au <- c("NT", "QLD", "WA", "SA", 
                  "NSW", "ACT", "VIC", "TAS")
au_df <- data.frame(id = id_region_au,
                    x = c(20, 40, 7, 10, 30, 25, 15, 10),
                    row.names = "id")
sp_au <- SpatialPointsDataFrame(cbind(c(1, 2, 0, 1, 2, 3, 2, 3),
                                      c(3, 3, 2, 2, 2, 2, 1, 0)),
                                au_df)
# OD data
# source : https://www.abs.gov.au/websitedbs/censushome.nsf/home/factsheetsim?opendocument&navpos=450
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
# transform into data.frame
au_flows_df <- data.frame(origin = rep(dimnames(au_flows)[[1]], each = 8),
                          dest = rep(dimnames(au_flows)[[1]], 8),
                          Y = as.vector(au_flows))

