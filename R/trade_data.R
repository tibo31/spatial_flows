# import data
wto_data <- read_dta("data/Chapter1Application1.dta")

# Filter year 2006
wto <- wto_data %>% 
  filter(year %in% 2006) 

wto[wto$exporter == wto$importer, "DIST"] <- 0

# data frame of spatial units
wto_country <- unique(wto[, "exporter"])
names(wto_country) <- "country"
wto_country$GDP <- c(637486, 1408675, 416835, 494763, 58222, 37508, 2055512,
                     1647120, 678938, 277080, 12234781, 34924, 309191, 57564,
                     22054, 3693204, 329865, 104295, 195135, 1314314, 252246,
                     2582492,  2631228, 203085, 341659, 139761, 1015539, 2575666,
                     331430, 460976, 24476, 353268, 1943835, 40708, 4872415,
                     74938, 1577524, 119534, 87356, 50361, 109708, 1158229,
                     12553, 67101, 13366, 6339, 314707, 8119, 375769, 
                     830572, 399470, 24870, 61838, 313595, 526211, 195041,
                     167605, 211803, 21126, 323901, 535607, 455302, 22104,
                     40068, 851541, 52090, 59180, 19485394, 348872)

# spatial polygons
# We now want to associate the spatial contours to the countries. 
# We first deal with Hong-Kong which is not included in our data basis. 
# For creating it, we decide to associate the same shape than Macau and transfer slightly the polygon. 

world <- readOGR("./data/World WGS84/Pays_WGS84.shp")
macau <- world[grep("Macau", world@data$NOM), ] 
crds <- macau@polygons[[1]]@Polygons[[1]]@coords + 0.04
Pl <- Polygon(crds)
# str(Pl)
ID <- "HKG"
Pls <- Polygons(list(Pl), ID=ID)
# str(Pls)
SPls <- SpatialPolygons(list(Pls))
# str(SPls)
df <- data.frame(NOM = "HKG", row.names=ID)
# str(df)
hkg <- SpatialPolygonsDataFrame(SPls, df)
proj4string(hkg) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

wto_spatial <- world[c(grep("Arg", world@data$NOM),
                       grep("Aus", world@data$NOM)[1],
                       grep("Aus", world@data$NOM)[2], 
                       grep("Bel", world@data$NOM)[1],
                       grep("Bul", world@data$NOM)[1],
                       grep("Bol", world@data$NOM), 
                       grep("Bra", world@data$NOM), 
                       grep("Can", world@data$NOM), 
                       grep("Swi", world@data$NOM),
                       grep("Chi", world@data$NOM)[2], 
                       grep("Chi", world@data$NOM)[1],
                       grep("Cam", world@data$NOM)[2], 
                       grep("Col", world@data$NOM), 
                       grep("Cos", world@data$NOM), 
                       grep("Cyp", world@data$NOM), 
                       grep("Ger", world@data$NOM), 
                       grep("Den", world@data$NOM), 
                       grep("Ecu", world@data$NOM),
                       grep("Egy", world@data$NOM), 
                       grep("Spa", world@data$NOM),
                       grep("Fin", world@data$NOM),
                       grep("Fra", world@data$NOM), 
                       grep("Kin", world@data$NOM), 
                       grep("Gre", world@data$NOM)[2]), ]

wto_spatial <- spRbind(wto_spatial, hkg)          
wto_spatial <- spRbind(wto_spatial, world[c(
  grep("Hun", world@data$NOM),
  grep("Indo", world@data$NOM),
  grep("India", world@data$NOM)[1],
  grep("Ire", world@data$NOM),
  grep("Iran", world@data$NOM),
  grep("Ice", world@data$NOM),
  grep("Isr", world@data$NOM), 
  grep("Ita", world@data$NOM), 
  grep("Jor", world@data$NOM), 
  grep("Jap", world@data$NOM), 
  grep("Ken", world@data$NOM),
  grep("Kor", world@data$NOM)[2],
  grep("Kuw", world@data$NOM),
  grep("Sri", world@data$NOM), 
  grep("Mac", world@data$NOM)[2],
  grep("Mor", world@data$NOM),
  grep("Mex", world@data$NOM),
  grep("Malta", world@data$NOM),
  grep("Mya", world@data$NOM),
  grep("Mau", world@data$NOM)[2], 
  grep("Malawi", world@data$NOM), 
  grep("Malay", world@data$NOM), 
  grep("Niger", world@data$NOM)[1],
  grep("Niger", world@data$NOM)[2], 
  grep("Net", world@data$NOM)[1], 
  grep("Nor", world@data$NOM)[1],
  grep("Nep", world@data$NOM), 
  grep("Pan", world@data$NOM), 
  grep("Phi", world@data$NOM), 
  grep("Pol", world@data$NOM)[1], 
  grep("Por", world@data$NOM),
  grep("Qat", world@data$NOM), 
  grep("Rom", world@data$NOM), 
  grep("Sen", world@data$NOM),
  grep("Sin", world@data$NOM), 
  grep("Swe", world@data$NOM), 
  grep("Tha", world@data$NOM),
  grep("Tri", world@data$NOM), 
  grep("Tun", world@data$NOM), 
  grep("Tur", world@data$NOM)[2],
  grep("Tan", world@data$NOM), 
  grep("Uru", world@data$NOM),
  grep("Uni", world@data$NOM)[3], 
  grep("Sou", world@data$NOM)[2]), ])
wto_spatial$NOM <- unique(wto$exporter)
wto_spatial$GDP <- wto_country$GDP 
