# source :
# https://journals.sagepub.com/doi/suppl/10.1177/0022002718823907
require(haven)
migration <- read_dta("./data/Migration_data_JCR.dta.dta")
mig_data <- migration %>%
  filter(year == 2014) %>%
  dplyr::select(iso3_d, iso3_o, country_d, country_o, 8, 29, 30, 34, 35, 41) %>%
  arrange(iso3_d, iso3_o)

mig_data$migi[is.na(mig_data$migi)] <- 0

# replace NA by zero
# number of origin and destination
n_o_mi <- length(unique(migration$iso3_o))
n_d_mi <- length(unique(migration$iso3_d))

# construct the variable for origin
mig_data_origin <- unique(mig_data[, c("iso3_o", "country_o", "popul_o", 
                                       "gdp_o", "civilconflict_o")])

mig_data_dest <- unique(mig_data[, c("iso3_d", "country_d", 
                                     "popul_d", "gdp_d")])
# correct some mistakes
mig_data_dest <- mig_data_dest[-c(4, 9, 13, 17, 19, 21, 23, 25,
                                  27, 29, 31, 33), ]  %>%
  arrange(country_d)
mig_data_origin <- mig_data_origin[-46, ]  %>%
  arrange(country_o)

# create the spatials
require("sf")
# https://gadm.org/download_country_v3.html
world <- st_read("./data/World WGS84/Pays_WGS84.shp")
# correction of morocco
moroco <- st_union(world[world$NOM %in% c("Morocco", "Western Sahara"),])
# drop wester
world <- world[!world$NOM == "Western Sahara", ]
st_geometry(world[world$NOM == "Morocco", ]) <- st_geometry(moroco)
world <- world[order(world$NOM), ]
europe <- world[world$NOM %in% 
      c("Austria", "Belgium", "Czech Republic", "Denmark", "Finland", "France", 
        "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy", "Luxembourg",
        "Netherlands", "Norway", "Poland", "Portugal", "Spain", "Sweden", 
        "Switzerland", "United Kingdom"), ]

europe <- merge(europe, mig_data_dest, by.x = "NOM", by.y = "country_d")


africa <- world[world$NOM %in% 
    c("Algeria", "Angola", "Benin", "Botswana", "Burundi",  "Burkina Faso",  
      "Cameroon", "Cape Verde", "Central African Republic", "Chad", "Comoros",
      "Congo", "Zaire", "Ivory Coast", "Djibouti",  
      "Egypt", "Equatorial Guinea", "Eritrea", "Ethiopia", "Gabon", "Gambia, The", 
      "Ghana", "Guinea", "Guinea-Bissau", "Kenya", "Lesotho", "Liberia", "Libya", 
      "Madagascar", "Malawi", "Mali", "Mauritania", "Mauritius", "Morocco", 
      "Mozambique", "Namibia", "Niger", "Nigeria", "Rwanda", "Senegal", "Seychelles",  
      "Sierra Leone", "Somalia",  "South Africa", "Sudan", # 44
      "Tanzania, United Republic of", "Togo", "Tunisia", "Uganda",  
      "Zambia", "Zimbabwe"),]

# keep the same names
africa$NOM <-  c("Algeria", "Angola", "Benin", "Botswana", 
                 "Burkina Faso", "Burindi",  "Cameroon", "Cape Verde",
                 "Central African Republic", "Chad", "Comoros", "Congo",
                 "Djibouti", "Egypt", "Equatorial Guinea", "Eritrea", 
                 "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Guinea-Bissau",                     
                 "Côte d'Ivoire", "Kenya", "Lesotho", "Liberia", 
                 "Libya", "Madagascar", "Malawi", "Mali", "Mauritania", 
                 "Mauritius", "Morocco", "Mozambique", "Namibia", 
                 "Niger", "Nigeria", "Rwanda", "Senegal", "Seychelles", "Sierra Leone",
                 "Somalia", "South Africa", "Sudan",  "Tanzania", "Togo", "Tunisia",
                 "Uganda", "Congo, Democratic Republic",  "Zambia", "Zimbabwe")  

africa <- merge(africa, mig_data_origin, by.x = "NOM", by.y = "country_o")

# re-order data 
mig_data_dest <- mig_data_dest  %>%
  arrange(country_d)
mig_data_origin <- mig_data_origin %>%
  arrange(country_o)

# create the variable Y
mig_data <- mig_data %>%
  arrange(country_o, country_d)
Y_mig <- matrix(mig_data$migi, n_d_mi, n_o_mi, byrow = F)
