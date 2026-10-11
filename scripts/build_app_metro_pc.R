# Principal-city markers for the Shiny metro map (1-2 min). Needs r1_principal_cities.csv
# (from reviewer1_urban_form_distance.R) and shiny_app/data/metro_list.csv.
if (!exists(".migrate_paths_loaded")) source(if (file.exists("scripts/_paths.R")) "scripts/_paths.R" else "_paths.R")
suppressPackageStartupMessages({ library(data.table); library(sf); library(tigris) }); options(tigris_use_cache=TRUE)
out_dir <- file.path(migrate_data_dir(),"pm25_avoidance_lag3_olag1"); app_dir <- file.path(repo_root,"shiny_app","data")
top <- fread(file.path(app_dir,"metro_list.csv"), colClasses=c(cbsa="character"))
pcs <- fread(file.path(out_dir,"r1_principal_cities.csv"), colClasses=c(cbsa="character",place_geoid="character"))[cbsa %in% top$cbsa]
sts <- unique(substr(pcs$place_geoid,1,2))
pl <- do.call(rbind, lapply(sts, function(st) places(state=st, cb=TRUE, year=2019)[,c("GEOID","NAME")]))
pl <- st_transform(pl, 4326); xy <- st_coordinates(st_centroid(pl))
pld <- data.table(place_geoid=pl$GEOID, lon=round(xy[,1],4), lat=round(xy[,2],4))
res <- merge(pcs, pld, by="place_geoid")[, .(cbsa, principal_city, lon, lat)]
fwrite(res, file.path(app_dir,"metro_principal_city.csv")); cat("wrote", nrow(res), "principal cities to", file.path(app_dir,"metro_principal_city.csv"), "\n")
