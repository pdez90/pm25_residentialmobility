# ============================================================
# Pack the MIGRATE origin-destination flows as Parquet for the Shiny
# "Trace a neighbourhood" tab (click any block group -> where its movers went).
#
# The app cannot bundle 23 GB of CSV. Instead each year becomes one Parquet
# file SORTED BY ORIGIN with 500k-row row groups, so DuckDB can answer
#   SELECT ... WHERE Origin = '...'
# over HTTP by reading only the row groups whose min/max span that origin
# (a few MB per query). Host the files anywhere that serves HTTP range
# requests - GitHub Releases is free and works (assets < 2 GB each).
#
# Also writes bg_attr.parquet (one row per block group: location, metro,
# 2011-2019 PM2.5 and racial/ethnic majority) which is small enough to ship
# inside the app (shiny_app/data/).
#
# All modelled flows are kept (MIN_FLOW = 0): at >=1 expected mover only ~11% of
# moves survive, because MIGRATE spreads each origin over thousands of fractional
# flows. Expect ~1 GB per year; GitHub Releases accepts files up to 2 GB.
# Runtime ~30-40 min; needs the 'arrow' package. Run after run_main.R.
#
# Then:
#   gh release create od-v1 app_od/od_*.parquet --repo pdez90/pm25_residentialmobility -t "OD flows (Parquet) for the app"
#   and set OD_PARQUET_BASE=https://github.com/pdez90/pm25_residentialmobility/releases/download/od-v1
#   in the Connect Cloud content's environment variables.
# ============================================================
if (!exists(".migrate_paths_loaded")) source(if (file.exists("scripts/_paths.R")) "scripts/_paths.R" else "_paths.R")
suppressPackageStartupMessages({ library(data.table); library(stringr); library(arrow); library(sf) })
base_dir <- migrate_data_dir(); mig_dir <- base_dir
out_od <- file.path(repo_root, "app_od"); dir.create(out_od, showWarnings=FALSE)
app_dir <- file.path(repo_root, "shiny_app", "data")
MIN_FLOW <- 0      # keep every modelled flow; MIGRATE mass is mostly in fractional flows (>=1 mover keeps only ~11% of moves)
pad_bg <- function(x) str_pad(as.character(x), 12, pad="0")

# ---- block-group attributes (wide by year) ----
poll <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character"), select=c("geoid_bg","year","pm25","no2"))
poll[, geoid_bg := pad_bg(geoid_bg)]; poll <- unique(poll, by=c("geoid_bg","year"))[year %in% 2011:2019]
rcs <- list.files(file.path(base_dir,"pollution_bg"), "^race_bg_.*csv$", full.names=TRUE)
rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}")); r<-fread(f,colClasses=c(geoid_bg="character")); r[,year:=y]; r }), fill=TRUE)
rc[, geoid_bg := pad_bg(geoid_bg)]
rc[, maj := fifelse(sh_black>0.5,"B", fifelse(sh_hisp>0.5,"H", fifelse(sh_white>0.5,"W","N")))]
att <- merge(poll, rc[year %in% 2011:2019, .(geoid_bg, year, maj)], by=c("geoid_bg","year"), all.x=TRUE)
# income and home value (for the 'reachable set': home value <= 3 x origin income, as in the paper)
acs <- fread(file.path(base_dir,"acs_gentrification","blockgroup_gentrification_2009_2019_lag3_with_HMDA.csv"), colClasses=c(geoid_bg="character"),
             select=c("year","geoid_bg","med_hh_income","med_home_value"))
acs[, geoid_bg := pad_bg(geoid_bg)]; acs <- unique(acs[year %in% 2011:2019], by=c("geoid_bg","year"))
att <- merge(att, acs[, .(geoid_bg, year, inc=round(med_hh_income), hv=round(med_home_value))], by=c("geoid_bg","year"), all.x=TRUE)
wide <- dcast(att, geoid_bg ~ year, value.var=c("pm25","maj","inc","hv"))
cxw <- fread(file.path(base_dir,"pollution_bg","county_cbsa_xwalk.csv"), colClasses=c(county="character",cbsa="character"))
bg_cache <- file.path(base_dir,"pollution_bg","cache_multi","bg_national.rds"); if (!file.exists(bg_cache)) bg_cache <- file.path(base_dir,"pollution_bg","cache","bg_national.rds")
bg <- st_transform(readRDS(bg_cache), 4326); cc <- st_coordinates(st_centroid(bg))
cent <- data.table(geoid_bg=pad_bg(bg$geoid_bg), lon=round(cc[,1],4), lat=round(cc[,2],4)); rm(bg)
cent[, county := substr(geoid_bg,1,5)]; cent <- merge(cent, unique(cxw,by="county")[, .(county, cbsa)], by="county", all.x=TRUE)
attr_ <- merge(cent[, .(geoid_bg, lon, lat, cbsa)], wide, by="geoid_bg", all.x=TRUE)
for (y in 2011:2019) set(attr_, j=paste0("pm25_",y), value=round(attr_[[paste0("pm25_",y)]], 2))
write_parquet(attr_, file.path(app_dir, "bg_attr.parquet"), compression="zstd")
cat(sprintf("bg_attr.parquet: %d block groups, %.1f MB\n", nrow(attr_), file.size(file.path(app_dir,"bg_attr.parquet"))/1e6))

# ---- OD flows, one file per year, sorted by origin ----
for (yr in 2011:2019) {
  cat("MIGRATE", yr, "... ")
  x <- fread(file.path(mig_dir, paste0(yr,".csv")), select=c("Origin","Destination","Flow"))
  tot <- sum(x$Flow[x$Origin != x$Destination])       # moves only: MIGRATE also carries stayers as Origin == Destination
  x <- x[Flow > MIN_FLOW & Origin != Destination]
  x[, `:=`(Origin=pad_bg(Origin), Destination=pad_bg(Destination), Flow=round(as.numeric(Flow), 2))]
  setorder(x, Origin)
  f <- file.path(out_od, sprintf("od_%d.parquet", yr))
  write_parquet(x, f, compression="zstd", chunk_size=500000L)
  cat(sprintf("%s rows kept (%.1f%% of moves), %.0f MB\n", format(nrow(x), big.mark=","), 100*sum(x$Flow)/tot, file.size(f)/1e6))
  rm(x); gc()
}
cat("\nDONE ->", out_od, "\nUpload the od_*.parquet files (e.g. GitHub Releases) and set OD_PARQUET_BASE to their folder URL.\n")
