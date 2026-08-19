# ============================================================
# BLOCK-GROUP MULTI-POLLUTANT PANEL
#   pm25           WashU/ACAG V6.NA total PM2.5 (var "PM25")
#   no2            figshare Global Surface NO2 (annual GeoTIFF)
#   wildfire_pm25  Childs et al. 2022 daily smokePM on a 10 km grid,
#                  aggregated to an ANNUAL mean per grid cell (including
#                  zero-smoke days) and mapped to block groups
#   ozone          Yale/Liu et al. 2026 warm-season (May-Oct) mean MDA8
#                  ozone, annual NetCDF (var "WSMDA8")
#
# Hard-fails if any configured input path/file is missing. Per-year
# gaps within a legitimately shorter series (e.g. ozone starts 2010)
# are warned, not fatal.
#
# Output: pollution_bg/bg_pollution_panel.csv
# ============================================================
suppressPackageStartupMessages({
  library(sf); library(data.table); library(terra); library(exactextractr)
  library(stringr); library(tigris); library(lubridate)
})
options(tigris_use_cache = TRUE); sf::sf_use_s2(FALSE)

base_dir  <- "/Users/priyanka/MIGRATE_dir"
nav_dir   <- "/Users/priyanka/Annual_PM_NAV6"                                  # WashU V6.NA (PM25/ subfolder)
no2_dir   <- "/Users/priyanka/Downloads/DHS_India/NO2"                         # figshare Global Surface NO2 annual GeoTIFFs
ozone_dir <- "/Users/priyanka/Ozone/annual"                                    # Yale annual MDA8 NetCDFs (var WSMDA8)
ozone_var <- "WSMDA8"
# --- wildfire smoke PM2.5: Childs et al. (grid shapefile + daily predictions RDS) ---
grid_shp  <- "/Users/priyanka/Downloads/10km_grid/10km_grid_wgs84/10km_grid_wgs84.shp"
pred_rds  <- "/Users/priyanka/Downloads/10km_grid/smokePM2pt5_predictions_daily_10km_20060101-20231231.rds"

out_dir   <- file.path(base_dir, "pollution_bg"); cache_dir <- file.path(out_dir, "cache_multi")
dir.create(cache_dir, showWarnings=FALSE, recursive=TRUE)
scratch   <- file.path(out_dir, "_terra_tmp_multi"); dir.create(scratch, showWarnings=FALSE, recursive=TRUE)
terraOptions(tempdir=scratch, memfrac=0.6, progress=0)
YEARS <- 2006:2019; BG_VINTAGE <- 2019; CHUNK_STATES <- TRUE
NO2_PATTERN <- "\\.tif$"

first_existing <- function(nms, cands){ h<-intersect(cands,nms); if(length(h)) h[1] else NA_character_ }

# ---- required inputs: fail early if missing ----
stopifnot("WashU PM25 folder not found" = dir.exists(file.path(nav_dir,"PM25")))
stopifnot("NO2 directory not found: fix no2_dir"      = dir.exists(no2_dir))
stopifnot("Ozone directory not found: fix ozone_dir"  = dir.exists(ozone_dir))
stopifnot("Wildfire grid shapefile not found: fix grid_shp" = file.exists(grid_shp))
stopifnot("Wildfire predictions RDS not found: fix pred_rds" = file.exists(pred_rds))

# ---- national block groups ----
sts <- unique(as.data.table(tigris::fips_codes)[as.integer(state_code)<=56]$state_code)
bg_cache <- file.path(cache_dir, "bg_national.rds")
if (file.exists(bg_cache)) bg <- readRDS(bg_cache) else {
  bg <- do.call(rbind, lapply(sts, function(st) tryCatch(tigris::block_groups(state=st, year=BG_VINTAGE, cb=TRUE), error=function(e) NULL)))
  bg <- st_transform(st_make_valid(bg), 4326); bg$geoid_bg <- str_pad(bg$GEOID,12,pad="0"); bg <- bg[,c("geoid_bg","STATEFP")]; saveRDS(bg, bg_cache)
}
message("Block groups: ", nrow(bg))

extract_to_bg <- function(r, nm) {
  if (is.na(terra::crs(r)) || terra::crs(r)=="") terra::crs(r) <- "EPSG:4326"
  if (!CHUNK_STATES) { v<-exact_extract(r,bg,"mean"); v[is.nan(v)]<-NA_real_; return(data.table(geoid_bg=bg$geoid_bg, value=as.numeric(v))) }
  out <- list()
  for (st in sort(unique(bg$STATEFP))) { sub<-bg[bg$STATEFP==st,]; rc<-tryCatch(terra::crop(r, terra::ext(vect(sub))+0.5), error=function(e) r)
    v<-exact_extract(rc,sub,"mean"); v[is.nan(v)]<-NA_real_; out[[st]]<-data.table(geoid_bg=sub$geoid_bg, value=as.numeric(v)) }
  setnames(rbindlist(out),"value",nm)[]
}
washu_year <- function(comp, nm, yr) {
  cf <- file.path(cache_dir, sprintf("%s_%d.csv", nm, yr)); if (file.exists(cf)) return(fread(cf, colClasses=c(geoid_bg="character")))
  f <- list.files(file.path(nav_dir, comp), pattern=sprintf("\\.%s\\.%d01-%d12\\.nc$", comp, yr, yr), full.names=TRUE)
  if (!length(f)) { message(comp," ",yr," missing"); return(NULL) }
  if (length(f) > 1) stop("Expected exactly 1 ", comp, " file for ", yr, "; found ", length(f), ":\n  ", paste(f, collapse="\n  "))
  d <- extract_to_bg(terra::rast(f[1], subds=comp), nm); d[, year:=yr]; fwrite(d, cf); d
}
annual_raster_year <- function(dir, pat, yr, nm, subds=NULL) {
  cf <- file.path(cache_dir, sprintf("%s_%d.csv", nm, yr)); if (file.exists(cf)) return(fread(cf, colClasses=c(geoid_bg="character")))
  fs <- list.files(dir, pattern=pat, full.names=TRUE, recursive=TRUE)
  fs <- fs[str_detect(basename(fs), as.character(yr))]; if (!length(fs)) { message(nm," ",yr," missing"); return(NULL) }
  if (length(fs) > 1) stop("Expected exactly 1 ", nm, " file for ", yr, "; found ", length(fs), ":\n  ",
                           paste(fs, collapse="\n  "), "\n  Tighten the pattern or move alternate/uncertainty/tile files out of the directory.")
  r <- if (!is.null(subds)) terra::rast(fs[1], subds=subds) else { rr<-terra::rast(fs[1]); if (terra::nlyr(rr)>1) rr[[1]] else rr }
  d <- extract_to_bg(r, nm); d[, year:=yr]; fwrite(d, cf); d
}

# ---- WILDFIRE (Childs): daily grid predictions -> annual per grid -> BG (centroid-in-grid) ----
extract_wildfire <- function() {
  cf <- file.path(cache_dir, "wildfire_pm25_ALL.csv"); if (file.exists(cf)) return(fread(cf, colClasses=c(geoid_bg="character")))
  message("Reading Childs smokePM predictions RDS (large) ...")
  pred <- as.data.table(readRDS(pred_rds))
  gcol <- first_existing(names(pred), c("grid_id_10km","grid_id","GRID_ID","ID","id"))
  dcol <- first_existing(names(pred), c("date","Date","DATE"))
  scol <- first_existing(names(pred), c("smokePM_pred","smokePM","smoke_pm","smokePM2pt5","pred","value"))
  stopifnot("could not detect grid/date/smokePM columns in RDS" = !any(is.na(c(gcol,dcol,scol))))
  pred[, yr := year(as.Date(get(dcol)))]; pred <- pred[yr %in% YEARS]
  # The Childs product predicts smoke PM2.5 on SMOKE days; non-smoke days are ~0 and are
  # typically absent from the RDS. Dividing the annual sum by CALENDAR days therefore treats
  # absent days as zero-smoke (the intended annual-mean definition). Verify no grid-year has
  # MORE day-records than calendar days (which would indicate duplicates), and report coverage.
  ndchk <- pred[, .(nd = uniqueN(as.Date(get(dcol)))), by=c(gcol,"yr")]
  ndchk[, cal := fifelse(yr %% 4 == 0 & (yr %% 100 != 0 | yr %% 400 == 0), 366L, 365L)]
  if (any(ndchk$nd > ndchk$cal)) stop("wildfire: some grid-years have more day-records than calendar days (duplicate records?).")
  message("wildfire: median smoke-day records/grid-year = ", round(median(ndchk$nd)),
          "; remaining days treated as zero-smoke (annual-mean definition).")
  ann <- pred[, .(ssum=sum(as.numeric(get(scol)), na.rm=TRUE)), by=c(gcol,"yr")]
  ann[, wildfire_pm25 := ssum / fifelse(yr %% 4 == 0 & (yr %% 100 != 0 | yr %% 400 == 0), 366, 365)]
  setnames(ann, gcol, "gid"); ann[, gid := as.character(gid)]
  message("Mapping 10 km grid to block groups (centroid-in-cell) ...")
  grid <- st_transform(st_make_valid(st_read(grid_shp, quiet=TRUE)), 4326)
  ggcol <- first_existing(names(grid), c("grid_id_10km","grid_id","GRID_ID","ID","id"))
  stopifnot("could not detect grid id column in shapefile" = !is.na(ggcol))
  cent <- st_join(st_centroid(bg[,"geoid_bg"]), grid[,ggcol], join=st_within)
  xw <- data.table(geoid_bg=cent$geoid_bg, gid=as.character(cent[[ggcol]]))[!is.na(gid)]
  xw <- unique(xw, by="geoid_bg")                       # one grid cell per BG (boundary centroids can multi-match)
  out <- merge(xw, ann[, .(gid, year=yr, wildfire_pm25)], by="gid", allow.cartesian=TRUE)[, .(geoid_bg, year, wildfire_pm25)]
  out <- unique(out, by=c("geoid_bg","year"))
  fwrite(out, cf); out
}

# ============================================================
# EXTRACT each pollutant; hard-fail if a configured source yields nothing
# ============================================================
pm  <- rbindlist(Filter(Negate(is.null), lapply(YEARS, function(y) washu_year("PM25","pm25",y))), fill=TRUE)
stopifnot("PM2.5: no years extracted (check nav_dir/PM25)" = nrow(pm)>0)
no2 <- rbindlist(Filter(Negate(is.null), lapply(YEARS, function(y) annual_raster_year(no2_dir, NO2_PATTERN, y, "no2"))), fill=TRUE)
stopifnot("NO2: no years extracted (check no2_dir / filenames contain the year)" = nrow(no2)>0)
oz  <- rbindlist(Filter(Negate(is.null), lapply(YEARS, function(y) annual_raster_year(ozone_dir, "\\.nc$", y, "ozone", subds=ozone_var))), fill=TRUE)
stopifnot("Ozone: no years extracted (check ozone_dir / var WSMDA8)" = nrow(oz)>0)
wf  <- extract_wildfire()
stopifnot("Wildfire: no rows produced (check grid_shp / pred_rds)" = nrow(wf)>0)

# CRITICAL: force geoid_bg to a 12-char zero-padded string on every
# component before merging. Cache files read with fread() can parse
# geoid_bg as numeric and drop the leading zero for FIPS 01-09 (AL, AK,
# AZ, AR, CA, CO, CT, ...), which makes those states fail to merge and
# duplicates their BG-year rows. Block-group GEOIDs are always 12 digits,
# so padding restores the correct key.
.padbg <- function(dt) { dt[, geoid_bg := str_pad(as.character(geoid_bg), 12, pad="0")]; dt }
pm <- .padbg(pm); no2 <- .padbg(no2); wf <- .padbg(wf); oz <- .padbg(oz)
panel <- Reduce(function(a,b) merge(a,b,by=c("geoid_bg","year"),all=TRUE), list(pm,no2,wf,oz))
panel <- unique(panel, by=c("geoid_bg","year"))   # guard against any residual duplicate keys
setcolorder(panel, c("geoid_bg","year", intersect(c("pm25","no2","wildfire_pm25","ozone"), names(panel))))
fwrite(panel, file.path(out_dir, "bg_pollution_panel.csv"))
.p <- c("pm25","no2","wildfire_pm25","ozone")
message("\nDONE. Multi-pollutant panel (", nrow(panel), " BG-years). Coverage: ",
        paste(sprintf("%s %.0f%%", .p, sapply(.p, function(c) if (c %in% names(panel)) 100*mean(!is.na(panel[[c]])) else 0)), collapse=" | "))
