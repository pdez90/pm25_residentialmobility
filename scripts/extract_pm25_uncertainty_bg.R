# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# BLOCK-GROUP PM2.5 UNCERTAINTY (WashU V6.NA AbsoluteUncertainty)
# ------------------------------------------------------------
# Extracts the per-grid absolute uncertainty (standard-error-like band)
# of the V6.NA PM2.5 product to the block-group-year, for use in the
# measurement-error robustness (inverse-variance weighting).
#
# Files: Annual_uncertainty_PM25/PM25/V6NA01.CNNPM25.AbsoluteUncertainty.YYYYMM-YYYYMM.nc (var "PM25")
# Output: pollution_bg/bg_pm25_uncertainty_panel.csv  (geoid_bg, year, pm25_se)
# ============================================================
suppressPackageStartupMessages({ library(sf); library(data.table); library(terra); library(exactextractr); library(stringr); library(tigris) })
options(tigris_use_cache=TRUE); sf::sf_use_s2(FALSE)
base_dir <- migrate_data_dir()
unc_dir  <- file.path(pm25_unc_dir(), "PM25")
out_dir  <- file.path(base_dir,"pollution_bg"); cache_dir <- file.path(out_dir,"cache_unc")
dir.create(cache_dir, showWarnings=FALSE, recursive=TRUE)
terraOptions(tempdir=file.path(out_dir,"_terra_tmp_unc"), memfrac=0.6, progress=0)
YEARS <- 2006:2019; BG_VINTAGE <- 2019; CHUNK_STATES <- TRUE
stopifnot("uncertainty PM25 folder not found" = dir.exists(unc_dir))

sts <- unique(as.data.table(tigris::fips_codes)[as.integer(state_code)<=56]$state_code)
bg_cache <- file.path(out_dir,"cache_multi","bg_national.rds")
if (!file.exists(bg_cache)) bg_cache <- file.path(out_dir,"cache","bg_national.rds")
if (file.exists(bg_cache)) bg <- readRDS(bg_cache) else {
  bg <- do.call(rbind, lapply(sts, function(st) tryCatch(tigris::block_groups(state=st, year=BG_VINTAGE, cb=TRUE), error=function(e) NULL)))
  bg <- st_transform(st_make_valid(bg),4326); bg$geoid_bg <- str_pad(bg$GEOID,12,pad="0"); bg <- bg[,c("geoid_bg","STATEFP")]
}
extract_to_bg <- function(r,nm){
  if (is.na(terra::crs(r))||terra::crs(r)=="") terra::crs(r) <- "EPSG:4326"
  if (!CHUNK_STATES){ v<-exact_extract(r,bg,"mean"); v[is.nan(v)]<-NA_real_; return(data.table(geoid_bg=bg$geoid_bg, value=as.numeric(v))) }
  out<-list(); for (st in sort(unique(bg$STATEFP))){ sub<-bg[bg$STATEFP==st,]; rc<-tryCatch(terra::crop(r, terra::ext(vect(sub))+0.5), error=function(e) r)
    v<-exact_extract(rc,sub,"mean"); v[is.nan(v)]<-NA_real_; out[[st]]<-data.table(geoid_bg=sub$geoid_bg, value=as.numeric(v)) }
  setnames(rbindlist(out),"value",nm)[]
}
unc_year <- function(yr){
  cf <- file.path(cache_dir, sprintf("pm25se_%d.csv",yr)); if (file.exists(cf)) return(fread(cf))
  f <- list.files(unc_dir, pattern=sprintf("AbsoluteUncertainty\\.%d01-%d12\\.nc$",yr,yr), full.names=TRUE)
  if (!length(f)) { message("uncertainty ",yr," missing"); return(NULL) }
  d <- extract_to_bg(terra::rast(f[1], subds="PM25"), "pm25_se"); d[, year:=yr]; fwrite(d,cf); d
}
panel <- rbindlist(Filter(Negate(is.null), lapply(YEARS, unc_year)), fill=TRUE)
stopifnot("no uncertainty years extracted" = nrow(panel)>0)
fwrite(panel, file.path(out_dir,"bg_pm25_uncertainty_panel.csv"))
message("DONE. PM2.5 uncertainty panel (", nrow(panel), " BG-years); mean SE ",
        round(mean(panel$pm25_se, na.rm=TRUE),2), " ug/m3")
