# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# BLOCK-GROUP SPECIATED PM2.5 (WashU / ACAG V6.NA components)
# ------------------------------------------------------------
# Extends extract_pm25_bg.R to the chemical components you already
# have, which proxy emission SOURCES and let the migration analysis
# speak to the sector-disparity literature (Tessum et al. 2021;
# Bekbulat et al. 2024) without new data:
#   SO4  sulfate     -> coal / power generation (+ secondary)
#   NO3  nitrate     -> traffic + agriculture
#   NH4  ammonium    -> agriculture (secondary)
#   OM   organic     -> combustion, biomass, cooking
#   BC   black carbon-> diesel / traffic / combustion
#   DUST crustal     -> road/windblown dust (largely natural)
#   SS   sea salt     -> marine (natural)
# Combustion/industry/traffic (BC, OM, SO4, NO3) are the anthropogenic,
# regulate-able, disparity-relevant fractions; DUST and SS are natural
# baselines.
#
# Output: pollution_bg/bg_pm25_components_panel.csv
#   (geoid_bg, year, pm25, bc, om, so4, no3, nh4, dust, ss)
# ============================================================
suppressPackageStartupMessages({
  library(sf); library(data.table); library(terra); library(exactextractr); library(stringr); library(tigris)
})
options(tigris_use_cache = TRUE); sf::sf_use_s2(FALSE)

base_dir  <- migrate_data_dir()
nav_dir   <- pm25_raw_dir()                     # parent folder with PM25/ BC/ OM/ ... subfolders
out_dir   <- file.path(base_dir, "pollution_bg"); cache_dir <- file.path(out_dir,"cache_comp")
dir.create(cache_dir, showWarnings=FALSE, recursive=TRUE)
scratch   <- file.path(out_dir,"_terra_tmp_comp"); dir.create(scratch, showWarnings=FALSE, recursive=TRUE)
terraOptions(tempdir=scratch, memfrac=0.6, progress=0)

YEARS <- 2006:2019; BG_VINTAGE <- 2019; CHUNK_STATES <- TRUE
COMPONENTS <- c(PM25="pm25", BC="bc", OM="om", SO4="so4", NO3="no3", NH4="nh4", DUST="dust", SS="ss")
# hard-fail if any component folder is missing
for (COMP in names(COMPONENTS))
  stopifnot(setNames(dir.exists(file.path(nav_dir, COMP)), paste0("component folder not found: ", file.path(nav_dir, COMP))))

sts <- unique(as.data.table(tigris::fips_codes)[as.integer(state_code)<=56]$state_code)
bg_cache <- file.path(cache_dir, "bg_national.rds")
if (file.exists(bg_cache)) bg <- readRDS(bg_cache) else {
  bg <- do.call(rbind, lapply(sts, function(st) tryCatch(tigris::block_groups(state=st, year=BG_VINTAGE, cb=TRUE), error=function(e) NULL)))
  bg <- st_transform(st_make_valid(bg), 4326); bg$geoid_bg <- str_pad(bg$GEOID,12,pad="0"); bg <- bg[,c("geoid_bg","STATEFP")]
  saveRDS(bg, bg_cache)
}
message("Block groups: ", nrow(bg))

extract_to_bg <- function(r, nm){
  if (is.na(terra::crs(r)) || terra::crs(r)=="") terra::crs(r) <- "EPSG:4326"
  if (!CHUNK_STATES){ v<-exact_extract(r,bg,"mean"); v[is.nan(v)]<-NA_real_; return(data.table(geoid_bg=bg$geoid_bg, value=as.numeric(v))) }
  out <- list()
  for (st in sort(unique(bg$STATEFP))){ sub<-bg[bg$STATEFP==st,]; rc<-tryCatch(terra::crop(r, terra::ext(vect(sub))+0.5), error=function(e) r)
    v<-exact_extract(rc,sub,"mean"); v[is.nan(v)]<-NA_real_; out[[st]]<-data.table(geoid_bg=sub$geoid_bg, value=as.numeric(v)) }
  setnames(rbindlist(out),"value",nm)[]
}
extract_comp_year <- function(COMP, nm, yr){
  cf <- file.path(cache_dir, sprintf("%s_%d.csv", nm, yr)); if (file.exists(cf)) return(fread(cf))
  f <- list.files(file.path(nav_dir, COMP), pattern=sprintf("\\.%s\\.%d01-%d12\\.nc$", COMP, yr, yr), full.names=TRUE)
  if (!length(f)) { message(COMP," ",yr," missing"); return(NULL) }
  r <- terra::rast(f[1], subds=COMP)
  d <- extract_to_bg(r, nm); d[, year:=yr]; fwrite(d, cf); d
}

panel <- NULL
for (i in seq_along(COMPONENTS)){
  COMP <- names(COMPONENTS)[i]; nm <- COMPONENTS[[i]]
  message("== ", COMP, " ==")
  dd <- rbindlist(Filter(Negate(is.null), lapply(YEARS, function(y) extract_comp_year(COMP, nm, y))), fill=TRUE)
  panel <- if (is.null(panel)) dd else merge(panel, dd, by=c("geoid_bg","year"), all=TRUE)
}
setcolorder(panel, c("geoid_bg","year", unname(COMPONENTS)))
fwrite(panel, file.path(out_dir, "bg_pm25_components_panel.csv"))
message("\nDONE. Speciated panel: ", file.path(out_dir,"bg_pm25_components_panel.csv"),
        " (", nrow(panel), " BG-years, ", length(COMPONENTS), " species)")
