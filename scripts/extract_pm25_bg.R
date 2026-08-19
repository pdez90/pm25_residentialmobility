# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# BLOCK-GROUP PM2.5 PANEL  (WashU / ACAG V6.NA)
# ------------------------------------------------------------
# Builds a block-group x year panel of annual mean PM2.5 from the
# North America V6.NA product (variable "PM25", one NetCDF per year).
# Rasters are extracted to block groups once per year, per-state
# chunked to bound memory, and cached (restartable).
#
# Output: $MIGRATE_DATA_DIR/pollution_bg/bg_pm25_panel.csv
# ============================================================
suppressPackageStartupMessages({
  library(sf); library(data.table); library(terra)
  library(exactextractr); library(stringr); library(tigris)
})
options(tigris_use_cache = TRUE); sf::sf_use_s2(FALSE)

# ---------------------------
# Paths + config
# ---------------------------
base_dir   <- migrate_data_dir()
pm_dir     <- file.path(pm25_raw_dir(), "PM25")     # WashU V6.NA total PM2.5 (mounted)
out_dir    <- file.path(base_dir, "pollution_bg")
cache_dir  <- file.path(out_dir, "cache"); dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
scratch    <- file.path(out_dir, "_terra_tmp"); dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
terraOptions(tempdir = scratch, memfrac = 0.6, progress = 0)

YEARS        <- 2006:2019       # full MIGRATE range (PM2.5 product covers 1998-2024)
BG_VINTAGE   <- 2019            # tigris block groups on 2010 boundaries
CHUNK_STATES <- TRUE            # extract per state to bound memory

# ---------------------------
# National block groups (2010 vintage) -> sf in EPSG:4326
# ---------------------------
message("Loading national block groups (", BG_VINTAGE, ") ...")
sts <- unique(as.data.table(tigris::fips_codes)[as.integer(state_code) <= 56]$state_code)
bg_cache <- file.path(cache_dir, "bg_national.rds")
if (file.exists(bg_cache)) {
  bg <- readRDS(bg_cache)
} else {
  bg <- do.call(rbind, lapply(sts, function(st)
    tryCatch(tigris::block_groups(state = st, year = BG_VINTAGE, cb = TRUE), error = function(e) NULL)))
  bg <- st_transform(st_make_valid(bg), 4326)
  bg$geoid_bg <- str_pad(bg$GEOID, 12, pad = "0")
  bg <- bg[, c("geoid_bg", "STATEFP")]
  saveRDS(bg, bg_cache)
}
message("Block groups: ", nrow(bg))

# ---------------------------
# Extract one raster to all BGs (per-state optional)
# ---------------------------
extract_to_bg <- function(r, value_name) {
  if (is.na(terra::crs(r)) || terra::crs(r) == "") terra::crs(r) <- "EPSG:4326"
  if (!CHUNK_STATES) {
    v <- exact_extract(r, bg, "mean"); v[is.nan(v)] <- NA_real_
    return(data.table(geoid_bg = bg$geoid_bg, value = as.numeric(v)))
  }
  out <- list()
  for (st in sort(unique(bg$STATEFP))) {
    sub <- bg[bg$STATEFP == st, ]
    rc  <- tryCatch(terra::crop(r, terra::ext(vect(sub)) + 0.5), error = function(e) r)
    v   <- exact_extract(rc, sub, "mean"); v[is.nan(v)] <- NA_real_
    out[[st]] <- data.table(geoid_bg = sub$geoid_bg, value = as.numeric(v))
  }
  setnames(rbindlist(out), "value", value_name)[]
}

# ---------------------------
# PM2.5 per year (cached)
# ---------------------------
extract_pm <- function(yr) {
  cf <- file.path(cache_dir, paste0("pm25_", yr, ".csv")); if (file.exists(cf)) return(fread(cf))
  f <- list.files(pm_dir, pattern = sprintf("PM25\\.%d01-%d12\\.nc$", yr, yr), full.names = TRUE)
  if (!length(f)) { message("PM2.5 ", yr, " missing"); return(NULL) }
  message("PM2.5 ", yr, ": ", basename(f[1]))
  r <- terra::rast(f[1], subds = "PM25")
  d <- extract_to_bg(r, "pm25"); d[, year := yr]
  fwrite(d, cf); d
}

panel <- rbindlist(Filter(Negate(is.null), lapply(YEARS, extract_pm)), fill = TRUE)
setcolorder(panel, c("geoid_bg", "year", "pm25"))
fwrite(panel, file.path(out_dir, "bg_pm25_panel.csv"))
message("\nDONE. Panel: ", file.path(out_dir, "bg_pm25_panel.csv"),
        "  (", nrow(panel), " BG-years; pm25 coverage ",
        round(100 * mean(!is.na(panel$pm25)), 1), "%)")
