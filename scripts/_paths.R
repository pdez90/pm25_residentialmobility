# ============================================================
# PATH CONFIGURATION -- the only place file locations are defined.
# ------------------------------------------------------------
# No script in this repository contains an absolute path. Every script
# sources this file and then refers to the helpers below.
#
# Set ONE environment variable to point at your data root:
#
#   export MIGRATE_DATA_DIR=/path/to/your/data
#
# The working root holds the MIGRATE flow files (<year>.csv), the derived
# panels (pollution_bg/), the covariate panel (acs_gentrification/), and the
# analysis outputs (pm25_avoidance_lag3_olag1/). If unset, it defaults to
# ./data inside the repository.
#
# Raw exposure surfaces are large and are usually kept outside the project.
# Each can be pointed somewhere else individually; otherwise they are looked
# for under $MIGRATE_DATA_DIR/raw/. See README.md for the expected layout.
# ============================================================

# --- repository root -------------------------------------------------
# here::here() locates the root from the .here sentinel regardless of the
# working directory. If `here` is not installed we fall back to getwd(),
# which is correct when scripts are run from the repository root.
repo_root <- if (requireNamespace("here", quietly = TRUE)) {
  here::here()
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

# --- working data root ------------------------------------------------
migrate_data_dir <- function() {
  Sys.getenv("MIGRATE_DATA_DIR", unset = file.path(repo_root, "data"))
}

# --- raw exposure surfaces -------------------------------------------
# raw_dir("PM25_RAW_DIR", "acag_pm25") returns $PM25_RAW_DIR if set,
# else $MIGRATE_DATA_DIR/raw/acag_pm25
raw_dir <- function(env_var, default_subdir) {
  Sys.getenv(env_var, unset = file.path(migrate_data_dir(), "raw", default_subdir))
}
raw_file <- function(env_var, ...) {
  Sys.getenv(env_var, unset = file.path(migrate_data_dir(), "raw", ...))
}

# Convenience accessors used by the extraction scripts
pm25_raw_dir     <- function() raw_dir("PM25_RAW_DIR",   "acag_pm25")
pm25_unc_dir     <- function() raw_dir("PM25_UNC_DIR",   "acag_pm25_uncertainty")
no2_raw_dir      <- function() raw_dir("NO2_RAW_DIR",    "no2")
ozone_raw_dir    <- function() raw_dir("OZONE_RAW_DIR",  "ozone")
smoke_rds        <- function() raw_file("SMOKE_PM25_RDS", "smoke",
                       "smokePM2pt5_predictions_daily_10km_20060101-20231231.rds")
smoke_grid_shp   <- function() raw_file("SMOKE_GRID_SHP", "smoke",
                       "10km_grid_wgs84", "10km_grid_wgs84.shp")
greenbook_xls    <- function() raw_file("EPA_GREENBOOK_XLS", "epa_greenbook", "phistory.xls")

.migrate_paths_loaded <- TRUE
