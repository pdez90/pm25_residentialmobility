# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# Download NHGIS ACS 5-year block-group tables via the IPUMS API
#
# Requested tables:
#   B01003 population        B19013 median hh income
#   B25077 median home value B25064 median gross rent
#   B15003 education         B17021 poverty
#   B03002 race/ethnicity    B25003 tenure
#   B25002 vacancy           B25034 year built
#   B19001 income distrib.   B25035 median year built
#   B25024 units in structure
#
# Two known availability wrinkles, handled automatically below:
#   * B15003 does not exist before endyear 2012 -> B15002
#     (Sex by Educational Attainment) is requested as the
#     pre-2012 education table.
#   * B17021 is not published at block-group level in all years
#     -> C17002 (Ratio of Income to Poverty) is requested as a
#     BG-safe poverty fallback in every year.
# The script checks each year's NHGIS metadata and only requests
# table/dataset combos actually offered at block-group level,
# printing whatever it had to skip.
#
# One-time setup:
#   1. Register (free): https://uma.pop.umn.edu/nhgis/user/new
#   2. Get an API key:  https://account.ipums.org/api_keys
#   3. Run: set_ipums_api_key("YOUR_KEY", save = TRUE)
#
# Note: ipumsr >= 0.9 renamed define_extract_nhgis() to
# define_extract_agg("nhgis", ...); both are handled below.
# ============================================================
library(ipumsr)
library(tidyverse)

base_dir <- migrate_data_dir()
nhgis_dir <- file.path(base_dir, "nhgis")
dir.create(nhgis_dir, showWarnings = FALSE, recursive = TRUE)

# Endyears to pull from NHGIS. 2009:2012 fills the pre-API gap;
# set to 2009:2019 if you want ALL years from a single source
# instead of mixing NHGIS (<=2012) with tidycensus (2013+).
endyears <- 2009:2012

tables_wanted <- c(
  "B01003", "B19013", "B25077", "B25064",
  "B15003", # education, endyear >= 2012 only
  "B15002", # education, pre-2012 substitute
  "B17021", # poverty (may be tract-only in some years)
  "C17002", # poverty ratio, BG-safe fallback
  "B03002", "B25003", "B25002", "B25034",
  "B19001", "B25035", "B25024"
)

# ----------------------------
# Discover which NHGIS dataset(s) hold each table at BG level
# ----------------------------
# Each ACS 5-year release is split into sub-datasets (e.g.
# 2008_2012_ACS5a, 2008_2012_ACS5b); typically only "a" goes down
# to block group. We query metadata rather than assume.
all_ds <- get_metadata_nhgis(type = "datasets")

build_specs_for_year <- function(endyr) {
  prefix <- paste0(endyr - 4, "_", endyr, "_ACS5")
  ds_names <- all_ds |>
    filter(str_starts(name, prefix)) |>
    pull(name)
  if (length(ds_names) == 0) {
    warning("No NHGIS ACS5 dataset found for endyear ", endyr)
    return(list())
  }
  specs <- list()
  found <- character()
  for (ds in ds_names) {
    meta <- get_metadata_nhgis(dataset = ds)
    if (!"blck_grp" %in% meta$geog_levels$name) next
    avail <- setdiff(intersect(tables_wanted, meta$data_tables$name), found)
    if (length(avail) == 0) next
    found <- c(found, avail)
    specs <- c(specs, list(
      ds_spec(ds, data_tables = avail, geog_levels = "blck_grp")
    ))
  }
  missing <- setdiff(tables_wanted, found)
  if (length(missing) > 0) {
    message(endyr, ": NOT available at block group, skipped: ",
            paste(missing, collapse = ", "))
  }
  message(endyr, ": requesting ", length(found), " tables from ",
          length(specs), " dataset(s)")
  specs
}

specs <- endyears |>
  map(build_specs_for_year) |>
  list_flatten()

# ----------------------------
# Define, submit, download
# ----------------------------
define_fn <- if (exists("define_extract_agg")) {
  function(...) define_extract_agg("nhgis", ...)
} else {
  define_extract_nhgis
}

ext <- define_fn(
  description = paste0(
    "Block-group gentrification tables, ACS5 endyears ",
    min(endyears), "-", max(endyears)
  ),
  datasets = specs,
  data_format = "csv_header" # csv with descriptive second header row
)

sub <- submit_extract(ext)
sub <- wait_for_extract(sub) # polls until the extract is ready
files <- download_extract(sub, download_dir = nhgis_dir)
message("Downloaded to:\n", paste(files, collapse = "\n"))

# ----------------------------
# Quick inspection of what arrived
# ----------------------------
zips <- list.files(nhgis_dir, pattern = "\\.zip$", full.names = TRUE)
walk(zips, function(z) {
  message("\n== ", basename(z), " ==")
  print(unzip(z, list = TRUE)$Name)
})
# Each csv covers one dataset-year; read with ipumsr::read_nhgis(),
# which attaches variable labels (used by the label-matching
# harmonizer in blockgroup_gentrification_earliest_acs.R).
