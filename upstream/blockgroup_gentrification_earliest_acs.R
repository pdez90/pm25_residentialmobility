# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# ACS block-group gentrification measures, earliest ACS -> 2019
#
# DATA SOURCES:
#   * tidycensus / Census API: endyears 2013-2019
#   * NHGIS extract zip(s): endyears 2009-2012 in `nhgis_dir`
#   * Endyear 2020 excluded because it uses 2020 boundaries
#
# Measures:
#   A. Multidimensional transition score
#   B. Bunten expectations score: home value percentile - income percentile
#   C. Classic Freeman/Bates-style binary measure
#
# CHANGES vs previous version:
#   * TOPCODE HANDLING: ACS topcodes median home value at
#     $1,000,001 (earlier releases) / $2,000,001 (later releases).
#     Changes involving a topcoded value are meaningless and the
#     topcode ceiling itself rose over time, so home_value_chg is
#     set to NA whenever either endpoint is topcoded
#     (value_topcoded flag kept in the output).
#   * WINSORIZATION: the change inputs to the transition score
#     are winsorized at p1/p99 within year (as *_w copies) before
#     z-scoring, so single extreme values cannot produce scores
#     of 30-45. Raw *_chg columns are unchanged in the output;
#     Measure C still uses raw changes.
# ============================================================
library(tidycensus)
library(tidyverse)
library(data.table)
library(janitor)
library(scales)

# ----------------------------
# Paths & switches
# ----------------------------
base_dir <- Sys.getenv("MIGRATE_DIR", unset=migrate_data_dir())
out_dir <- file.path(base_dir, "acs_gentrification")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
nhgis_dir <- file.path(base_dir, "nhgis")
walk(c(out_dir, fig_dir, tab_dir, nhgis_dir),
     ~ dir.create(.x, showWarnings = FALSE, recursive = TRUE))

use_nhgis <- TRUE
# SENSITIVITY: set change_lag <- 5 and re-run this script plus the
# HMDA add-on for the fully non-overlapping robustness version
# (5-year ACS windows share no sample years at lag 5+; analysis
# window becomes 2015-2019). All filenames and provenance columns
# adjust automatically, and the add-on's section 10 will compare
# the two lag versions once both exist.
change_lag <- 3
years_api <- 2013:2019
states <- c(state.abb, "DC")

clear_cached_acs <- FALSE
if (clear_cached_acs) {
  file.remove(list.files(out_dir, pattern = "^acs_bg_51_20.*\\.csv$", full.names = TRUE))
}

# ----------------------------
# ACS variables, tidycensus endyears 2013+
# ----------------------------
vars <- c(
  total_pop      = "B01003_001",
  med_hh_income  = "B19013_001",
  med_gross_rent = "B25064_001",
  med_home_value = "B25077_001",
  total_25plus   = "B15003_001",
  bach_plus_1    = "B15003_022",
  bach_plus_2    = "B15003_023",
  bach_plus_3    = "B15003_024",
  bach_plus_4    = "B15003_025",
  poverty_den    = "B17021_001",
  poverty_num    = "B17021_002"
)

# ----------------------------
# Helpers
# ----------------------------
safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}
safe_percent_rank <- function(x) {
  out <- rep(NA_real_, length(x))
  idx <- !is.na(x)
  if (sum(idx) < 2) return(out)
  vals <- x[idx]
  if (length(unique(vals)) < 2) return(out)
  out[idx] <- percent_rank(vals)
  out
}
z_eligible <- function(x, eligible) {
  out <- rep(NA_real_, length(x))
  idx <- !is.na(eligible) & eligible == TRUE & !is.na(x)
  if (sum(idx) < 2) return(out)
  vals <- x[idx]
  sdx <- sd(vals, na.rm = TRUE)
  if (is.na(sdx) || !is.finite(sdx) || sdx == 0) return(out)
  out[idx] <- (vals - mean(vals, na.rm = TRUE)) / sdx
  out
}
pct_rank_eligible <- function(x, eligible) {
  out <- rep(NA_real_, length(x))
  idx <- !is.na(eligible) & eligible == TRUE & !is.na(x)
  if (sum(idx) < 2) return(out)
  vals <- x[idx]
  if (length(unique(vals)) < 2) return(out)
  out[idx] <- percent_rank(vals)
  out
}
jaccard <- function(a, b) {
  both <- sum(a == 1 & b == 1, na.rm = TRUE)
  either <- sum(a == 1 | b == 1, na.rm = TRUE)
  if (either == 0) return(NA_real_)
  both / either
}
# cor() that returns NA instead of erroring when a year has no
# complete pairs (matches the helper in the HMDA add-on)
safe_cor <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok])
}
# Winsorize at p / 1-p (returns x unchanged if too few values).
wins <- function(x, p = 0.01) {
  if (sum(!is.na(x)) < 2) return(x)
  qs <- quantile(x, c(p, 1 - p), na.rm = TRUE, names = FALSE)
  pmin(pmax(x, qs[1]), qs[2])
}

# ----------------------------
# tidycensus download/cache
# ----------------------------
clean_bg_acs <- function(dat, yr) {
  dat |>
    transmute(
      year = yr,
      geoid_bg = GEOID,
      total_pop = total_popE,
      med_hh_income = med_hh_incomeE,
      med_gross_rent = med_gross_rentE,
      med_home_value = med_home_valueE,
      total_25plus = total_25plusE,
      bach_plus = bach_plus_1E + bach_plus_2E + bach_plus_3E + bach_plus_4E,
      poverty_den = poverty_denE,
      poverty_num = poverty_numE
    )
}
get_bg_acs_state <- function(st, yr) {
  message("Downloading ACS block groups: ", st, " ", yr)
  tryCatch(
    get_acs(
      geography = "block group",
      variables = vars,
      state = st,
      year = yr,
      survey = "acs5",
      output = "wide",
      geometry = FALSE
    ) |>
      clean_bg_acs(yr),
    error = function(e) {
      message("  failed: ", st, " ", yr, " | ", e$message)
      tibble()
    }
  )
}
acs_bg <- map_dfr(years_api, function(yr) {
  yr_file <- file.path(out_dir, paste0("acs_bg_51_", yr, ".csv"))
  if (file.exists(yr_file) && file.size(yr_file) > 0) {
    fread(yr_file, colClasses = c(geoid_bg = "character")) |>
      as_tibble()
  } else {
    dat <- map_dfr(states, ~ get_bg_acs_state(.x, yr))
    if (nrow(dat) > 0) fwrite(dat, yr_file)
    dat
  }
})

# ----------------------------
# NHGIS ingestion, endyears 2009-2012
# ----------------------------
read_nhgis_extract <- function(zip_path) {
  contents <- unzip(zip_path, list = TRUE)$Name
  csvs <- str_subset(contents, "blck_grp\\.csv$")
  map_dfr(csvs, function(csv_name) {
    message("Reading NHGIS: ", basename(csv_name))
    cb_name <- str_replace(csv_name, "\\.csv$", "_codebook.txt")
    tmp <- tempdir()
    unzip(zip_path, files = c(csv_name, cb_name), exdir = tmp, overwrite = TRUE)
    cb <- readLines(file.path(tmp, cb_name), warn = FALSE)
    src <- str_match(cb, "Source code:\\s*(\\S+)")[, 2]
    nhg <- str_match(cb, "NHGIS code:\\s*(\\S+)")[, 2]
    src_idx <- which(!is.na(src))
    lookup <- tibble(
      source = src[src_idx],
      prefix = nhg[src_idx + 1]
    )
    stopifnot("codebook parse failed" = !anyNA(lookup$prefix))
    dat <- fread(file.path(tmp, csv_name), colClasses = "character")
    if ("GISJOIN" %in% names(dat) && dat$GISJOIN[1] == "GIS Join Match Code") {
      dat <- dat[-1]
    }
    col <- function(tbl, i) {
      pfx <- lookup$prefix[match(tbl, lookup$source)]
      if (is.na(pfx)) return(rep(NA_real_, nrow(dat)))
      suppressWarnings(as.numeric(dat[[paste0(pfx, "E", sprintf("%03d", i))]]))
    }
    col_sum <- function(tbl, is) {
      # do.call(cbind, ...) always yields a matrix, even with a
      # single column - sapply() would collapse that case to a
      # bare vector (its simplify= argument doesn't prevent this)
      # and rowSums() would then fail.
      mat <- do.call(cbind, lapply(is, function(i) col(tbl, i)))
      out <- rowSums(mat, na.rm = TRUE)
      ifelse(rowSums(!is.na(mat)) == 0, NA_real_, out)
    }
    endyear <- as.integer(str_match(csv_name, "_(20\\d{2})5_blck_grp")[, 2])
    if ("B15003" %in% lookup$source) {
      edu_tot <- col("B15003", 1)
      edu_bach <- col_sum("B15003", 22:25)
    } else {
      edu_tot <- col("B15002", 1)
      edu_bach <- col_sum("B15002", c(15:18, 32:35))
    }
    tibble(
      year = endyear,
      geoid_bg = paste0(
        str_pad(dat$STATEA, 2, pad = "0"),
        str_pad(dat$COUNTYA, 3, pad = "0"),
        str_pad(dat$TRACTA, 6, pad = "0"),
        str_pad(dat$BLKGRPA, 1, pad = "0")
      ),
      total_pop = col("B01003", 1),
      med_hh_income = col("B19013", 1),
      med_gross_rent = col("B25064", 1),
      med_home_value = col("B25077", 1),
      total_25plus = edu_tot,
      bach_plus = edu_bach,
      poverty_den = col("B17021", 1),
      poverty_num = col("B17021", 2)
    ) |>
      filter(str_sub(geoid_bg, 1, 2) != "72")
  })
}
if (use_nhgis) {
  nhgis_zips <- list.files(nhgis_dir, pattern = "\\.zip$", full.names = TRUE)
  if (length(nhgis_zips) > 0) {
    nhgis_bg <- map_dfr(nhgis_zips, read_nhgis_extract) |>
      filter(year <= 2012)
    acs_bg <- bind_rows(nhgis_bg, acs_bg)
    message(
      "Appended NHGIS endyears: ",
      paste(sort(unique(nhgis_bg$year)), collapse = ", ")
    )
  } else {
    message(
      "No NHGIS zips found in ",
      nhgis_dir,
      " - panel starts at endyear ",
      min(years_api)
    )
  }
}

# ----------------------------
# Clean ACS panel
# ----------------------------
acs_clean <- acs_bg |>
  mutate(
    geoid_bg = str_pad(geoid_bg, 12, pad = "0"),
    # Harmonize county FIPS changes.
    geoid_bg = str_replace(geoid_bg, "^46113", "46102"),
    geoid_bg = str_replace(geoid_bg, "^51515", "51019"),
    geoid_tract = str_sub(geoid_bg, 1, 11),
    statefp = str_sub(geoid_bg, 1, 2),
    countyfp = str_sub(geoid_bg, 3, 5),
    state_county = str_sub(geoid_bg, 1, 5),
    # 2005-2009 ACS is on Census 2000 block-group boundaries.
    bg_vintage = if_else(year <= 2009, "2000", "2010"),
    across(
      c(med_hh_income, med_gross_rent, med_home_value),
      ~ ifelse(.x < 0, NA_real_, .x)
    ),
    bach_pct = 100 * bach_plus / total_25plus,
    poverty_pct = 100 * poverty_num / poverty_den,
    across(
      c(bach_pct, poverty_pct),
      ~ ifelse(is.nan(.x) | is.infinite(.x) | .x < 0 | .x > 100, NA_real_, .x)
    )
  ) |>
  distinct(year, geoid_bg, .keep_all = TRUE)

# ----------------------------
# County-year levels and percentiles (needed before baseline join)
# ----------------------------
acs_lvl <- acs_clean |>
  group_by(state_county, year) |>
  mutate(
    county_med_income = safe_median(med_hh_income),
    county_med_home_value = safe_median(med_home_value),
    income_pctile_county = safe_percent_rank(med_hh_income),
    home_value_pctile_county = safe_percent_rank(med_home_value),
    # Measure B, level version.
    expectation_gap = home_value_pctile_county - income_pctile_county
  ) |>
  ungroup()

# ----------------------------
# Gap-safe baselines: explicit self-join on year - change_lag
# ----------------------------
baseline <- acs_lvl |>
  transmute(
    geoid_bg,
    bg_vintage,
    year = year + change_lag, # attach to the comparison year
    base_med_hh_income = med_hh_income,
    base_med_gross_rent = med_gross_rent,
    base_med_home_value = med_home_value,
    base_bach_pct = bach_pct,
    base_poverty_pct = poverty_pct,
    base_county_med_income = county_med_income,
    base_county_med_home_value = county_med_home_value
  )

first_2010_vintage_year <- min(acs_clean$year[acs_clean$bg_vintage == "2010"], na.rm = TRUE)
analysis_start_year <- first_2010_vintage_year + change_lag

# ----------------------------
# Gentrification measures
# ----------------------------
acs_gent <- acs_lvl |>
  left_join(
    baseline,
    by = c("geoid_bg", "bg_vintage", "year"),
    relationship = "one-to-one" # errors loudly on duplicate keys
  ) |>
  mutate(
    income_chg = med_hh_income - base_med_hh_income,
    rent_chg = med_gross_rent - base_med_gross_rent,
    home_value_chg = med_home_value - base_med_home_value,
    bach_chg = bach_pct - base_bach_pct,
    poverty_chg = poverty_pct - base_poverty_pct,
    # --- TOPCODE HANDLING -----------------------------------
    # ACS topcodes B25077 at $1,000,001 (earlier releases) or
    # $2,000,001 (later releases). A change with a topcoded
    # endpoint is meaningless (and the ceiling itself rose over
    # time), so those changes are set to NA. Using >= 1,000,001
    # captures both historical ACS topcode thresholds
    # ($1,000,001 and $2,000,001).
    value_topcoded =
      (!is.na(med_home_value) & med_home_value >= 1000001) |
      (!is.na(base_med_home_value) & base_med_home_value >= 1000001),
    home_value_chg = ifelse(value_topcoded, NA_real_, home_value_chg),
    # ---------------------------------------------------------
    low_income_base = base_med_hh_income < base_county_med_income,
    low_home_value_base = base_med_home_value < base_county_med_home_value,
    eligible_gentrifiable = case_when(
      low_income_base == TRUE | low_home_value_base == TRUE ~ TRUE,
      low_income_base == FALSE & low_home_value_base == FALSE ~ FALSE,
      TRUE ~ NA
    ),
    # ELIGIBILITY SENSITIVITY (additive; the PRIMARY rule is the OR above).
    # Alternatives used only in robustness analyses, primary classification unchanged:
    eligible_both        = low_income_base == TRUE & low_home_value_base == TRUE,   # stricter (both below county median)
    eligible_income_only = low_income_base == TRUE,                                 # income only
    eligible_value_only  = low_home_value_base == TRUE                              # home value only
  ) |>
  group_by(state_county, year) |>
  mutate(
    county_income_chg = safe_median(income_chg),
    county_rent_chg = safe_median(rent_chg),
    county_home_value_chg = safe_median(home_value_chg),
    county_bach_chg = safe_median(bach_chg)
  ) |>
  ungroup() |>
  mutate(
    rising_rent = rent_chg > county_rent_chg,
    rising_home_value = home_value_chg > county_home_value_chg,
    rising_ses = income_chg > county_income_chg | bach_chg > county_bach_chg,
    # Measure C (uses raw changes).
    gentrifying_classic_binary = as.integer(
      eligible_gentrifiable == TRUE &
        ((rising_rent == TRUE) | (rising_home_value == TRUE)) &
        rising_ses == TRUE
    )
  ) |>
  group_by(year) |>
  mutate(
    # --- WINSORIZED COPIES for Measure A ---------------------
    # p1/p99 within year, so single extreme values cannot
    # dominate the z-scores. Raw *_chg columns stay untouched.
    income_chg_w = wins(income_chg),
    rent_chg_w = wins(rent_chg),
    home_value_chg_w = wins(home_value_chg),
    bach_chg_w = wins(bach_chg),
    poverty_chg_w = wins(poverty_chg),
    # ---------------------------------------------------------
    # Measure A (z-scores on winsorized changes).
    z_rent = z_eligible(rent_chg_w, eligible_gentrifiable),
    z_home_value = z_eligible(home_value_chg_w, eligible_gentrifiable),
    z_income = z_eligible(income_chg_w, eligible_gentrifiable),
    z_bach = z_eligible(bach_chg_w, eligible_gentrifiable),
    z_poverty = z_eligible(poverty_chg_w, eligible_gentrifiable),
    gent_score_transition =
      z_rent + z_home_value + z_income + z_bach - z_poverty,
    gent_score_transition_no_rent =
      z_home_value + z_income + z_bach - z_poverty,
    gent_score_transition_pctile =
      pct_rank_eligible(gent_score_transition, eligible_gentrifiable),
    gent_score_transition_no_rent_pctile =
      pct_rank_eligible(gent_score_transition_no_rent, eligible_gentrifiable),
    gentrifying_transition_top10 = as.integer(
      eligible_gentrifiable == TRUE &
        !is.na(gent_score_transition_pctile) &
        gent_score_transition_pctile >= 0.90
    ),
    gentrifying_transition_top05 = as.integer(
      eligible_gentrifiable == TRUE &
        !is.na(gent_score_transition_pctile) &
        gent_score_transition_pctile >= 0.95
    ),
    gentrifying_transition_top10_no_rent = as.integer(
      eligible_gentrifiable == TRUE &
        !is.na(gent_score_transition_no_rent_pctile) &
        gent_score_transition_no_rent_pctile >= 0.90
    ),
    gentrifying_transition_top05_no_rent = as.integer(
      eligible_gentrifiable == TRUE &
        !is.na(gent_score_transition_no_rent_pctile) &
        gent_score_transition_no_rent_pctile >= 0.95
    ),
    # Measure B.
    expectation_gap_pctile =
      pct_rank_eligible(expectation_gap, eligible_gentrifiable),
    gentrifying_expectation_top10 = as.integer(
      eligible_gentrifiable == TRUE &
        !is.na(expectation_gap_pctile) &
        expectation_gap_pctile >= 0.90
    ),
    gentrifying_expectation_top05 = as.integer(
      eligible_gentrifiable == TRUE &
        !is.na(expectation_gap_pctile) &
        expectation_gap_pctile >= 0.95
    )
  ) |>
  ungroup() |>
  mutate(
    # Backward-compatible names.
    gent_score = gent_score_transition,
    gent_score_no_rent = gent_score_transition_no_rent,
    gentrifying_top10 = gentrifying_transition_top10,
    gentrifying_top05 = gentrifying_transition_top05,
    gentrifying_top10_no_rent = gentrifying_transition_top10_no_rent,
    gentrifying_top05_no_rent = gentrifying_transition_top05_no_rent,
    gentrifying_binary_county = gentrifying_classic_binary,
    # Provenance: change horizon baked into every exported row.
    change_lag_years = change_lag,
    # Explicit analysis window.
    analysis_sample_main = year >= analysis_start_year & year <= 2019
  )

# ----------------------------
# ROBUSTNESS (additive): CBSA-year re-ranking of the transition score.
# The PRIMARY classification above ranks nationally within year; this adds a
# LOCAL-MARKET alternative that re-ranks the SAME gent_score_transition within
# each CBSA x year among eligible neighborhoods (top decile of the metro rather
# than of the nation). All primary columns are unchanged. Skipped without error
# if the county->CBSA crosswalk is absent.
# ----------------------------
.xwf <- file.path(base_dir, "pollution_bg", "county_cbsa_xwalk.csv")
if (file.exists(.xwf)) {
  .cbsa_xw <- fread(.xwf, colClasses = c(county = "character", cbsa = "character"))
  acs_gent <- acs_gent |>
    left_join(.cbsa_xw, by = c("state_county" = "county")) |>
    mutate(cbsa = ifelse(is.na(cbsa), paste0("nonmetro_", state_county), cbsa)) |>
    group_by(cbsa, year) |>
    mutate(
      gent_score_transition_pctile_cbsa = pct_rank_eligible(gent_score_transition, eligible_gentrifiable),
      gentrifying_transition_top10_cbsa = as.integer(eligible_gentrifiable == TRUE &
          !is.na(gent_score_transition_pctile_cbsa) & gent_score_transition_pctile_cbsa >= 0.90),
      gentrifying_transition_top05_cbsa = as.integer(eligible_gentrifiable == TRUE &
          !is.na(gent_score_transition_pctile_cbsa) & gent_score_transition_pctile_cbsa >= 0.95)
    ) |>
    ungroup()
  cat("CBSA-year robustness ranking added:",
      sum(acs_gent$gentrifying_transition_top10_cbsa, na.rm = TRUE), "CBSA-top-decile vs",
      sum(acs_gent$gentrifying_transition_top10, na.rm = TRUE), "national-top-decile gentrifying rows\n")
} else message("county_cbsa_xwalk.csv not found; CBSA-year robustness columns skipped.")

# ----------------------------
# Save combined and separate measure files
# ----------------------------
yr_rng <- paste0(min(acs_gent$year), "_", max(acs_gent$year))
lag_tag <- paste0("lag", change_lag)
out_file <- file.path(
  out_dir,
  paste0("blockgroup_gentrification_", yr_rng, "_", lag_tag, ".csv")
)
fwrite(acs_gent, out_file)
transition_file <- file.path(
  out_dir,
  paste0("transition_gentrification_", yr_rng, "_", lag_tag, ".csv")
)
expectation_file <- file.path(
  out_dir,
  paste0("expectation_gap_gentrification_", yr_rng, "_", lag_tag, ".csv")
)
classic_file <- file.path(
  out_dir,
  paste0("classic_binary_gentrification_", yr_rng, "_", lag_tag, ".csv")
)
acs_gent |>
  select(
    year, geoid_bg, geoid_tract, statefp, countyfp, state_county,
    bg_vintage, change_lag_years, analysis_sample_main, eligible_gentrifiable,
    value_topcoded, # so transition-file users can see when the
                    # home-value component was excluded
    gent_score_transition, gent_score_transition_no_rent,
    gent_score_transition_pctile, gent_score_transition_no_rent_pctile,
    gentrifying_transition_top10, gentrifying_transition_top05,
    gentrifying_transition_top10_no_rent, gentrifying_transition_top05_no_rent
  ) |>
  fwrite(transition_file)
acs_gent |>
  select(
    year, geoid_bg, geoid_tract, statefp, countyfp, state_county,
    bg_vintage, change_lag_years, analysis_sample_main, eligible_gentrifiable,
    expectation_gap, expectation_gap_pctile,
    gentrifying_expectation_top10, gentrifying_expectation_top05
  ) |>
  fwrite(expectation_file)
acs_gent |>
  select(
    year, geoid_bg, geoid_tract, statefp, countyfp, state_county,
    bg_vintage, change_lag_years, analysis_sample_main, eligible_gentrifiable,
    gentrifying_classic_binary
  ) |>
  fwrite(classic_file)
message("Saved combined gentrification file to: ", out_file)
message("Saved transition file to: ", transition_file)
message("Saved expectation-gap file to: ", expectation_file)
message("Saved classic binary file to: ", classic_file)

# ============================================================
# DIAGNOSTICS
# ============================================================
options(pillar.sigfig = 4)

# --- 1. Panel coverage: rows, unique BGs, states per year ----
diag_panel <- acs_clean |>
  group_by(year, bg_vintage) |>
  summarise(
    n_rows = n(),
    n_bg = n_distinct(geoid_bg),
    n_states = n_distinct(statefp),
    .groups = "drop"
  )
message("\n---- 1. Panel coverage by year ----")
print(diag_panel, n = Inf)

# --- 2. Year-gap report (informational; join is gap-safe) ----
gap_report <- acs_clean |>
  group_by(geoid_bg, bg_vintage) |>
  summarise(
    n_years = n(),
    has_gap = n() > 1 && any(diff(sort(year)) != 1),
    .groups = "drop"
  )
message("\n---- 2. BGs with gaps in year sequence ----")
message("BGs with gaps: ", sum(gap_report$has_gap),
        " of ", nrow(gap_report),
        " (gap-safe join means these get NA changes, never wrong-horizon ones)")

# --- 3. Missingness by variable and year ----
diag_missing <- acs_clean |>
  group_by(year) |>
  summarise(
    across(
      c(med_hh_income, med_gross_rent, med_home_value, bach_pct, poverty_pct),
      ~ mean(is.na(.x)),
      .names = "pctNA_{.col}"
    ),
    .groups = "drop"
  )
message("\n---- 3. Share missing by variable and year ----")
print(diag_missing, n = Inf, width = Inf)

# --- 4. Baseline join success, eligibility, topcoding ----
diag_baseline <- acs_gent |>
  group_by(year) |>
  summarise(
    n_bg = n(),
    pct_with_baseline = mean(!is.na(base_med_hh_income)),
    pct_with_income_chg = mean(!is.na(income_chg)),
    pct_value_topcoded = mean(value_topcoded, na.rm = TRUE),
    n_eligible = sum(eligible_gentrifiable == TRUE, na.rm = TRUE),
    pct_eligible = mean(eligible_gentrifiable == TRUE, na.rm = TRUE),
    .groups = "drop"
  )
message("\n---- 4. Baseline (t-", change_lag, ") join success, topcoding by year ----")
print(diag_baseline, n = Inf, width = Inf)

# --- 5. Flag rates by year (all BGs and among eligible) ----
diagnostics <- acs_gent |>
  filter(analysis_sample_main) |>
  group_by(year) |>
  summarise(
    n_bg = n(),
    n_eligible = sum(eligible_gentrifiable == TRUE, na.rm = TRUE),
    pct_eligible = mean(eligible_gentrifiable == TRUE, na.rm = TRUE),
    pct_transition_top10 = mean(gentrifying_transition_top10 == 1, na.rm = TRUE),
    pct_transition_top10_no_rent = mean(gentrifying_transition_top10_no_rent == 1, na.rm = TRUE),
    pct_expectation_top10 = mean(gentrifying_expectation_top10 == 1, na.rm = TRUE),
    pct_classic_binary = mean(gentrifying_classic_binary == 1, na.rm = TRUE),
    # Same flags, but as a share of ELIGIBLE BGs.
    pct_transition_top10_elig =
      mean(gentrifying_transition_top10[eligible_gentrifiable == TRUE] == 1, na.rm = TRUE),
    pct_expectation_top10_elig =
      mean(gentrifying_expectation_top10[eligible_gentrifiable == TRUE] == 1, na.rm = TRUE),
    pct_classic_binary_elig =
      mean(gentrifying_classic_binary[eligible_gentrifiable == TRUE] == 1, na.rm = TRUE),
    .groups = "drop"
  )
message("\n---- 5. Gentrification flag rates by year ----")
print(diagnostics, n = Inf, width = Inf)
fwrite(diagnostics, file.path(tab_dir, paste0("diagnostics_by_year_", lag_tag, ".csv")))

# --- 6. Distribution of the transition score by year ----
diag_score_dist <- acs_gent |>
  filter(analysis_sample_main, !is.na(gent_score_transition)) |>
  group_by(year) |>
  summarise(
    n = n(),
    p05 = quantile(gent_score_transition, 0.05),
    p25 = quantile(gent_score_transition, 0.25),
    p50 = quantile(gent_score_transition, 0.50),
    p75 = quantile(gent_score_transition, 0.75),
    p95 = quantile(gent_score_transition, 0.95),
    p99 = quantile(gent_score_transition, 0.99),
    max = max(gent_score_transition),
    .groups = "drop"
  )
message("\n---- 6. Transition score distribution by year (max should now be single digits) ----")
print(diag_score_dist, n = Inf, width = Inf)
fwrite(diag_score_dist, file.path(tab_dir, paste0("score_distribution_by_year_", lag_tag, ".csv")))

# --- 7. Agreement between measures: correlations + Jaccard ----
measure_correlations <- acs_gent |>
  filter(analysis_sample_main) |>
  group_by(year) |>
  summarise(
    cor_transition_expectation =
      safe_cor(gent_score_transition_no_rent, expectation_gap),
    cor_transition_classic =
      safe_cor(gentrifying_transition_top10_no_rent, gentrifying_classic_binary),
    cor_expectation_classic =
      safe_cor(gentrifying_expectation_top10, gentrifying_classic_binary),
    jac_transition_expectation =
      jaccard(gentrifying_transition_top10, gentrifying_expectation_top10),
    jac_transition_classic =
      jaccard(gentrifying_transition_top10, gentrifying_classic_binary),
    jac_expectation_classic =
      jaccard(gentrifying_expectation_top10, gentrifying_classic_binary),
    .groups = "drop"
  )
message("\n---- 7. Measure agreement by year (cor = correlation, jac = Jaccard overlap of flags) ----")
print(measure_correlations, n = Inf, width = Inf)
fwrite(measure_correlations,
       file.path(tab_dir, paste0("measure_correlations_by_year_", lag_tag, ".csv")))

# --- 8. Flag rates by state (pooled over analysis years) ----
state_lookup <- tidycensus::fips_codes |>
  distinct(statefp = state_code, state_abbr = state)
diag_state <- acs_gent |>
  filter(analysis_sample_main) |>
  group_by(statefp) |>
  summarise(
    n_bg_yrs = n(),
    pct_eligible = mean(eligible_gentrifiable == TRUE, na.rm = TRUE),
    pct_transition_top10_elig =
      mean(gentrifying_transition_top10[eligible_gentrifiable == TRUE] == 1, na.rm = TRUE),
    pct_classic_binary_elig =
      mean(gentrifying_classic_binary[eligible_gentrifiable == TRUE] == 1, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(state_lookup, by = "statefp") |>
  relocate(state_abbr) |>
  arrange(desc(pct_transition_top10_elig))
message("\n---- 8. Flag rates by state, pooled ", analysis_start_year, "-2019 ----")
print(diag_state, n = Inf, width = Inf)
fwrite(diag_state, file.path(tab_dir, paste0("diagnostics_by_state_", lag_tag, ".csv")))

# --- 9. Sanity spot-check: 10 highest-scoring BGs in latest year ----
message("\n---- 9. Top-10 transition-score BGs, 2019 (spot-check against places you know) ----")
acs_gent |>
  filter(year == 2019, gentrifying_transition_top10 == 1) |>
  arrange(desc(gent_score_transition)) |>
  select(geoid_bg, state_county, gent_score_transition, value_topcoded,
         med_hh_income, income_chg, med_home_value, home_value_chg,
         bach_pct, bach_chg) |>
  head(10) |>
  print(width = Inf)

message("\nAll diagnostics saved to: ", tab_dir)
