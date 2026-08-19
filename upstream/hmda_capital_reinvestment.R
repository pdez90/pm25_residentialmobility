# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# ADD-ON: HMDA capital reinvestment measure (Measure D)
# Run AFTER blockgroup_gentrification_earliest_acs.R
#
# TIMEFRAME & CONSISTENCY WITH THE OTHER MEASURES:
#   * Input: blockgroup_gentrification_2009_2019_lag<k>.csv
#   * HMDA years 2012-2019. Two reasons 2012 is the floor:
#       - HMDA reports 2000-vintage census tracts before 2012 and
#         2010-vintage tracts from 2012 on; the panel's 2010-2019
#         rows are 2010-vintage, so pre-2012 HMDA would silently
#         join wrong geographies.
#       - The analysis window is 2013-2019 (lag 3), so 2012-2019
#         covers every analysis year.
#   * Sources: 2012-2017 from CFPB historic national LAR files
#     (legacy schema); 2018-2019 from the FFIEC/CFPB data browser
#     API (new schema). Both are harmonized to the same tract-year
#     indicators; a seam diagnostic checks 2017 vs 2018.
#   * All 50 states + DC (PR dropped), matching the panel.
#   * Scoring mirrors Measures A-C: inputs winsorized at p1/p99
#     within year, z-scored among ELIGIBLE block groups only,
#     percentile-ranked within year, top10/top05 flags, and
#     change_lag_years provenance column in every output.
#
# HARMONIZATION NOTES (legacy vs 2018+ schema):
#   * "Investor share" is really NON-PRINCIPAL-RESIDENCE share:
#     legacy owner_occupancy == 2; 2018+ occupancy_type in (2,3).
#     Legacy HMDA cannot separate second homes from investors.
#   * Applications = purchase loans (loan_purpose == 1) with
#     action_taken 1-5 (excludes purchased loans code 6 and
#     preapproval codes 7-8, which are not applications).
#   * 1-4 unit properties only (legacy property_type 1-2;
#     2018+ total_units in 1-4) so multifamily deals don't
#     distort loan-amount medians.
#   * Loan amounts: legacy is $000s -> converted to dollars.
#     Income is $000s in both eras -> dollars.
#
# RELIABILITY THRESHOLD: all tract indicators are suppressed
# (set to NA) for tracts with fewer than 10 purchase
# applications in a year (min_apps below).
#
# DOWNLOAD SIZE WARNING: the 2012-2017 national LAR zips are
# ~1-2 GB EACH (~10 GB total). Each is processed to a small
# tract-year cache once; you can delete the zips afterwards.
# ============================================================
library(tidyverse)
library(data.table)
library(janitor)
library(scales)
library(tidycensus) # for the B25004 seasonal-home screen

# ----------------------------
# Paths & switches
# ----------------------------
base_dir <- Sys.getenv("MIGRATE_DIR", unset=migrate_data_dir())
gent_dir <- file.path(base_dir, "acs_gentrification")
tab_dir <- file.path(gent_dir, "tables")
hmda_dir <- file.path(base_dir, "hmda")
hud_dir <- file.path(base_dir, "hud_usps_vacancy")
walk(c(hmda_dir, hud_dir, tab_dir),
     ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))

change_lag <- 3                  # must match the main pipeline
lag_tag <- paste0("lag", change_lag)
legacy_years <- 2012:2017        # CFPB historic national files
api_years <- 2018:2019           # FFIEC data browser API
hmda_years <- c(legacy_years, api_years)
states <- c(state.abb, "DC")     # matches the panel
min_apps <- 10                   # reliability screen: indicators -> NA
                                 # for tracts with fewer purchase apps
# Seasonal/second-home screen: BGs where seasonal units (B25004_006)
# exceed this share of housing stock are resort markets, not
# gentrifying neighborhoods; they are excluded from the SCREENED
# capital flags. One reference year suffices - seasonal character
# is highly persistent.
seasonal_share_max <- 0.30
seasonal_ref_year <- 2019
options(timeout = 3600)          # legacy files are huge

gent_file <- file.path(
  gent_dir,
  paste0("blockgroup_gentrification_2009_2019_", lag_tag, ".csv")
)
stopifnot("run the main pipeline first" = file.exists(gent_file))
acs_gent <- fread(
  gent_file,
  colClasses = c(
    geoid_bg = "character", geoid_tract = "character",
    statefp = "character", countyfp = "character",
    state_county = "character"
  )
) |>
  as_tibble() |>
  mutate(
    geoid_bg = str_pad(geoid_bg, 12, pad = "0"),
    geoid_tract = str_pad(geoid_tract, 11, pad = "0")
  )
# Provenance: guarantee change_lag_years exists and matches this
# script's change_lag setting (robust to older panel files that
# predate the provenance column; errors if the two scripts were
# run with different lags).
if ("change_lag_years" %in% names(acs_gent)) {
  stopifnot(
    "panel was built with a different change_lag" =
      all(acs_gent$change_lag_years == change_lag, na.rm = TRUE)
  )
}
acs_gent <- acs_gent |> mutate(change_lag_years = change_lag)

# ----------------------------
# Helpers (identical to main pipeline)
# ----------------------------
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
wins <- function(x, p = 0.01) {
  if (sum(!is.na(x)) < 2) return(x)
  qs <- quantile(x, c(p, 1 - p), na.rm = TRUE, names = FALSE)
  pmin(pmax(x, qs[1]), qs[2])
}
# cor() that returns NA instead of erroring when a year has no
# complete pairs (e.g. a variable is entirely NA in that year)
safe_cor <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok])
}
# Rank-based (Mann-Whitney) AUC: probability a randomly chosen
# positive outranks a randomly chosen negative. label must be 0/1.
auc_score <- function(score, label) {
  ok <- stats::complete.cases(score, label)
  s <- score[ok]
  l <- label[ok]
  if (length(l) < 3 || !all(l %in% c(0, 1))) return(NA_real_)
  # doubles, not integers: n1 * n0 overflows 32-bit ints when both
  # groups exceed ~46k (as they do in national BG samples)
  n1 <- as.numeric(sum(l == 1))
  n0 <- as.numeric(sum(l == 0))
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(s)
  (sum(r[l == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ============================================================
# 1. Tract-year indicator builder (common to both eras)
# ============================================================
# Expects a tibble with: geoid_tract, is_originated, is_denied,
# is_non_principal, loan_amount (dollars), applicant_income (dollars).
summarise_hmda <- function(x, yr) {
  x |>
    group_by(geoid_tract) |>
    summarise(
      hmda_apps_purchase = n(),
      hmda_originations_purchase = sum(is_originated, na.rm = TRUE),
      hmda_denials_purchase = sum(is_denied, na.rm = TRUE),
      hmda_denial_rate = hmda_denials_purchase / hmda_apps_purchase,
      hmda_origination_rate = hmda_originations_purchase / hmda_apps_purchase,
      hmda_non_principal_apps = sum(is_non_principal, na.rm = TRUE),
      hmda_non_principal_share = hmda_non_principal_apps / hmda_apps_purchase,
      hmda_median_applicant_income = median(applicant_income, na.rm = TRUE),
      hmda_median_loan_amount = median(loan_amount, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      year = yr,
      across(
        c(hmda_denial_rate, hmda_origination_rate, hmda_non_principal_share,
          hmda_median_applicant_income, hmda_median_loan_amount),
        ~ ifelse(is.nan(.x) | is.infinite(.x), NA_real_, .x)
      )
    ) |>
    relocate(year)
}

# ============================================================
# 2. Legacy era, 2012-2017: national LAR files (codes version)
# ============================================================
# URL pattern from consumerfinance.gov/data-research/hmda/historic-data
# If a download 404s, check that page for the current file naming.
process_hmda_legacy_year <- function(yr) {
  cache <- file.path(hmda_dir, paste0("hmda_tract_", yr, ".csv"))
  if (file.exists(cache)) {
    cached <- fread(cache, colClasses = c(geoid_tract = "character")) |> as_tibble()
    # schema check: caches written by older script versions (e.g.
    # pre-rename column names) are deleted and rebuilt, not trusted
    if ("hmda_non_principal_share" %in% names(cached)) {
      message("HMDA tract cache found: ", yr)
      return(cached)
    }
    message("Stale cache schema for ", yr, " - deleting and rebuilding.")
    file.remove(cache)
  }
  zip_file <- file.path(hmda_dir, paste0("hmda_", yr, "_nationwide_all-records_codes.zip"))
  # A zip is only trusted if it is plausibly large (legacy national
  # files run ~170-800 MB compressed) AND its central directory is
  # readable. Anything else is deleted and re-fetched, so a corrupt
  # cache can never poison later runs.
  zip_ok <- function(f) {
    file.exists(f) &&
      file.size(f) > 5e7 &&
      !inherits(try(unzip(f, list = TRUE), silent = TRUE), "try-error")
  }
  if (!zip_ok(zip_file)) {
    if (file.exists(zip_file)) {
      message("Corrupt/partial zip for ", yr, " - deleting and re-downloading.")
      file.remove(zip_file)
    }
    url <- paste0(
      "https://files.consumerfinance.gov/hmda-historic-loan-data/",
      "hmda_", yr, "_nationwide_all-records_codes.zip"
    )
    message("Downloading legacy HMDA ", yr, " (a few hundred MB) ...")
    tryCatch(
      download.file(url, zip_file, mode = "wb"),
      error = function(e) {
        if (file.exists(zip_file)) file.remove(zip_file) # no partial cache
        stop("legacy HMDA download failed for ", yr, ": ", e$message,
             "\nCheck filenames at consumerfinance.gov/data-research/hmda/historic-data",
             call. = FALSE)
      }
    )
    stopifnot(
      "downloaded zip is invalid - re-run, or check URL/network" = zip_ok(zip_file)
    )
  }
  message("Processing legacy HMDA ", yr, " ...")
  tmp <- tempdir()
  csv_inside <- unzip(zip_file, list = TRUE)$Name[1]
  unzip(zip_file, files = csv_inside, exdir = tmp, overwrite = TRUE)
  x <- fread(
    file.path(tmp, csv_inside),
    select = c("state_code", "county_code", "census_tract_number",
               "loan_purpose", "property_type", "owner_occupancy",
               "action_taken", "loan_amount_000s", "applicant_income_000s"),
    showProgress = FALSE
  )
  file.remove(file.path(tmp, csv_inside))
  out <- x |>
    as_tibble() |>
    filter(
      loan_purpose == 1,             # home purchase
      action_taken %in% 1:5,         # true applications only
      property_type %in% c(1, 2),    # 1-4 family + manufactured
      !is.na(state_code), !is.na(county_code), !is.na(census_tract_number)
    ) |>
    mutate(
      geoid_tract = paste0(
        sprintf("%02d", as.integer(state_code)),
        sprintf("%03d", as.integer(county_code)),
        # "0201.00" / "9601" both -> 6-digit tract code
        sprintf("%06d", round(suppressWarnings(as.numeric(census_tract_number)) * 100))
      ),
      is_originated = action_taken == 1,
      is_denied = action_taken == 3,
      is_non_principal = owner_occupancy == 2,
      loan_amount = suppressWarnings(as.numeric(loan_amount_000s)) * 1000,
      applicant_income = suppressWarnings(as.numeric(applicant_income_000s)) * 1000
    ) |>
    filter(str_detect(geoid_tract, "^\\d{11}$"),
           !str_sub(geoid_tract, 1, 2) %in% c("72", "66", "60", "69", "78")) |>
    summarise_hmda(yr)
  fwrite(out, cache)
  out
}

# ============================================================
# 3. New era, 2018-2019: data browser API, per state
# ============================================================
process_hmda_api_year <- function(yr) {
  cache <- file.path(hmda_dir, paste0("hmda_tract_", yr, ".csv"))
  if (file.exists(cache)) {
    cached <- fread(cache, colClasses = c(geoid_tract = "character")) |> as_tibble()
    # schema check: caches written by older script versions (e.g.
    # pre-rename column names) are deleted and rebuilt, not trusted
    if ("hmda_non_principal_share" %in% names(cached)) {
      message("HMDA tract cache found: ", yr)
      return(cached)
    }
    message("Stale cache schema for ", yr, " - deleting and rebuilding.")
    file.remove(cache)
  }
  per_state <- map_dfr(states, function(st) {
    raw_file <- file.path(hmda_dir, paste0("hmda_", st, "_", yr, ".csv"))
    if (!file.exists(raw_file)) {
      url <- paste0(
        "https://ffiec.cfpb.gov/v2/data-browser-api/view/csv?",
        "states=", st, "&years=", yr
      )
      message("Downloading HMDA: ", st, " ", yr)
      ok <- tryCatch({
        download.file(url, raw_file, mode = "wb", quiet = TRUE)
        TRUE
      }, error = function(e) {
        message("  failed: ", st, " ", yr, " | ", e$message)
        FALSE
      })
      if (!ok) return(tibble())
    }
    x <- fread(
      raw_file,
      select = c("census_tract", "loan_purpose", "action_taken",
                 "occupancy_type", "total_units", "loan_amount", "income"),
      # fread type-guesses per file; census_tract/total_units/income
      # come back integer64 in some states and character in others,
      # which breaks bind_rows. Force stable types up front.
      colClasses = c(census_tract = "character",
                     total_units = "character",
                     income = "character"),
      showProgress = FALSE
    )
    x |>
      as_tibble() |>
      filter(
        loan_purpose == 1,
        action_taken %in% 1:5,
        total_units %in% c("1", "2", "3", "4")
      ) |>
      mutate(
        geoid_tract = str_pad(census_tract, 11, pad = "0"),
        is_originated = action_taken == 1,
        is_denied = action_taken == 3,
        # occupancy_type: 1 principal, 2 second home, 3 investment.
        # (2,3) matches the legacy "not owner-occupied" definition.
        is_non_principal = occupancy_type %in% c(2, 3),
        loan_amount = suppressWarnings(as.numeric(loan_amount)),
        applicant_income = suppressWarnings(as.numeric(income)) * 1000
      ) |>
      filter(str_detect(geoid_tract, "^\\d{11}$")) |>
      # return ONLY the harmonized columns so per-state schema
      # quirks can never reach bind_rows
      select(geoid_tract, is_originated, is_denied, is_non_principal,
             loan_amount, applicant_income)
  })
  out <- summarise_hmda(per_state, yr)
  fwrite(out, cache)
  out
}

hmda_tract_year <- bind_rows(
  map_dfr(legacy_years, process_hmda_legacy_year),
  map_dfr(api_years, process_hmda_api_year)
) |>
  # reliability screen: tiny-application tracts -> NA indicators
  mutate(
    across(
      c(hmda_denial_rate, hmda_origination_rate, hmda_non_principal_share,
        hmda_median_applicant_income, hmda_median_loan_amount),
      ~ ifelse(hmda_apps_purchase < min_apps, NA_real_, .x)
    )
  )
fwrite(
  hmda_tract_year,
  file.path(hmda_dir, paste0("hmda_tract_year_indicators_",
                             min(hmda_years), "_", max(hmda_years), ".csv"))
)

# ============================================================
# 3b. Seasonal/second-home screen (ACS B25004, block group)
# ============================================================
# The non-principal-residence share cannot separate investors from
# vacation homes, so unscreened resort markets (Nantucket, the
# Hamptons) top the capital measure. This pulls seasonal-unit
# shares once and flags BGs above seasonal_share_max.
seasonal_file <- file.path(
  gent_dir, paste0("acs_bg_seasonal_", seasonal_ref_year, ".csv")
)
if (file.exists(seasonal_file) && file.size(seasonal_file) > 0) {
  seasonal_bg <- fread(seasonal_file, colClasses = c(geoid_bg = "character")) |>
    as_tibble()
} else {
  seasonal_bg <- map_dfr(states, function(st) {
    message("Downloading ACS seasonal-vacancy: ", st, " ", seasonal_ref_year)
    tryCatch(
      get_acs(
        geography = "block group",
        variables = c(total_units = "B25002_001",
                      seasonal_units = "B25004_006"),
        state = st, year = seasonal_ref_year, survey = "acs5",
        output = "wide", geometry = FALSE
      ) |>
        transmute(
          geoid_bg = GEOID,
          seasonal_share = ifelse(total_unitsE > 0,
                                  seasonal_unitsE / total_unitsE, NA_real_)
        ),
      error = function(e) {
        message("  failed: ", st, " | ", e$message)
        tibble()
      }
    )
  })
  if (nrow(seasonal_bg) > 0) fwrite(seasonal_bg, seasonal_file)
}
seasonal_bg <- seasonal_bg |>
  mutate(geoid_bg = str_pad(geoid_bg, 12, pad = "0")) |>
  distinct(geoid_bg, .keep_all = TRUE)

# ============================================================
# 4. Join to the BG panel; score consistently with Measures A-C
# ============================================================
# Income-to-AMI proxy uses the panel's own county median income
# (already computed per county-year in the main pipeline).
acs_capital <- acs_gent |>
  left_join(hmda_tract_year, by = c("year", "geoid_tract")) |>
  left_join(seasonal_bg, by = "geoid_bg") |>
  mutate(
    # BGs without seasonal data are treated as non-seasonal
    seasonal_market = replace_na(seasonal_share > seasonal_share_max, FALSE),
    hmda_income_to_county_ami =
      hmda_median_applicant_income / county_med_income,
    hmda_income_to_county_ami =
      ifelse(is.infinite(hmda_income_to_county_ami), NA_real_,
             hmda_income_to_county_ami)
  ) |>
  group_by(year) |>
  mutate(
    # winsorized inputs, consistent with Measure A
    hmda_non_principal_share_w = wins(hmda_non_principal_share),
    hmda_income_to_ami_w = wins(hmda_income_to_county_ami),
    hmda_log_loan_w = wins(log1p(hmda_median_loan_amount)),
    hmda_denial_rate_w = wins(hmda_denial_rate),
    # z-scores among ELIGIBLE BGs only, like the transition score
    z_hmda_non_principal = z_eligible(hmda_non_principal_share_w, eligible_gentrifiable),
    z_hmda_income_to_ami = z_eligible(hmda_income_to_ami_w, eligible_gentrifiable),
    z_hmda_loan = z_eligible(hmda_log_loan_w, eligible_gentrifiable),
    z_hmda_denial = z_eligible(hmda_denial_rate_w, eligible_gentrifiable),
    # Measure D: capital flowing in = non-owner + high-income
    # applicants + big loans + easy credit (low denial)
    capital_reinvestment_score =
      z_hmda_non_principal + z_hmda_income_to_ami + z_hmda_loan - z_hmda_denial,
    capital_reinvestment_score_pctile =
      pct_rank_eligible(capital_reinvestment_score, eligible_gentrifiable),
    gentrifying_capital_top10 = as.integer(
      eligible_gentrifiable == TRUE &
        !is.na(capital_reinvestment_score_pctile) &
        capital_reinvestment_score_pctile >= 0.90
    ),
    gentrifying_capital_top05 = as.integer(
      eligible_gentrifiable == TRUE &
        !is.na(capital_reinvestment_score_pctile) &
        capital_reinvestment_score_pctile >= 0.95
    ),
    # SCREENED variant: seasonal/resort markets are removed from
    # the percentile pool AND cannot be flagged. The score itself
    # is unchanged - only the ranking universe differs.
    capital_pctile_screened = pct_rank_eligible(
      capital_reinvestment_score,
      eligible_gentrifiable & !seasonal_market
    ),
    gentrifying_capital_top10_screened = as.integer(
      eligible_gentrifiable == TRUE & !seasonal_market &
        !is.na(capital_pctile_screened) &
        capital_pctile_screened >= 0.90
    ),
    gentrifying_capital_top05_screened = as.integer(
      eligible_gentrifiable == TRUE & !seasonal_market &
        !is.na(capital_pctile_screened) &
        capital_pctile_screened >= 0.95
    )
  ) |>
  ungroup()

# ============================================================
# 5. Optional HUD-USPS vacancy loader (unchanged; restricted data)
# ============================================================
hud_raw_dir <- file.path(hud_dir, "raw")
dir.create(hud_raw_dir, recursive = TRUE, showWarnings = FALSE)
hud_files <- list.files(hud_raw_dir, pattern = "\\.csv$", full.names = TRUE)
process_hud_usps_file <- function(file) {
  message("Processing HUD-USPS vacancy: ", basename(file))
  x <- fread(file, showProgress = FALSE) |> clean_names()
  yr <- str_extract(basename(file), "20\\d{2}") |> as.integer()
  x <- x |>
    rename_with(~ "geoid_tract", matches("tract|geoid|fips")) |>
    mutate(
      geoid_tract = str_pad(str_remove_all(as.character(geoid_tract), "\\."), 11, pad = "0"),
      year = yr
    )
  res_total <- intersect(names(x), c("residential_addresses", "res_total", "resaddr", "total_residential"))
  res_vac <- intersect(names(x), c("residential_vacant", "res_vacant", "vacant_residential", "resvac"))
  if (length(res_total) == 0 || length(res_vac) == 0) {
    message("Could not identify vacancy columns; names are: ",
            paste(names(x), collapse = ", "))
    return(tibble())
  }
  x |>
    transmute(
      year, geoid_tract,
      hud_res_total = as.numeric(.data[[res_total[1]]]),
      hud_res_vacant = as.numeric(.data[[res_vac[1]]]),
      hud_res_vacancy_rate = hud_res_vacant / hud_res_total
    )
}
if (length(hud_files) > 0) {
  hud_usps_year <- map_dfr(hud_files, process_hud_usps_file) |>
    group_by(year, geoid_tract) |>
    summarise(across(starts_with("hud_"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  acs_capital <- acs_capital |>
    left_join(hud_usps_year, by = c("year", "geoid_tract"))
} else {
  message("No HUD-USPS files in ", hud_raw_dir, " - skipping vacancy join.")
}

# ============================================================
# 6. Save (provenance column included, same naming scheme)
# ============================================================
capital_file <- file.path(
  gent_dir,
  paste0("capital_reinvestment_gentrification_2009_2019_", lag_tag, ".csv")
)
acs_capital |>
  select(
    year, geoid_bg, geoid_tract, statefp, countyfp, state_county,
    bg_vintage, change_lag_years, analysis_sample_main, eligible_gentrifiable,
    hmda_apps_purchase, hmda_denial_rate, hmda_origination_rate,
    hmda_non_principal_share, hmda_median_applicant_income,
    hmda_median_loan_amount, hmda_income_to_county_ami,
    seasonal_share, seasonal_market,
    capital_reinvestment_score, capital_reinvestment_score_pctile,
    gentrifying_capital_top10, gentrifying_capital_top05,
    capital_pctile_screened,
    gentrifying_capital_top10_screened, gentrifying_capital_top05_screened
  ) |>
  fwrite(capital_file)
combined_file <- file.path(
  gent_dir,
  paste0("blockgroup_gentrification_2009_2019_", lag_tag, "_with_HMDA.csv")
)
fwrite(acs_capital, combined_file)
message("Saved capital reinvestment file to: ", capital_file)
message("Saved combined panel to: ", combined_file)

# ============================================================
# DIAGNOSTICS
# ============================================================
options(pillar.sigfig = 4)

# --- 1. HMDA raw coverage by year ----
diag_hmda_cov <- hmda_tract_year |>
  group_by(year) |>
  summarise(
    n_tracts = n(),
    total_apps = sum(hmda_apps_purchase),
    pct_tracts_below_min_apps = mean(hmda_apps_purchase < min_apps),
    .groups = "drop"
  )
message("\n---- 1. HMDA tract coverage by year ----")
print(diag_hmda_cov, n = Inf, width = Inf)

# --- 2. Legacy/new schema seam check (2017 vs 2018) ----
diag_seam <- hmda_tract_year |>
  group_by(year) |>
  summarise(
    mean_denial_rate = mean(hmda_denial_rate, na.rm = TRUE),
    mean_non_principal_share = mean(hmda_non_principal_share, na.rm = TRUE),
    med_loan_amount = median(hmda_median_loan_amount, na.rm = TRUE),
    med_applicant_income = median(hmda_median_applicant_income, na.rm = TRUE),
    .groups = "drop"
  )
message("\n---- 2. Era-seam check (watch for a jump between 2017 and 2018) ----")
print(diag_seam, n = Inf, width = Inf)
fwrite(diag_seam, file.path(tab_dir, "hmda_seam_check_by_year.csv"))

# --- 3. Join coverage into the BG panel ----
diag_join <- acs_capital |>
  filter(year %in% hmda_years) |>
  group_by(year) |>
  summarise(
    n_bg = n(),
    pct_bg_with_hmda = mean(!is.na(hmda_apps_purchase)),
    pct_bg_with_score = mean(!is.na(capital_reinvestment_score)),
    pct_elig_with_score =
      mean(!is.na(capital_reinvestment_score[eligible_gentrifiable == TRUE])),
    .groups = "drop"
  )
message("\n---- 3. HMDA join coverage in BG panel by year ----")
print(diag_join, n = Inf, width = Inf)
fwrite(diag_join, file.path(tab_dir, paste0("hmda_join_coverage_", lag_tag, ".csv")))

# --- 4. Capital score distribution by year ----
diag_cap_dist <- acs_capital |>
  filter(analysis_sample_main, !is.na(capital_reinvestment_score)) |>
  group_by(year) |>
  summarise(
    n = n(),
    p05 = quantile(capital_reinvestment_score, 0.05),
    p50 = quantile(capital_reinvestment_score, 0.50),
    p95 = quantile(capital_reinvestment_score, 0.95),
    p99 = quantile(capital_reinvestment_score, 0.99),
    max = max(capital_reinvestment_score),
    .groups = "drop"
  )
message("\n---- 4. Capital score distribution by year ----")
print(diag_cap_dist, n = Inf, width = Inf)

# --- 5. Flag rates by year ----
diag_cap_flags <- acs_capital |>
  filter(analysis_sample_main) |>
  group_by(year) |>
  summarise(
    n_eligible = sum(eligible_gentrifiable == TRUE, na.rm = TRUE),
    pct_capital_top10 = mean(gentrifying_capital_top10 == 1, na.rm = TRUE),
    pct_capital_top10_elig =
      mean(gentrifying_capital_top10[eligible_gentrifiable == TRUE] == 1, na.rm = TRUE),
    pct_capital_top05_elig =
      mean(gentrifying_capital_top05[eligible_gentrifiable == TRUE] == 1, na.rm = TRUE),
    .groups = "drop"
  )
message("\n---- 5. Capital flag rates by year ----")
print(diag_cap_flags, n = Inf, width = Inf)
fwrite(diag_cap_flags, file.path(tab_dir, paste0("hmda_capital_flags_", lag_tag, ".csv")))

# --- 6. Agreement with Measures A-C ----
diag_cap_agree <- acs_capital |>
  filter(analysis_sample_main) |>
  group_by(year) |>
  summarise(
    cor_capital_transition = safe_cor(capital_reinvestment_score, gent_score_transition),
    cor_capital_expectation = safe_cor(capital_reinvestment_score, expectation_gap),
    jac_capital_transition =
      jaccard(gentrifying_capital_top10, gentrifying_transition_top10),
    jac_capital_expectation =
      jaccard(gentrifying_capital_top10, gentrifying_expectation_top10),
    jac_capital_classic =
      jaccard(gentrifying_capital_top10, gentrifying_classic_binary),
    .groups = "drop"
  )
message("\n---- 6. Agreement of capital measure with Measures A-C ----")
print(diag_cap_agree, n = Inf, width = Inf)
fwrite(diag_cap_agree, file.path(tab_dir, paste0("hmda_capital_agreement_", lag_tag, ".csv")))

# --- 6b. Pooled agreement (one sentence for the paper) ----
# Note: the capital-transition correlation trends upward over the
# decade (recovery-era noise in early baselines), so report the
# range alongside the mean.
diag_cap_agree_pooled <- diag_cap_agree |>
  summarise(
    mean_cor_transition = mean(cor_capital_transition, na.rm = TRUE),
    min_cor_transition = min(cor_capital_transition, na.rm = TRUE),
    max_cor_transition = max(cor_capital_transition, na.rm = TRUE),
    mean_cor_expectation = mean(cor_capital_expectation, na.rm = TRUE),
    mean_jaccard_transition = mean(jac_capital_transition, na.rm = TRUE)
  )
message("\n---- 6b. Pooled agreement across analysis years ----")
print(diag_cap_agree_pooled, width = Inf)
fwrite(diag_cap_agree_pooled,
       file.path(tab_dir, paste0("hmda_capital_agreement_pooled_", lag_tag, ".csv")))

# --- 6c. Seasonal/second-home screen effect ----
diag_seasonal <- acs_capital |>
  filter(analysis_sample_main) |>
  summarise(
    pct_bg_seasonal = mean(seasonal_market),
    pct_elig_seasonal =
      mean(seasonal_market[eligible_gentrifiable == TRUE], na.rm = TRUE),
    n_top10_unscreened = sum(gentrifying_capital_top10 == 1, na.rm = TRUE),
    n_top10_screened = sum(gentrifying_capital_top10_screened == 1, na.rm = TRUE),
    # share of UNSCREENED top-10 flags that sit in resort markets -
    # the size of the problem the screen removes
    pct_unscreened_flags_seasonal =
      mean(seasonal_market[gentrifying_capital_top10 == 1], na.rm = TRUE),
    jac_screened_unscreened =
      jaccard(gentrifying_capital_top10, gentrifying_capital_top10_screened)
  )
message("\n---- 6c. Seasonal screen effect (threshold ", seasonal_share_max, ") ----")
print(diag_seasonal, width = Inf)
fwrite(diag_seasonal, file.path(tab_dir, paste0("hmda_seasonal_screen_", lag_tag, ".csv")))

# --- 7. Spot-check: top-10 SCREENED capital-score BGs, 2019 ----
# (Nantucket/Hamptons-type rows should now be gone; urban
# gentrification BGs should remain.)
message("\n---- 7. Top-10 screened capital-score BGs, 2019 ----")
acs_capital |>
  filter(year == 2019, gentrifying_capital_top10_screened == 1) |>
  arrange(desc(capital_reinvestment_score)) |>
  select(geoid_bg, state_county, capital_reinvestment_score,
         seasonal_share, hmda_apps_purchase, hmda_non_principal_share,
         hmda_denial_rate, hmda_median_loan_amount, hmda_income_to_county_ami) |>
  head(10) |>
  print(width = Inf)

# --- 8. Missingness summary (for the Methods section) ----
diag_missing_cap <- acs_capital |>
  filter(analysis_sample_main) |>
  summarise(
    pct_missing_score = mean(is.na(capital_reinvestment_score)),
    pct_missing_income = mean(is.na(hmda_income_to_county_ami)),
    pct_missing_non_principal = mean(is.na(hmda_non_principal_share)),
    # among eligible BGs only - the denominator that matters for flags
    pct_missing_score_eligible =
      mean(is.na(capital_reinvestment_score[eligible_gentrifiable == TRUE]))
  )
message("\n---- 8. Missingness in analysis sample (Methods numbers) ----")
print(diag_missing_cap, width = Inf)
fwrite(diag_missing_cap,
       file.path(tab_dir, paste0("hmda_missingness_", lag_tag, ".csv")))

# ============================================================
# 9. VALIDATION SUMMARY (Table 2 scaffold)
# ============================================================
# One row per measure: validation against the Pratt metro-tract
# gentrification scores, plus the mean within-year correlation
# with the capital reinvestment score.
#
# Pratt file: metrotracts_gentscores_2010tr.csv
#   * 2010-VINTAGE tracts - matches the panel's geography, so a
#     plain GEOID join is valid (do NOT use the 2020tr files).
#   * id column tr2010gj is an NHGIS GISJOIN
#     ("G" + SS + "0" + CCC + "0" + TTTTTT) -> 11-digit GEOID.
#   * pratt_r   = correlation with GentIntensity_2010to2020, the
#     continuous factor-score intensity of 2010->2020 change.
#   * pratt_auc = AUC predicting Gentrifying_10to20 (0/1: tract
#     factor score above its metro average).
# Table 2 caption notes: Pratt covers METRO tracts only, and its
# window is the full 2010->2020 decade vs our lag-3 changes; BGs
# inherit their tract's Pratt values in the pratt_year cross-section.
pratt_dir <- file.path(base_dir, "Pratt")
pratt_file <- file.path(pratt_dir, "metrotracts_gentscores_2010tr.csv")
pratt_year <- 2019

measures <- tribble(
  ~measure,                        ~score_col,
  "Neighborhood transition (A)",   "gent_score_transition",
  "Market expectations (B)",       "expectation_gap",
  "Classic binary (C)",            "gentrifying_classic_binary",
  "Capital reinvestment (D)",      "capital_reinvestment_score"
)
is_binary_measure <- measures$score_col == "gentrifying_classic_binary"

# mean within-year correlation with the capital score (analysis years)
mean_cor_with_capital <- map_dbl(measures$score_col, function(v) {
  acs_capital |>
    filter(analysis_sample_main) |>
    group_by(year) |>
    summarise(r = safe_cor(.data[[v]], capital_reinvestment_score),
              .groups = "drop") |>
    pull(r) |>
    mean(na.rm = TRUE)
})
mean_cor_with_capital[measures$score_col == "capital_reinvestment_score"] <- NA_real_

pratt_r <- rep(NA_real_, nrow(measures))
pratt_r_n <- rep(NA_real_, nrow(measures))
pratt_auc <- rep(NA_real_, nrow(measures))
pratt_r_tract <- rep(NA_real_, nrow(measures))
pratt_auc_tract <- rep(NA_real_, nrow(measures))
classic_sens <- NA_real_
classic_spec <- NA_real_
if (file.exists(pratt_file)) {
  pratt <- fread(
    pratt_file,
    select = c("tr2010gj", "Gentrifying_10to20", "Gentrified_10to20",
               "GentIntensity_2010to2020")
  ) |>
    as_tibble() |>
    mutate(
      d = str_remove_all(tr2010gj, "\\D"), # 13 digits: SS 0 CCC 0 TTTTTT
      geoid_tract = paste0(str_sub(d, 1, 2), str_sub(d, 4, 6), str_sub(d, 8, 13)),
      pratt_intensity = suppressWarnings(as.numeric(GentIntensity_2010to2020)),
      pratt_gentrifying = suppressWarnings(as.numeric(Gentrifying_10to20)),
      pratt_gentrified = suppressWarnings(as.numeric(Gentrified_10to20))
    ) |>
    filter(nchar(d) == 13, str_detect(geoid_tract, "^\\d{11}$")) |>
    distinct(geoid_tract, .keep_all = TRUE) |>
    select(geoid_tract, pratt_intensity, pratt_gentrifying, pratt_gentrified)
  val <- acs_capital |>
    filter(year == pratt_year) |>
    inner_join(pratt, by = "geoid_tract")
  message("Pratt validation: ", nrow(val), " BG rows matched ",
          "(metro tracts only, ", pratt_year, " cross-section)")
  pratt_r <- map_dbl(measures$score_col,
                     ~ safe_cor(val[[.x]], val$pratt_intensity))
  # n behind each correlation, for Fisher-z CIs in the figures
  pratt_r_n <- map_dbl(measures$score_col,
                       ~ sum(stats::complete.cases(val[[.x]], val$pratt_intensity)))
  pratt_auc <- map_dbl(measures$score_col,
                       ~ auc_score(val[[.x]], val$pratt_gentrifying))
  # secondary benchmark: AUC vs the stricter class-upgrade flag
  pratt_auc_gentrified <- map_dbl(measures$score_col,
                                  ~ auc_score(val[[.x]], val$pratt_gentrified))
  message("AUC vs Gentrified_10to20 (stricter class-upgrade flag): ",
          paste(measures$measure, round(pratt_auc_gentrified, 3),
                sep = " = ", collapse = "; "))
  # Classic binary is a 0/1 predictor: AUC is not meaningful, so
  # report sensitivity/specificity against Pratt instead.
  classic_sens <- mean(
    val$gentrifying_classic_binary[val$pratt_gentrifying == 1] == 1,
    na.rm = TRUE
  )
  classic_spec <- mean(
    val$gentrifying_classic_binary[val$pratt_gentrifying == 0] == 0,
    na.rm = TRUE
  )
  # Tract-collapsed robustness: BGs within a tract share their
  # Pratt values (and the tract-level capital score), so BG rows
  # are not independent. Repeat validation on one row per tract,
  # averaging BG-level measure scores within tract.
  val_tract <- val |>
    group_by(geoid_tract) |>
    summarise(
      across(all_of(measures$score_col), ~ mean(.x, na.rm = TRUE)),
      pratt_intensity = first(pratt_intensity),
      pratt_gentrifying = first(pratt_gentrifying),
      .groups = "drop"
    )
  message("Tract-collapsed validation: ", nrow(val_tract), " tracts")
  pratt_r_tract <- map_dbl(measures$score_col,
                           ~ safe_cor(val_tract[[.x]], val_tract$pratt_intensity))
  pratt_auc_tract <- map_dbl(measures$score_col,
                             ~ auc_score(val_tract[[.x]], val_tract$pratt_gentrifying))
} else {
  message("Pratt file not found: ", pratt_file, " - Pratt columns left NA.")
}

# AUC is suppressed for the binary measure (not meaningful for a
# 0/1 predictor); it gets sensitivity/specificity instead. The
# rank column orders measures by alignment with Pratt intensity.
validation_summary <- measures |>
  transmute(
    measure,
    pratt_r = pratt_r,
    pratt_r_n = pratt_r_n,
    pratt_auc = ifelse(is_binary_measure, NA_real_, pratt_auc),
    pratt_sensitivity = ifelse(is_binary_measure, classic_sens, NA_real_),
    pratt_specificity = ifelse(is_binary_measure, classic_spec, NA_real_),
    pratt_r_tract = pratt_r_tract,
    pratt_auc_tract = ifelse(is_binary_measure, NA_real_, pratt_auc_tract),
    mean_cor_with_capital = mean_cor_with_capital,
    rank_pratt_r = dense_rank(desc(pratt_r))
  )
message("\n---- 9. Validation summary (Table 2) ----")
print(validation_summary, width = Inf)
fwrite(validation_summary,
       file.path(tab_dir, paste0("validation_summary_table2_", lag_tag, ".csv")))

# ============================================================
# 10. LAG SENSITIVITY (lag-3 vs lag-5)
# ============================================================
# Lag 5 uses fully non-overlapping ACS windows. To generate it:
#   1. Set change_lag <- 5 in blockgroup_gentrification_earliest_acs.R
#      and re-run (analysis window becomes 2015-2019).
#   2. Set change_lag <- 5 at the top of THIS script and re-run.
# When both lag versions of the transition file exist, this block
# compares them on the overlapping years automatically.
other_lag <- if (change_lag == 3) 5 else 3
this_transition_file <- file.path(
  gent_dir, paste0("transition_gentrification_2009_2019_lag", change_lag, ".csv")
)
other_transition_file <- file.path(
  gent_dir, paste0("transition_gentrification_2009_2019_lag", other_lag, ".csv")
)
if (file.exists(this_transition_file) && file.exists(other_transition_file)) {
  read_transition <- function(f) {
    fread(f, colClasses = c(geoid_bg = "character"),
          select = c("year", "geoid_bg", "gent_score_transition",
                     "gentrifying_transition_top10")) |>
      as_tibble()
  }
  lag_compare <- inner_join(
    read_transition(this_transition_file),
    read_transition(other_transition_file),
    by = c("year", "geoid_bg"),
    suffix = c("_a", "_b") # a = this lag, b = other lag
  ) |>
    group_by(year) |>
    summarise(
      n_both = sum(stats::complete.cases(gent_score_transition_a,
                                         gent_score_transition_b)),
      cor_scores = safe_cor(gent_score_transition_a, gent_score_transition_b),
      jac_top10 = jaccard(gentrifying_transition_top10_a,
                          gentrifying_transition_top10_b),
      .groups = "drop"
    ) |>
    filter(!is.na(cor_scores))
  message("\n---- 10. Transition measure, lag-", change_lag, " (a) vs lag-",
          other_lag, " (b), overlapping years ----")
  print(lag_compare, n = Inf, width = Inf)
  fwrite(lag_compare,
         file.path(tab_dir, paste0("lag_sensitivity_transition_lag",
                                   change_lag, "_vs_lag", other_lag, ".csv")))
} else {
  message("\n---- 10. Lag sensitivity skipped ----\n",
          "Re-run the main pipeline and this script with change_lag <- ",
          other_lag, "; this block will then compare the two lag versions.")
}

message("\nAll HMDA diagnostics saved to: ", tab_dir)
