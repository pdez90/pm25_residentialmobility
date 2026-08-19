# ============================================================
# IDEA 2 (PM2.5 ONLY): DOES MIGRATION EQUALIZE OR ENTRENCH
#         PM2.5 EXPOSURE?  The pollution-avoidance privilege.
# ------------------------------------------------------------
# Conditional on moving, do residents of disadvantaged and
# majority-minority neighborhoods achieve SMALLER PM2.5 reductions
# than those of advantaged/white neighborhoods, and does
# gentrification-driven displacement blunt pollution avoidance?
# Anchored on Keller et al. (NBER w35317); extended to block-group
# destination sorting and the exposure-reallocation health framework.
#
# Estimand: flow-weighted destination-minus-origin PM2.5 change for
# out-movers, aggregated to the origin-BG-year, compared across origin
# composition (race, income) and gentrification status via flow-weighted
# WLS with county + year FE and county-clustered SEs. Negative change =
# reduction = avoidance. An exposure-reallocation signature by subgroup is
# reported for context; deliberately NOT converted to averted deaths (see below).
#
# RACE CAVEAT: MIGRATE flows are not disaggregated by race, so the race
# contrast is ECOLOGICAL, at the level of ORIGIN-neighborhood composition
# (movers FROM majority-Black/Hispanic/white origins), not individuals.
#
# Input : pollution_bg/bg_pm25_panel.csv  (from extract_pm25_bg.R)
# Output: /Users/priyanka/MIGRATE_dir/pm25_avoidance_<tag>/
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr) })
if (!requireNamespace("fixest", quietly = TRUE)) install.packages("fixest")
library(fixest)

# ---------------------------
# Config
# ---------------------------
base_dir <- "/Users/priyanka/MIGRATE_dir"; mig_dir <- base_dir
gent_dir <- file.path(base_dir, "acs_gentrification")
change_lag <- 3; origin_lag <- 1
tag <- paste0("lag", change_lag, "_olag", origin_lag)
out_dir <- file.path(base_dir, paste0("pm25_avoidance_", tag)); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
migrate_years <- 2011:2019          # full available MIGRATE flow range (2006-2010 files absent)
panel_file <- file.path(base_dir, "pollution_bg", "bg_pm25_panel.csv")
combined_file <- file.path(gent_dir, paste0("blockgroup_gentrification_2009_2019_lag", change_lag, "_with_HMDA.csv"))

RACE_PULL <- TRUE                                # pull ACS race if race_bg.csv absent
CENSUS_KEY <- Sys.getenv("CENSUS_API_KEY")
# NOTE: the PM2.5 mortality CRF and averted-death conversion were removed; the
# design (summed multi-year move-events, adjusted group contrast, no exposure
# duration/age-specific rates) does not support attributable-mortality claims.

pad_bg <- function(x) str_pad(as.character(x), 12, pad = "0")
wmean <- function(x, w){ ok <- !is.na(x) & !is.na(w); if (!any(ok)) return(NA_real_); sum(x[ok]*w[ok])/sum(w[ok]) }

# ---------------------------
# 1. PM2.5 panel
# ---------------------------
poll <- fread(panel_file, colClasses = c(geoid_bg = "character")); poll[, geoid_bg := pad_bg(geoid_bg)]
stopifnot("pm25 missing from panel" = "pm25" %in% names(poll))

# ---------------------------
# 2. Origin composition: income/poverty + gentrification + race (pull)
# ---------------------------
hdr <- names(fread(combined_file, nrows = 0))
cap_col <- if ("gentrifying_capital_top10_screened" %in% hdr) "gentrifying_capital_top10_screened" else "gentrifying_capital_top10"
acs <- fread(combined_file, colClasses = c(geoid_bg = "character", state_county = "character"),
  select = unique(c("year","geoid_bg","state_county","eligible_gentrifiable",
                    "gentrifying_transition_top10", "med_hh_income","poverty_pct","bach_pct","total_pop")))
acs[, geoid_bg := pad_bg(geoid_bg)]
acs[, inc_tercile := cut(med_hh_income, quantile(med_hh_income, c(0,1/3,2/3,1), na.rm=TRUE),
                         labels = c("Low income","Middle","High income"), include.lowest = TRUE), by = year]
acs[, high_poverty := as.integer(poverty_pct > median(poverty_pct, na.rm=TRUE)), by = year]
acs[, gent_group := fifelse(gentrifying_transition_top10 == 1, "Gentrifying",
        fifelse(gentrifying_transition_top10 == 0 & eligible_gentrifiable == TRUE, "Eligible, not gentrifying", NA_character_))]

# Race is matched to each flow YEAR from its own ACS 5-year vintage:
# NHGIS block-group B03002 for endyears <=2012 (2011, 2012) and
# tidycensus for 2013+. Each vintage is loaded once and cached.
race_dir <- file.path(base_dir, "pollution_bg"); dir.create(race_dir, showWarnings = FALSE, recursive = TRUE)
# True year match: NHGIS block-group ACS5 covers endyears <=2012 (2011 -> 2007-2011,
# 2012 -> 2008-2012); tidycensus covers 2013+. Each flow year uses its own vintage.
vintage_of <- function(y) as.integer(y)
nhgis_zip <- { z <- list.files(file.path(base_dir,"nhgis"), pattern="nhgis.*_csv\\.zip$", full.names=TRUE)
               if (length(z)) z[order(file.mtime(z), decreasing=TRUE)][1] else NA_character_ }

# NHGIS loader for endyears <=2012: reads B03002 straight from the extract,
# parsing the per-dataset NHGIS code prefix (e.g. MN2 in 2011, QSY in 2012).
load_race_nhgis <- function(v) {
  if (is.na(nhgis_zip)) return(NULL)
  members <- unzip(nhgis_zip, list=TRUE)$Name
  tok <- sprintf("_20%02d5_blck_grp", v %% 100)                 # 2011 -> _20115_blck_grp
  csvm <- members[grepl(paste0(tok,"\\.csv$"), members)]
  cbm  <- members[grepl(paste0(tok,"_codebook\\.txt$"), members)]
  if (!length(csvm) || !length(cbm)) return(NULL)
  cb <- readLines(unz(nhgis_zip, cbm[1]), warn=FALSE)
  si <- grep("Source code:\\s*B03002", cb); nc <- grep("NHGIS code:", cb)
  if (!length(si) || !length(nc)) return(NULL)
  pref <- trimws(sub(".*NHGIS code:\\s*", "", cb[nc[nc > si[1]][1]]))   # e.g. "MN2"
  cc <- function(n) paste0(pref, "E", n)
  vcols <- c(cc("001"), cc("003"), cc("004"), cc("012"))
  tmp <- tempfile("nhgis_"); dir.create(tmp)
  unzip(nhgis_zip, files=csvm[1], exdir=tmp, junkpaths=TRUE)
  f <- list.files(tmp, pattern="\\.csv$", full.names=TRUE)[1]
  # NHGIS CSVs carry a SECOND header row of labels; read all as character,
  # drop that non-numeric row, then coerce the value columns.
  d <- fread(f, colClasses="character",
             select=c("STATEA","COUNTYA","TRACTA","BLKGRPA", vcols))
  d <- d[!is.na(suppressWarnings(as.numeric(get(cc("001")))))]
  d[, (vcols) := lapply(.SD, as.numeric), .SDcols=vcols]
  d[, geoid_bg := paste0(str_pad(STATEA,2,pad="0"), str_pad(COUNTYA,3,pad="0"),
                         str_pad(TRACTA,6,pad="0"), str_pad(BLKGRPA,1,pad="0"))]
  setnames(d, vcols, c("tot","white","black","hisp"))
  d[tot > 0, `:=`(sh_white=white/tot, sh_black=black/tot, sh_hisp=hisp/tot)]
  d[, .(geoid_bg, sh_white, sh_black, sh_hisp)]
}

pull_race_vintage <- function(v) {
  cf <- file.path(race_dir, paste0("race_bg_", v, ".csv"))
  if (file.exists(cf)) { r <- fread(cf, colClasses = c(geoid_bg="character")); r[, geoid_bg := pad_bg(geoid_bg)]; return(r) }
  r <- NULL
  if (v <= 2012) {
    r <- load_race_nhgis(v)
    if (is.null(r)) message("NHGIS race for ", v, " not found in nhgis/ extract; run nhgis_download.R (B03002).")
  } else if (RACE_PULL && requireNamespace("tidycensus", quietly = TRUE) && nzchar(CENSUS_KEY)) {
    library(tidycensus); census_api_key(CENSUS_KEY, install = FALSE)
    vars <- c(tot="B03002_001", white="B03002_003", black="B03002_004", hisp="B03002_012")
    sts <- unique(substr(acs$geoid_bg, 1, 2))
    message("Pulling ACS race (B03002) block groups for ", v, " via tidycensus ...")
    rl <- lapply(sts, function(st) tryCatch(
      get_acs(geography="block group", variables=vars, state=st, year=v, survey="acs5", output="wide"),
      error = function(e) NULL))
    rr <- as.data.table(do.call(rbind, Filter(Negate(is.null), rl)))
    if (nrow(rr)) { rr[, geoid_bg := pad_bg(GEOID)]
      rr[, `:=`(sh_white=whiteE/totE, sh_black=blackE/totE, sh_hisp=hispE/totE)]
      r <- rr[, .(geoid_bg, sh_white, sh_black, sh_hisp)] }
  }
  if (!is.null(r)) { r[, geoid_bg := pad_bg(geoid_bg)]; fwrite(r, cf) }
  r
}
race <- rbindlist(lapply(migrate_years, function(y){
  rv <- pull_race_vintage(vintage_of(y)); if (is.null(rv)) return(NULL)
  copy(rv)[, year := y]
}), fill = TRUE)
race <- if (nrow(race)) race else NULL
if (!is.null(race)) {
  race[, race_majority := fifelse(sh_black > 0.5, "Majority Black",
              fifelse(sh_hisp > 0.5, "Majority Hispanic",
              fifelse(sh_white > 0.5, "Majority white (NH)", "No majority")))]
  acs <- merge(acs, race[, .(year, geoid_bg, race_majority)], by = c("year","geoid_bg"), all.x = TRUE)
} else message("Race not available (no key/tidycensus/cache). Income/poverty subgroups only.")
acs <- merge(acs, poll, by = c("geoid_bg","year"), all.x = TRUE)

# ---------------------------
# 3. Origin / destination lookups
# ---------------------------
sub_cols <- c("inc_tercile","high_poverty","gent_group", if (!is.null(race)) "race_majority")
o_keep <- c("year","geoid_bg","state_county","med_hh_income","poverty_pct","bach_pct","total_pop", sub_cols, "pm25")
oa <- acs[, ..o_keep]; setnames(oa, setdiff(o_keep, c("year","geoid_bg")), paste0("o_", setdiff(o_keep, c("year","geoid_bg"))))
setnames(oa, c("year","geoid_bg"), c("cy","origin_bg")); setkey(oa, cy, origin_bg)
da <- acs[, .(year, geoid_bg, pm25)]; setnames(da, c("pm25","year","geoid_bg"), c("d_pm25","cy","dest_bg")); setkey(da, cy, dest_bg)

# ---------------------------
# 4. MIGRATE pass: flow-weighted PM2.5 reallocation to origin-BG-year
# ---------------------------
acc <- list()
for (yr in migrate_years) {
  x <- fread(file.path(mig_dir, paste0(yr, ".csv")), select = c("Origin","Destination","Flow"))
  x <- x[, .(origin_bg = pad_bg(Origin), dest_bg = pad_bg(Destination), flow = as.numeric(Flow))][flow > 0 & origin_bg != dest_bg]
  x[, cy := yr]
  setkey(x, cy, origin_bg); x <- oa[x, on = .(cy, origin_bg)]
  setkey(x, cy, dest_bg);   x <- da[x, on = .(cy, dest_bg)]
  x[, pm25_re := d_pm25 - o_pm25]
  byv <- c("origin_bg","o_state_county","o_med_hh_income","o_poverty_pct","o_bach_pct", paste0("o_", sub_cols))
  acc[[as.character(yr)]] <- x[, .(year = yr, out_flow = sum(flow, na.rm=TRUE),
                                   re_pm25 = wmean(pm25_re, flow), orig_pm25 = wmean(o_pm25, flow)), by = byv]
}
oy <- rbindlist(acc, fill = TRUE)
fwrite(oy, file.path(out_dir, "origin_year_pm25_reallocation.csv"))

# ---------------------------
# 5. Who avoids PM2.5? Flow-weighted reduction by subgroup
#    (avoidance = -reallocation; positive = exposure fell after moving)
# ---------------------------
avoid <- rbindlist(lapply(paste0("o_", sub_cols), function(g){
  s <- oy[!is.na(get(g)), .(mean_change = wmean(re_pm25, out_flow),
                            avoidance = -wmean(re_pm25, out_flow),
                            orig_pm25 = wmean(orig_pm25, out_flow),
                            flow = sum(out_flow, na.rm=TRUE)), by = g]
  setnames(s, g, "subgroup"); s[, dimension := sub("^o_", "", g)][]
}), fill = TRUE)
fwrite(avoid, file.path(out_dir, "avoidance_by_subgroup.csv"))
cat("\n=== PM2.5 reduction achieved by moving (flow-weighted), by subgroup ===\n")
print(avoid[order(dimension, -avoidance)])

# ---------------------------
# 6. Adjusted differential-avoidance models
# ---------------------------
run_diff <- function(dim_col, ref) {
  d <- oy[!is.na(get(dim_col)) & out_flow > 0 & !is.na(o_med_hh_income) & o_med_hh_income > 0]
  d[, grp := relevel(factor(get(dim_col)), ref = ref)]
  f <- as.formula("re_pm25 ~ grp + log(o_med_hh_income) + o_poverty_pct + o_bach_pct | o_state_county + year")
  fit <- tryCatch(feols(f, data = d, weights = ~out_flow, cluster = ~o_state_county), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ct <- as.data.table(coeftable(fit), keep.rownames = "term")[grepl("^grp", term)]
  ct[, dimension := dim_col][]
}
models <- rbindlist(list(
  if (!is.null(race)) run_diff("o_race_majority", "Majority white (NH)"),
  run_diff("o_inc_tercile", "High income"),
  run_diff("o_gent_group", "Eligible, not gentrifying")
), fill = TRUE)
fwrite(models, file.path(out_dir, "differential_avoidance_models.csv"))
cat("\n=== Differential PM2.5 reallocation vs reference (adjusted; + = reached dirtier air) ===\n"); print(models)

# ---------------------------
# 7. Sorting vs local cleanup
# ---------------------------
cleanup <- merge(poll[year == min(migrate_years), .(geoid_bg, pm0 = pm25)],
                 poll[year == max(migrate_years), .(geoid_bg, pm1 = pm25)], by = "geoid_bg")
cleanup[, own_cleanup := pm1 - pm0]
dim_use <- if (!is.null(race)) "o_race_majority" else "o_inc_tercile"
decomp <- merge(oy[, c("origin_bg","re_pm25","out_flow", dim_use), with = FALSE],
                cleanup[, .(origin_bg = geoid_bg, own_cleanup)], by = "origin_bg", all.x = TRUE)
dtab <- decomp[!is.na(get(dim_use)), .(sorting = wmean(re_pm25, out_flow),
               own_cleanup = wmean(own_cleanup, out_flow)), by = dim_use]
fwrite(dtab, file.path(out_dir, "sorting_vs_cleanup_pm25.csv"))
cat("\n=== PM2.5: exposure change via MOVING (sorting) vs origin's own CLEANUP, by subgroup ===\n"); print(dtab)

# ---------------------------
# 8. Exposure-reallocation signature by subgroup (NOT a health-impact assessment).
#    We do NOT convert to averted deaths: out_flow is summed move-events across
#    years (double-counting individuals), the exposure term is an adjusted group
#    contrast rather than an intervention-induced individual change, and there is
#    no exposure duration, age-specific baseline rate, lag, or uncertainty
#    propagation. We therefore report only the flow-weighted exposure reduction
#    and the exposure reallocated (person-microgram signature), as environmental-
#    health CONTEXT, not mortality. A proper HIA would require mover person-years,
#    duration, age-specific rates, validated CRFs, and a causal contrast.
# ---------------------------
exsig <- rbindlist(lapply(paste0("o_", sub_cols), function(g){
  s <- oy[!is.na(get(g)), .(delta = wmean(re_pm25, out_flow), move_events = sum(out_flow, na.rm=TRUE)), by = g]
  setnames(s, g, "subgroup")
  s[, exposure_reallocated := delta * move_events]   # person-microgram signature (negative = reduction); NOT deaths
  s[, dimension := sub("^o_", "", g)][]
}), fill = TRUE)
fwrite(exsig, file.path(out_dir, "exposure_signature_by_subgroup.csv"))
cat("\n=== Exposure-reallocation signature by subgroup (flow-weighted delta; NOT averted deaths) ===\n"); print(exsig[order(dimension, delta)])

cat("\nDONE. Outputs in ", out_dir, "\n", sep = "")
