# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# DIAGNOSTICS / VERIFICATION  (stage 11)
# ------------------------------------------------------------
# Verifies the pipeline actually worked and HARD-FAILS on any problem.
# It collects EVERY problem it finds, prints a full report, then stops if
# the list is non-empty (so you see all issues at once, not just the first).
#
# Checks:
#   A  Panels exist, are non-empty, and have NO duplicate BG-year rows
#   B  PM2.5 coverage ~100% every migration year, and no whole state
#      dropped pre-2015  (catches the geoid zero-padding bug)
#   C  The two PM2.5 sources (bg_pm25_panel vs bg_pollution_panel) AGREE
#   D  Flow-join and race-join coverage across ALL migration years 2011-2019
#   E  Pollutant value ranges are physically plausible
#   F  Every expected core output file exists and is non-empty
#   (optional outputs -- ACS / NAAQS -- are checked only if their stage ran)
#
# Run standalone:  Rscript scripts/run_diagnostics.R
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr) })
base_dir <- migrate_data_dir(); mig_dir <- base_dir
poll_dir <- file.path(base_dir,"pollution_bg"); out_dir <- file.path(base_dir,"pm25_avoidance_lag3_olag1")
migrate_years <- 2011:2019
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")
PROB <- character(0)
bad  <- function(msg) PROB[[length(PROB)+1]] <<- msg
ok   <- function(msg) cat("  [ok] ", msg, "\n", sep="")
sect <- function(s) cat("\n== ", s, " ==\n", sep="")
rd   <- function(f, ...) fread(f, colClasses=c(geoid_bg="character"), ...)

# ---- shared inputs ----
popf <- file.path(base_dir,"acs_gentrification","blockgroup_gentrification_2009_2019_lag3_with_HMDA.csv")
pop  <- rd(popf, select=c("year","geoid_bg","total_pop")); pop[, geoid_bg:=pad_bg(geoid_bg)]; pop[is.na(total_pop), total_pop:=0]

# ============================================================
sect("A/B  PM2.5 core panel (bg_pm25_panel.csv)")
f_core <- file.path(poll_dir,"bg_pm25_panel.csv")
if (!file.exists(f_core) || file.info(f_core)$size==0) bad("bg_pm25_panel.csv missing/empty") else {
  pm <- rd(f_core); pm[, geoid_bg:=pad_bg(geoid_bg)]
  ndup <- nrow(pm) - uniqueN(pm[, .(geoid_bg, year)])
  if (ndup>0) bad(sprintf("bg_pm25_panel.csv has %d duplicate BG-year rows", ndup)) else ok("no duplicate BG-year rows")
  m <- merge(pm[year %in% migrate_years, .(geoid_bg,year,pm25)], pop, by=c("geoid_bg","year"))
  cov <- m[, .(bg=round(100*mean(!is.na(pm25)),1), pop=round(100*sum(total_pop*!is.na(pm25))/sum(total_pop),1)), by=year][order(year)]
  print(cov)
  if (any(cov$pop < 99)) bad(sprintf("bg_pm25_panel.csv pop-weighted PM2.5 coverage <99%% in year(s): %s",
                                      paste(cov[pop<99]$year, collapse=", "))) else ok("PM2.5 pop-weighted coverage >=99% all years")
}

# ============================================================
sect("B  Multi-pollutant panel (bg_pollution_panel.csv) + zero-padding-bug guard")
f_multi <- file.path(poll_dir,"bg_pollution_panel.csv")
mp <- NULL
if (!file.exists(f_multi) || file.info(f_multi)$size==0) bad("bg_pollution_panel.csv missing/empty") else {
  mp <- rd(f_multi); mp[, geoid_bg:=pad_bg(geoid_bg)]
  ndup <- nrow(mp) - uniqueN(mp[, .(geoid_bg, year)])
  if (ndup>0) bad(sprintf("bg_pollution_panel.csv has %d duplicate BG-year rows (padding bug not fixed)", ndup)) else ok("no duplicate BG-year rows")
  pols <- intersect(c("pm25","no2","wildfire_pm25","ozone"), names(mp))
  cov <- mp[year %in% migrate_years, c(list(year=year), lapply(.SD, function(x) round(100*mean(!is.na(x)),1))), .SDcols=pols][
            , lapply(.SD, mean), by=year][order(year)]
  print(cov)
  if (any(cov$pm25 < 99)) bad(sprintf("bg_pollution_panel.csv PM2.5 coverage <99%% in year(s): %s",
                                       paste(cov[pm25<99]$year, collapse=", "))) else ok("PM2.5 coverage >=99% all years")
  # whole-state drop check (the bug signature) for a low year
  st <- mp[year==2012, .(pop_cov = 100*sum((geoid_bg %in% pop[year==2012 & total_pop>0]$geoid_bg) & !is.na(pm25))/.N,
                         cov=100*mean(!is.na(pm25))), by=.(state=substr(geoid_bg,1,2))]
  zero <- st[cov < 5]$state
  if (length(zero)) bad(sprintf("bg_pollution_panel.csv: whole states with ~0%% PM2.5 in 2012: %s (zero-padding bug)", paste(zero, collapse=", "))) else ok("no whole state dropped in 2012")
}

# ============================================================
sect("C  The two PM2.5 sources agree")
if (exists("pm") && !is.null(mp)) {
  cmp <- merge(pm[, .(geoid_bg,year,pm_core=pm25)], mp[, .(geoid_bg,year,pm_multi=pm25)], by=c("geoid_bg","year"))
  cmp <- cmp[!is.na(pm_core) & !is.na(pm_multi)]
  if (nrow(cmp)==0) bad("no overlap between the two PM2.5 panels to compare") else {
    r <- cor(cmp$pm_core, cmp$pm_multi); mad <- max(abs(cmp$pm_core - cmp$pm_multi))
    cat(sprintf("  overlap=%d  correlation=%.4f  max|diff|=%.3f ug/m3\n", nrow(cmp), r, mad))
    if (r < 0.999 || mad > 0.5) bad(sprintf("PM2.5 panels disagree (cor=%.4f, max|diff|=%.3f)", r, mad)) else ok("PM2.5 sources agree")
  }
}

# ============================================================
sect("D  Flow-join and race-join coverage, ALL years")
rcs <- list.files(poll_dir, "^race_bg_.*csv$", full.names=TRUE)
rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}")); r<-rd(f); r[,year:=y]; r }), fill=TRUE)
rc[, geoid_bg:=pad_bg(geoid_bg)]
if (!is.null(mp)) {
  jr <- data.table()
  for (yr in migrate_years) {
    fx <- file.path(mig_dir, paste0(yr,".csv")); if (!file.exists(fx)) { bad(paste0("MIGRATE flow file missing: ",yr,".csv")); next }
    x <- fread(fx, select=c("Origin","Destination","Flow"))
    x <- x[, .(o=pad_bg(Origin), d=pad_bg(Destination), flow=as.numeric(Flow))][flow>0]
    keys <- mp[year==yr & !is.na(pm25), geoid_bg]
    oc <- 100*sum(x$flow*(x$o %in% keys))/sum(x$flow)
    dc <- 100*sum(x$flow*(x$d %in% keys))/sum(x$flow)
    rkeys <- rc[year==yr, geoid_bg]; racec <- 100*mean(keys %in% rkeys)
    jr <- rbind(jr, data.table(year=yr, origin_flow_cov=round(oc,1), dest_flow_cov=round(dc,1), race_join=round(racec,1)))
  }
  print(jr)
  if (any(jr$origin_flow_cov < 99 | jr$dest_flow_cov < 99))
    bad(sprintf("flow-join <99%% (flow-weighted) in year(s): %s", paste(jr[origin_flow_cov<99 | dest_flow_cov<99]$year, collapse=", "))) else ok("flow-join >=99% all years")
  if (any(jr$race_join < 95))
    bad(sprintf("race-join <95%% in year(s): %s", paste(jr[race_join<95]$year, collapse=", "))) else ok("race-join >=95% all years")
}

# ============================================================
sect("E  Pollutant value ranges")
if (!is.null(mp)) {
  rng <- list(pm25=c(1,40), no2=c(0,60), ozone=c(15,85), wildfire_pm25=c(0,60))
  for (p in names(rng)) if (p %in% names(mp)) {
    v <- mp[[p]][!is.na(mp[[p]])]; lo<-min(v); hi<-max(v); md<-median(v)
    cat(sprintf("  %-14s median %.2f  [%.2f, %.2f]\n", p, md, lo, hi))
    if (lo < rng[[p]][1]-1 || hi > rng[[p]][2]) bad(sprintf("%s out of plausible range [%.0f,%.0f]: observed [%.2f,%.2f]", p, rng[[p]][1], rng[[p]][2], lo, hi))
  }
}

# ============================================================
sect("F  Expected output files")
core_out <- c("disparity_trends.csv","exposure_decomposition.csv","gap_decomposition.csv",
  "conditional_models_all_fe.csv","avoidance_signature_by_subgroup.csv","siting_pollution_dynamics.csv",
  "siting_robustness.csv","sorting_destination_pull.csv","sorting_by_pm25_quintile.csv",
  "cross_pollutant_estimates.csv","affordability_decomposition.csv","gelbach_decomposition.csv","selection_robustness.csv",
  "measurement_error_robustness.csv","pollutant_correlations.csv","opportunity_by_race.csv",
  "segregation_interaction.csv","gap_by_segregation_quartile.csv","affordability_plus.csv")
for (f in core_out) { p <- file.path(out_dir,f)
  if (!file.exists(p) || file.info(p)$size==0) bad(paste0("missing/empty output: ", f)) }
miss <- sum(sapply(core_out, function(f) !file.exists(file.path(out_dir,f)) || file.info(file.path(out_dir,f))$size==0))
if (miss==0) ok(sprintf("all %d core outputs present and non-empty", length(core_out)))
# Gelbach decomposition: additivity identity must hold (b_base-b_full = sum_k delta_k*gamma_k)
gp <- file.path(out_dir,"gelbach_decomposition.csv")
if (file.exists(gp) && file.info(gp)$size>0) {
  g <- fread(gp)
  if (!"identity_gap" %in% names(g)) bad("gelbach_decomposition.csv missing identity_gap column") else {
    mg <- max(abs(g$identity_gap), na.rm=TRUE)
    if (!is.finite(mg) || mg > 1e-5) bad(sprintf("Gelbach additivity identity violated: max|identity_gap|=%.3e (>1e-5)", mg))
    else ok(sprintf("Gelbach decomposition additivity holds (max|identity_gap|=%.1e)", mg))
  }
}
# optional (only if that stage ran)
for (f in c("acs_household_gap.csv","naaqs_capture_origin.csv")) { p<-file.path(out_dir,f)
  if (file.exists(p)) ok(paste0("optional output present: ", f)) }

# ============================================================
cat("\n=================== DIAGNOSTICS SUMMARY ===================\n")
if (length(PROB)==0) {
  cat("ALL DIAGNOSTICS PASSED.\n")
} else {
  cat(sprintf("FOUND %d PROBLEM(S):\n", length(PROB)))
  for (i in seq_along(PROB)) cat(sprintf("  %d. %s\n", i, PROB[[i]]))
  stop("Diagnostics failed -- see the numbered problems above.")
}
