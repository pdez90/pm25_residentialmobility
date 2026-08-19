# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# MASTER DRIVER (reproducible): PM2.5 x MIGRATION analysis, end to end
# ------------------------------------------------------------
# Runs the whole pipeline reproducibly and self-verifies:
#   * fixed RNG seed + kind                         (reproducible randomness)
#   * full run LOG teed to a timestamped file       (reproducible record)
#   * R version, platform, and package versions      logged at the top
#   * every stage ASSERTS its expected outputs exist (embedded per-stage checks)
#   * a comprehensive DIAGNOSTICS stage at the end   (run_diagnostics.R)
#   * an output MANIFEST with row counts + md5        (compare runs bit-for-bit)
#   * any failure STOPS the run with a clear message (nothing silently skipped;
#     a stage you don't want is turned OFF via its switch)
#
# Usage:  Rscript run_main.R   (paper pipeline, no credentials required)
#         Rscript run_full.R   (adds ACS microdata replication + exploratory NAAQS)
# Outputs: $MIGRATE_DATA_DIR/pm25_avoidance_lag3_olag1/
#   logs/run_<timestamp>.log        full console+message log
#   logs/manifest_<timestamp>.csv   file, bytes, nrow, md5 for every output
#   logs/sessioninfo_<timestamp>.txt sessionInfo()
#
# STAGES: 0a PM2.5 | 0b components | 0c multi-pollutant | 0c* padding-fix rebuild |
#   0d uncertainty | 1 core avoidance | 2 figures/specs | 3 sorting-vs-siting |
#   4 cross-pollutant | 5 diagnostics-maps | 6 sankey | 7 affordability robustness |
#   8 opportunity/segregation | 9 ACS household validation | 10 NAAQS quasi-experiment |
#   11 DIAGNOSTICS (hard-fail verification)
# Stages 9-10 need data you supply; if their switch is ON and the input is
# missing the run STOPS (they auto-build from raw inputs when possible).
# ============================================================

scripts_dir <- file.path(repo_root, "scripts")
base_dir    <- migrate_data_dir()
poll_dir    <- file.path(base_dir, "pollution_bg")
out_dir     <- file.path(base_dir, "pm25_avoidance_lag3_olag1")
log_dir     <- file.path(out_dir, "logs")

# ---- switches ----
RUN_EXTRACTION      <- TRUE     # 0a PM2.5 panel
RUN_COMPONENTS      <- TRUE     # 0b speciated SO4/NO3/BC/OM
RUN_MULTIPOLLUTANT  <- TRUE     # 0c PM2.5+NO2+wildfire+ozone caches/panel
RUN_REBUILD_PANEL   <- TRUE     # 0c* re-merge caches with geoid padding fix (fast; always)
RUN_UNCERTAINTY     <- TRUE     # 0d PM2.5 uncertainty
RUN_DIAGNOSTICS     <- TRUE     # 5-6 maps + sankey, and 11 verification
RUN_ROBUSTNESS      <- TRUE     # 7 affordability / selection / ME / mean reversion
RUN_OPPORTUNITY     <- TRUE     # 8 opportunity set + segregation + affordability+
RUN_VALIDATION      <- as.logical(Sys.getenv("MIGRATE_RUN_VALIDATION", "FALSE"))   # 9 ACS microdata replication (needs IPUMS extract)
RUN_QUASIEXPERIMENT <- as.logical(Sys.getenv("MIGRATE_RUN_QUASIEXPERIMENT", "FALSE")) # 10 NAAQS, exploratory (needs EPA Green Book)
FORCE_REEXTRACT     <- FALSE    # re-run heavy extraction even if a panel exists
SEED                <- 20260709 # fixed seed for reproducible randomness

# ============================================================
# MAIN (wrapped so logging/sinks always close, even on error)
# ============================================================
main <- function() {
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  stamp    <- format(Sys.time(), "%Y%m%d_%H%M%S")
  logfile  <- file.path(log_dir, paste0("run_", stamp, ".log"))
  manifest <- file.path(log_dir, paste0("manifest_", stamp, ".csv"))
  sifile   <- file.path(log_dir, paste0("sessioninfo_", stamp, ".txt"))
  logcon   <- file(logfile, open = "wt")
  sink(logcon, split = TRUE)             # stdout -> console + log
  sink(logcon, type = "message")         # messages/warnings/errors -> log
  options(warn = 1)                      # surface warnings as they occur
  t_start <- Sys.time()
  on.exit({ sink(type = "message"); sink(); close(logcon) }, add = TRUE)

  # ---- reproducibility header ----
  set.seed(SEED); RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
  cat("============================================================\n")
  cat("REPRODUCIBLE RUN  ", stamp, "\n", sep = "")
  cat("============================================================\n")
  cat("R:        ", R.version.string, "\n", sep = "")
  cat("platform: ", R.version$platform, "\n", sep = "")
  cat("seed:     ", SEED, "  RNGkind: ", paste(RNGkind(), collapse=","), "\n", sep = "")
  cat("log:      ", logfile, "\n\n", sep = "")
  writeLines(capture.output(sessionInfo()), sifile)

  sfile       <- function(f) file.path(scripts_dir, f)
  rscript_bin <- file.path(R.home("bin"), "Rscript")
  exists_ne   <- function(f) file.exists(f) && file.info(f)$size > 0

  # ---- one-time dependency check + version log ----
  need <- c("data.table","stringr","fixest","splines","ggplot2","ggrepel","sf","terra",
            "exactextractr","tigris","tidycensus","readxl","ggalluvial","foreign")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) { cat("Installing: ", paste(miss, collapse=", "), "\n")
                      install.packages(miss, repos = "https://cloud.r-project.org") }
  cat("Package versions:\n")
  for (p in need) cat(sprintf("  %-14s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e) "MISSING")))
  cat("\n")

  # ---- stage runners (both tee into the log via the active sink) ----
  run_sub <- function(f) {
    cat("\n========== SUBPROCESS: ", f, " ==========\n", sep=""); t0 <- Sys.time()
    res <- system2(rscript_bin, shQuote(sfile(f)), stdout = TRUE, stderr = TRUE)
    cat(res, sep = "\n"); cat("\n")
    st <- attr(res, "status")
    if (!is.null(st) && st != 0) stop("STAGE FAILED (exit ", st, "): ", f)
    cat("========== DONE: ", f, "  (", round(difftime(Sys.time(), t0, units="mins"),1), " min) ==========\n", sep="")
  }
  run_ses <- function(f) {
    cat("\n========== STAGE: ", f, " ==========\n", sep=""); t0 <- Sys.time()
    source(sfile(f), local = FALSE, echo = FALSE)
    cat("========== DONE: ", f, "  (", round(difftime(Sys.time(), t0, units="mins"),1), " min) ==========\n", sep="")
  }
  # embedded per-stage assertion: required outputs must exist and be non-empty
  assert <- function(stage, files) {
    missing <- files[!vapply(files, exists_ne, logical(1))]
    if (length(missing)) stop("STAGE '", stage, "' did not produce: ", paste(basename(missing), collapse=", "))
    cat("  [check] ", stage, ": ", length(files), " expected output(s) present.\n", sep="")
  }
  P <- function(...) file.path(poll_dir, ...); O <- function(...) file.path(out_dir, ...)

  # ============================================================
  # STAGE 0: extraction (heavy; separate processes; cached)
  # ============================================================
  if (RUN_EXTRACTION && (FORCE_REEXTRACT || !exists_ne(P("bg_pm25_panel.csv")))) run_sub("extract_pm25_bg.R")
  assert("0a PM2.5 panel", P("bg_pm25_panel.csv"))

  if (RUN_COMPONENTS && (FORCE_REEXTRACT || !exists_ne(P("bg_pm25_components_panel.csv")))) run_sub("extract_pm25_components_bg.R")
  if (RUN_COMPONENTS) assert("0b components", P("bg_pm25_components_panel.csv"))

  if (RUN_MULTIPOLLUTANT && (FORCE_REEXTRACT || !exists_ne(P("bg_pollution_panel.csv")))) run_sub("extract_pollution_multi_bg.R")
  if (RUN_REBUILD_PANEL && dir.exists(P("cache_multi"))) run_sub("rebuild_pollution_panel.R")   # 0c*: padding fix
  if (RUN_MULTIPOLLUTANT) assert("0c multi-pollutant panel", P("bg_pollution_panel.csv"))

  if (RUN_UNCERTAINTY && (FORCE_REEXTRACT || !exists_ne(P("bg_pm25_uncertainty_panel.csv")))) run_sub("extract_pm25_uncertainty_bg.R")
  if (RUN_UNCERTAINTY) assert("0d uncertainty", P("bg_pm25_uncertainty_panel.csv"))

  # ============================================================
  # STAGES 1-3: core analysis chain (one session; reads bg_pm25_panel.csv)
  # ============================================================
  run_ses("idea2_pm25_avoidance.R");    assert("1 core avoidance", O("avoidance_by_subgroup.csv"))
  run_ses("idea2_pm25_avoidance_v3.R"); assert("2 figures/specs",  c(O("disparity_trends.csv"), O("conditional_models_all_fe.csv")))
  run_ses("pm25_sorting_vs_siting.R");  assert("3 sorting/siting", c(O("siting_pollution_dynamics.csv"), O("sorting_destination_pull.csv")))

  # ============================================================
  # STAGE 4: cross-pollutant
  # ============================================================
  if (RUN_MULTIPOLLUTANT) { run_ses("idea2_pollution_crosscut.R"); assert("4 cross-pollutant", O("cross_pollutant_estimates.csv")) }
# 4b: descriptive disparity trends + shift-share decomposition for NO2 and ozone (SI)
if (RUN_MULTIPOLLUTANT) { run_ses("multipollutant_disparity.R")
  assert("4b multi-pollutant disparity", c(O("multipollutant_disparity_trends.csv"), O("multipollutant_decomposition.csv"))) }

  # ============================================================
  # STAGE 5-6: diagnostics maps + Sankey (SI)
  # ============================================================
  if (RUN_DIAGNOSTICS) {
    run_ses("pollution_diagnostics_maps.R"); assert("5 diagnostics-maps", O("diagnostics_coverage.csv"))
    run_ses("pollution_sankey.R");           assert("6 sankey", O("sankey_flows_pm25.csv"))
  }

  # ============================================================
  # STAGE 7: methodological hardening
  # ============================================================
  if (RUN_ROBUSTNESS) { run_ses("affordability_robustness.R")
    assert("7 robustness", c(O("affordability_decomposition.csv"), O("gelbach_decomposition.csv"),
                             O("selection_robustness.csv"),
                             O("measurement_error_robustness.csv"), O("pollutant_correlations.csv"),
                             O("origin_year_affordability.csv"))) }

  # ============================================================
  # STAGE 8: opportunity + segregation + affordability+
  # ============================================================
  if (RUN_OPPORTUNITY) { run_ses("opportunity_heterogeneity.R")
    assert("8 opportunity", c(O("opportunity_by_race.csv"), O("segregation_interaction.csv"), O("affordability_plus.csv"))) }

  # ============================================================
  # STAGE 8b: multi-pollutant robustness (affordability, selection,
  #           opportunity, segregation for ALL eight pollutants)
  # ============================================================
  if (RUN_OPPORTUNITY) { run_ses("multipollutant_robustness.R")
    assert("8b multipollutant robustness", c(O("multipollutant_affordability.csv"), O("multipollutant_selection.csv"),
                                             O("multipollutant_opportunity.csv"), O("multipollutant_segregation.csv"))) }

  # ============================================================
  # STAGE 9: ACS household validation (ON => input required, else STOP)
  # ============================================================
  if (RUN_VALIDATION) {
    if (!exists_ne(file.path(base_dir,"ipums","usa_acs_migration.csv")))
      stop("STAGE 9 ON but ipums/usa_acs_migration.csv missing. Run: Rscript scripts/acs_ipums_pull.R (needs IPUMS_API_KEY), or set RUN_VALIDATION <- FALSE.")
    if (!exists_ne(file.path(base_dir,"ipums","puma_migpuma_xwalk.csv"))) run_sub("build_puma_migpuma_xwalk.R")
    run_ses("acs_micro_validation.R"); assert("9 ACS validation", O("acs_household_gap.csv"))
  }

  # ============================================================
  # STAGE 10: NAAQS quasi-experiment (ON => input required, else STOP)
  # ============================================================
  if (RUN_QUASIEXPERIMENT) {
    gb <- file.path(base_dir,"naaqs","greenbook_nonattainment.csv")
    if (!exists_ne(gb) && any(file.exists(c(greenbook_xls(), file.path(base_dir,"naaqs","phistory.xls")))))
      run_sub("build_greenbook.R")
    if (!exists_ne(gb)) stop("STAGE 10 ON but ", gb, " missing and no phistory.xls to build it. Set EPA_GREENBOOK_XLS to your Green Book phistory.xls, or leave the optional stages off (run_main.R).")
    run_ses("naaqs_nonattainment.R"); assert("10 NAAQS", O("naaqs_capture_origin.csv"))
  }

  # ============================================================
  # STAGE 11: comprehensive DIAGNOSTICS (hard-fails on any problem)
  # ============================================================
  if (RUN_DIAGNOSTICS) run_ses("run_diagnostics.R")

  # ============================================================
  # OUTPUT MANIFEST (reproducibility record)
  # ============================================================
  cat("\n== writing output manifest ==\n")
  csvs <- list.files(out_dir, pattern="\\.csv$", full.names=TRUE)
  man <- data.frame(
    file  = basename(csvs),
    bytes = file.info(csvs)$size,
    nrow  = vapply(csvs, function(f) tryCatch(nrow(data.table::fread(f, showProgress=FALSE)), error=function(e) NA_integer_), integer(1)),
    md5   = unname(tools::md5sum(csvs)),
    stringsAsFactors = FALSE)
  data.table::fwrite(man, manifest)
  cat("  manifest: ", manifest, "  (", nrow(man), " output tables)\n", sep="")

  # ============================================================
  # SUMMARY
  # ============================================================
  cat("\n============================================================\n")
  cat("ALL STAGES COMPLETE + DIAGNOSTICS PASSED\n")
  cat("  total time: ", round(difftime(Sys.time(), t_start, units="mins"),1), " min\n", sep="")
  cat("  outputs:    ", out_dir, "\n", sep="")
  cat("  log:        ", logfile, "\n", sep="")
  cat("  manifest:   ", manifest, "\n", sep="")
  cat("  sessioninfo:", sifile, "\n", sep="")
  cat("============================================================\n")
}

# run; on any error, report and exit non-zero (so CI / shell sees the failure)
tryCatch(main(), error = function(e) {
  message("\n*** RUN FAILED: ", conditionMessage(e))
  quit(save = "no", status = 1)
})
