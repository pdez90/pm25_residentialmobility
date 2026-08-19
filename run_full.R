# ============================================================
# Reproduce the paper pipeline PLUS the two ancillary analyses.
#
#   Rscript run_full.R
#
# Additional requirements beyond run_main.R:
#   * IPUMS USA extract for the individual-level ACS replication.
#     Build it with:  Rscript scripts/acs_ipums_pull.R   (needs IPUMS_API_KEY)
#   * EPA Green Book phistory.xls for the exploratory nonattainment
#     analysis. Point EPA_GREENBOOK_XLS at it.
#
# These stages fail loudly rather than silently skipping, so a missing
# input stops the run with an explicit message.
# ============================================================
Sys.setenv(MIGRATE_RUN_VALIDATION      = "TRUE")
Sys.setenv(MIGRATE_RUN_QUASIEXPERIMENT = "TRUE")
source(file.path("scripts", "run_all_pm25.R"))
