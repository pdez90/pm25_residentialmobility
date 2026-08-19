# ============================================================
# Reproduce the paper pipeline.
#
#   Rscript run_main.R
#
# Requires only the public inputs listed in README.md; no credentials.
# Set MIGRATE_DATA_DIR to your data root first, e.g.
#   export MIGRATE_DATA_DIR=/path/to/data
#
# The two ancillary analyses that need external credentials or a manually
# downloaded file (the ACS microdata replication and the exploratory
# nonattainment analysis) are OFF here. Use run_full.R for those.
# ============================================================
source(file.path("scripts", "run_all_pm25.R"))
