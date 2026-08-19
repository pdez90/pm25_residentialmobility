# ============================================================
# PULL the ACS migration extract from the IPUMS API (ipumsr)
# ------------------------------------------------------------
# Defines, submits, waits for, and downloads the exact IPUMS USA extract
# that acs_micro_validation.R consumes -- no manual extract-builder needed.
#
# ONE-TIME SETUP:
#   1. Register at https://usa.ipums.org and create an API key at
#      https://account.ipums.org/api_keys
#   2. Put the key in ~/.Renviron (NOT in this file / not in chat):
#        IPUMS_API_KEY=xxxxxxxxxxxxxxxxxxxx
#      then restart R (or run readRenviron("~/.Renviron")).
#
# Output: MIGRATE_dir/ipums/usa_acs_migration.csv  (feeds the validation)
# NOTE: still obtain the small PUMA->MIGPUMA correspondence separately
#       (ipums/puma_migpuma_xwalk.csv) -- the API does not provide it.
# ============================================================
if (!requireNamespace("ipumsr", quietly=TRUE)) install.packages("ipumsr", repos="https://cloud.r-project.org")
suppressPackageStartupMessages({ library(ipumsr); library(data.table) })
base_dir <- "/Users/priyanka/MIGRATE_dir"
ip_dir <- file.path(base_dir,"ipums"); dir.create(ip_dir, showWarnings=FALSE, recursive=TRUE)

key <- Sys.getenv("IPUMS_API_KEY")
stopifnot("Set IPUMS_API_KEY in ~/.Renviron (see header)" = nzchar(key))
set_ipums_api_key(key)

# ACS 1-year samples 2012-2019 (sample id convention: us<year>a = ACS 1-yr)
samples <- paste0("us", 2012:2019, "a")
vars <- c("YEAR","SAMPLE","SERIAL","PERNUM","PERWT",
          "STATEFIP","PUMA","MIGRATE1","MIGPLAC1","MIGPUMA1",
          "RACE","RACED","HISPAN","HISPAND","AGE","INCTOT","EDUCD",
          "OWNERSHP")   # OWNERSHP enables the renter-vs-owner subsample (Table S19)

cat("Defining extract (", length(samples), " samples, ", length(vars), " variables)...\n", sep="")
ext <- define_extract_micro(
  collection  = "usa",
  description = "ACS 2012-2019 1yr: migration (MIGPUMA1) x race for pollution-avoidance validation",
  samples     = samples,
  variables   = vars,
  data_format = "csv")

cat("Submitting to IPUMS API...\n")
sub <- submit_extract(ext)
cat("Extract number:", sub$number, " -- waiting (this can take several minutes)...\n")
sub <- wait_for_extract(sub)               # polls until ready
stopifnot("extract did not complete" = sub$status == "completed")

paths <- download_extract(sub, download_dir = ip_dir, overwrite = TRUE)
cat("Downloaded:\n"); print(paths)

# read (ipumsr returns the labelled microdata) and write the plain CSV the
# validation script expects (numeric codes, one row per person-wave)
ddi <- paths[grepl("\\.xml$", paths)][1]
d <- read_ipums_micro(ddi, verbose = FALSE)
setDT(d); setnames(d, toupper(names(d)))
keep <- intersect(toupper(vars), names(d))
d <- d[, ..keep]
for (c in keep) d[[c]] <- as.numeric(zap_labels(d[[c]]))   # strip value labels -> plain codes
fwrite(d, file.path(ip_dir,"usa_acs_migration.csv"))
cat("\nDONE. Wrote ", file.path(ip_dir,"usa_acs_migration.csv"),
    " (", nrow(d), " person-wave records). Now run acs_micro_validation.R.\n", sep="")
