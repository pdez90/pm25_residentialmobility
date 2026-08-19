# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# BUILD the EPA Green Book nonattainment file for the NAAQS stage
# ------------------------------------------------------------
# Reshapes the Green Book PHISTORY export (county x NAAQS x year; cells
# W/P/blank = whole/partial/not nonattainment, 1992-present) into
#     naaqs/greenbook_nonattainment.csv   (county, pollutant, design_year)
# where design_year = the first year the county was nonattainment for the
# chosen PM2.5 standard.
#
# PHISTORY columns (confirmed): pollutant, revoked_naaqs, state_name,
#   st_abbr, fips_state (2), fips_cnty (3), countyname, pw_1992..pw_2026, exportdt
#
# STANDARD selects the shock:
#   "2012" (default) -> 20 counties designated 2015  = IN-WINDOW exogenous shock
#                       (populous: CA San Joaquin Valley + South Coast, etc.)
#   "2006"           -> ~121 counties designated ~2009 (persistent, mostly pre-window)
#   "all"            -> any PM2.5 nonattainment (broadest treatment group)
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr) })
STANDARD <- "2012"
base_dir <- migrate_data_dir()
nd <- file.path(base_dir,"naaqs"); dir.create(nd, showWarnings=FALSE, recursive=TRUE)
dest <- file.path(nd,"greenbook_nonattainment.csv")

# ---- locate phistory.xls (your EPA_greenbook folder first), else download ----
cands <- c(greenbook_xls(), file.path(nd,"phistory.xls"))
src <- cands[file.exists(cands)][1]
if (is.na(src)) {
  src <- file.path(nd,"phistory.xls")
  ok <- FALSE
  try({ download.file("https://www3.epa.gov/airquality/greenbook/downld/phistory.xls", src, quiet=TRUE, mode="wb")
        ok <- file.exists(src) && file.info(src)$size>5000 }, silent=TRUE)
  if (!ok) stop("phistory.xls not found at $EPA_GREENBOOK_XLS or ", nd,
                " and could not be downloaded. Place the file and re-run.")
}
if (!requireNamespace("readxl", quietly=TRUE)) install.packages("readxl", repos="https://cloud.r-project.org")
d <- as.data.table(readxl::read_excel(src, col_names=TRUE, .name_repair="minimal"))
setnames(d, tolower(names(d)))
req <- c("pollutant","fips_state","fips_cnty")
stopifnot("PHISTORY missing expected columns (pollutant, fips_state, fips_cnty)" = all(req %in% names(d)))

# ---- county FIPS + PM2.5 filter ----
d[, county := paste0(str_pad(gsub("[^0-9]","",fips_state),2,pad="0"),
                     str_pad(gsub("[^0-9]","",fips_cnty),3,pad="0"))]
d <- d[grepl("PM-?2\\.?5", pollutant, ignore.case=TRUE)]
if (STANDARD != "all") d <- d[grepl(STANDARD, pollutant)]
stopifnot("no PM2.5 rows for the chosen STANDARD" = nrow(d)>0)

# ---- first nonattainment year per county ----
ycols <- grep("^pw_[0-9]{4}$", names(d), value=TRUE)
stopifnot("no pw_YYYY year columns found" = length(ycols) >= 5)
long <- melt(d, id.vars="county", measure.vars=ycols, variable.name="yc", value.name="status")
long <- long[toupper(trimws(status)) %in% c("W","P")]
long[, year := as.integer(str_extract(as.character(yc), "[0-9]{4}"))]
naa <- long[, .(design_year = min(year, na.rm=TRUE)), by=county]
naa[, pollutant := "PM2.5"]
naa <- naa[grepl("^[0-9]{5}$", county), .(county, pollutant, design_year)]
setorder(naa, design_year, county)
fwrite(naa, dest)
cat(sprintf("DONE. Wrote %s\n  %d PM2.5 (%s standard) nonattainment counties; design_year %d-%d (median %d)\n",
            dest, nrow(naa), STANDARD, min(naa$design_year), max(naa$design_year), as.integer(median(naa$design_year))))
if (uniqueN(naa$county) < 30)
  cat("  NOTE: few treated counties -> the NAAQS event study has <30 clusters;\n",
      "  cluster-robust SEs are approximate. Consider STANDARD <- \"all\" for the\n",
      "  capture tests, or interpret the first stage as suggestive.\n", sep="")
