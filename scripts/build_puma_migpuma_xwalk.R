# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# BUILD the PUMA -> MIGPUMA crosswalk (2010 vintage)
# ------------------------------------------------------------
# Downloads the OFFICIAL equivalence file (2010 MIGPUMA1 = PWPUMA00),
# hosted by IPUMS, and writes the small crosswalk acs_micro_validation.R
# needs:  ipums/puma_migpuma_xwalk.csv  (columns: statefip, puma, migpuma)
#
# No API key needed. Runs on your machine (network access required).
# Source: https://usa.ipums.org/usa/resources/volii/puma_migpuma1_pwpuma00_2010.xls
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr) })
if (!requireNamespace("readxl", quietly=TRUE)) install.packages("readxl", repos="https://cloud.r-project.org")
library(readxl)
base_dir <- migrate_data_dir()
ip_dir <- file.path(base_dir,"ipums"); dir.create(ip_dir, showWarnings=FALSE, recursive=TRUE)
dest <- file.path(ip_dir,"puma_migpuma_xwalk.csv")
raw  <- file.path(ip_dir,"puma_migpuma1_pwpuma00_2010.xls")

urls <- c(
  "https://usa.ipums.org/usa/resources/volii/puma_migpuma1_pwpuma00_2010.xls",
  "https://assets.ipums.org/_files/volii/puma_migpuma1_pwpuma00_2010.xls"
)
ok <- file.exists(raw) && file.info(raw)$size > 5000
if (!ok) for (u in urls) {
  message("downloading ", u)
  try({ download.file(u, raw, quiet=TRUE, mode="wb")
        if (file.exists(raw) && file.info(raw)$size > 5000) { ok <- TRUE; break } }, silent=TRUE)
}
if (!ok) stop("Could not download the crosswalk. Download it manually from\n",
              "  https://usa.ipums.org/usa/resources/volii/puma_migpuma1_pwpuma00_2010.xls\n",
              "and save it at ", raw, ", then re-run this script.")

# Layout (2010-based file): row1 title, row2 blank, row3 header, then data.
# Header columns:
#   1 "State of Residence (ST)"                         -> statefip
#   2 "PUMA"                                            -> puma
#   3 "Place of Work State (PWSTATE2) or Migration ..." -> (ignored; = state)
#   4 "PWPUMA00 or MIGPUMA1"                            -> migpuma
# Auto-locate the header row, then take columns by name with a positional
# fallback (readxl strips leading zeros, so we re-pad from the integer).
tab <- as.data.table(read_excel(raw, col_names=FALSE, .name_repair="minimal"))
# the TITLE row also contains "PUMA ... MIGPUMA", so require an exact "PUMA"
# cell (only the real header row has a standalone PUMA column)
hdr <- which(apply(tab, 1, function(r){
  r <- toupper(trimws(as.character(r))); any(r=="PUMA", na.rm=TRUE) & any(grepl("PWPUMA|MIGPUMA", r)) }))[1]
stopifnot("could not find the header row in the xls" = !is.na(hdr))
nm  <- toupper(trimws(as.character(unlist(tab[hdr]))))
dat <- tab[(hdr+1):.N]
pick <- function(pat, pos){ i <- which(grepl(pat, nm)); if (length(i)) i[1] else pos }
i_st <- pick("STATE OF RESIDENCE|STATEFIP|^ST$", 1)
i_pu <- which(nm=="PUMA"); i_pu <- if (length(i_pu)) i_pu[1] else 2
i_mg <- pick("PWPUMA00 OR MIGPUMA1|PWPUMA|MIGPUMA1", 4)
clean_int <- function(x) suppressWarnings(as.integer(gsub("[^0-9]","", as.character(x))))
xw <- data.table(statefip=str_pad(clean_int(dat[[i_st]]),2,pad="0"),
                 puma    =str_pad(clean_int(dat[[i_pu]]),5,pad="0"),
                 migpuma =str_pad(clean_int(dat[[i_mg]]),5,pad="0"))
xw <- xw[!is.na(statefip) & !is.na(puma) & !is.na(migpuma) & statefip!="0NA" & grepl("^[0-9]{2}$",statefip)]
xw <- unique(xw, by=c("statefip","puma"))
stopifnot("crosswalk came out empty; inspect the xls layout" = nrow(xw) > 1000)
fwrite(xw, dest)
cat("DONE. Wrote ", dest, " (", nrow(xw), " PUMA->MIGPUMA rows across ",
    length(unique(xw$statefip)), " states).\n", sep="")
