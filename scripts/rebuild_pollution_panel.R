# ============================================================
# FAST REBUILD of bg_pollution_panel.csv from cached extractions
# ------------------------------------------------------------
# Fixes the zero-padding merge bug (FIPS 01-09, incl. California, were
# duplicated / dropped for pre-2015 PM2.5) WITHOUT re-extracting any
# raster. Reads the per-pollutant cache CSVs, forces a 12-char geoid_bg
# on every component, then merges. Seconds to run; needs only the caches.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr) })
base_dir <- "/Users/priyanka/MIGRATE_dir"; poll_dir <- file.path(base_dir,"pollution_bg")
cache_dir <- file.path(poll_dir,"cache_multi"); if (!dir.exists(cache_dir)) cache_dir <- file.path(poll_dir,"cache")
padbg <- function(dt){ dt[, geoid_bg := str_pad(as.character(geoid_bg),12,pad="0")]; dt }

read_comp <- function(prefix, valname){
  fs <- list.files(cache_dir, pattern=sprintf("^%s_[0-9]{4}\\.csv$", prefix), full.names=TRUE)
  if (!length(fs)) { message("no cache for ", prefix); return(NULL) }
  d <- rbindlist(lapply(fs, function(f) fread(f, colClasses=c(geoid_bg="character"))), fill=TRUE)
  padbg(d); setnames(d, setdiff(names(d), c("geoid_bg","year")), valname)   # value col -> canonical name
  unique(d[!is.na(geoid_bg)], by=c("geoid_bg","year"))
}
pm  <- read_comp("pm25","pm25")
no2 <- read_comp("no2","no2")
oz  <- read_comp("ozone","ozone")
wf_f <- file.path(cache_dir,"wildfire_pm25_ALL.csv")
wf  <- if (file.exists(wf_f)) unique(padbg(fread(wf_f, colClasses=c(geoid_bg="character"))), by=c("geoid_bg","year")) else NULL

parts <- Filter(Negate(is.null), list(pm,no2,wf,oz))
panel <- Reduce(function(a,b) merge(a,b,by=c("geoid_bg","year"),all=TRUE), parts)
panel <- unique(panel, by=c("geoid_bg","year"))
setcolorder(panel, c("geoid_bg","year", intersect(c("pm25","no2","wildfire_pm25","ozone"), names(panel))))
fwrite(panel, file.path(poll_dir,"bg_pollution_panel.csv"))

# report coverage to confirm the fix
.p <- intersect(c("pm25","no2","wildfire_pm25","ozone"), names(panel))
cat("Rebuilt bg_pollution_panel.csv:", nrow(panel), "BG-years (was ~3.40M with duplicates)\n")
cov <- panel[year %in% 2011:2019, c(list(year=year), lapply(.SD, function(x) round(100*mean(!is.na(x)),1))), .SDcols=.p][
  , lapply(.SD, mean), by=year][order(year)]
cat("PM2.5 coverage by year (should be ~100% for all, incl. 2011-2014):\n"); print(cov[, .(year, pm25=round(pm25,1))])
ca <- panel[substr(geoid_bg,1,2)=="06" & year==2012]
ca_cov <- 100*mean(!is.na(ca$pm25))
cat(sprintf("California 2012: %d rows, %d unique BGs, pm25 non-missing %.0f%%\n",
            nrow(ca), uniqueN(ca$geoid_bg), ca_cov))
# HARD FAIL if the padding fix did not take (CA / low-FIPS states still dropped)
if (nrow(panel) != uniqueN(panel[, .(geoid_bg, year)]))
  stop("REBUILD FAILED: duplicate BG-year rows remain in the panel.")
if (ca_cov < 99)
  stop(sprintf("REBUILD FAILED: California 2012 PM2.5 coverage is %.0f%% (expected ~100%%). The geoid padding fix did not take.", ca_cov))
cat("REBUILD OK: no duplicate BG-years; low-FIPS states (incl. CA) fully covered.\n")
