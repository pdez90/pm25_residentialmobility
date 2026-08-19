# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# TARGETED COVERAGE DIAGNOSTIC (pre-2015 PM2.5)
# ------------------------------------------------------------
# The coverage table shows PM2.5 at 82% of block groups in 2006-2014 and
# 100% from 2015. This checks whether that 18% gap matters for the
# 2011-2019 analysis, two ways:
#   (1) POPULATION-WEIGHTED coverage by year -- if missing BGs are
#       unpopulated, pop-weighted coverage stays ~100% and results are
#       unaffected.
#   (2) Coverage by STATE for a low year -- flags spatial clustering
#       (regional bias) vs. a scatter of tiny/empty BGs.
# Instant (no MIGRATE pass).
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr) })
base_dir <- migrate_data_dir()
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")
poll <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character"))
poll[, geoid_bg:=pad_bg(geoid_bg)]; poll <- unique(poll[,.(geoid_bg,year,pm25)], by=c("geoid_bg","year"))
pop <- fread(file.path(base_dir,"acs_gentrification","blockgroup_gentrification_2009_2019_lag3_with_HMDA.csv"),
             colClasses=c(geoid_bg="character"), select=c("year","geoid_bg","total_pop")); pop[, geoid_bg:=pad_bg(geoid_bg)]
d <- merge(poll, pop, by=c("geoid_bg","year"), all.x=TRUE); d[is.na(total_pop), total_pop:=0]
d[, state:=substr(geoid_bg,1,2)]

# (1) BG-count vs population-weighted coverage, by year
cov <- d[year %in% 2011:2019, .(
  bg_cov  = round(100*mean(!is.na(pm25)),1),
  pop_cov = round(100*sum(total_pop*(!is.na(pm25)))/sum(total_pop),1)), by=year][order(year)]
cat("\n=== PM2.5 coverage 2011-2019: BG-count vs POPULATION-weighted ===\n")
cat("(if pop_cov ~100 while bg_cov ~82, the missing BGs are unpopulated -> analysis unaffected)\n")
print(cov)

# (2) coverage by state, 2012 (a low-coverage year) -- flag clustering
st <- d[year==2012, .(bg_cov=round(100*mean(!is.na(pm25)),1),
                      pop_cov=round(100*sum(total_pop*(!is.na(pm25)))/sum(total_pop),1),
                      n_bg=.N), by=state][order(bg_cov)]
cat("\n=== 2012 PM2.5 coverage by STATE (lowest first; watch for whole states missing) ===\n")
print(st[1:15])
cat("\nStates with <90% BG coverage in 2012:\n")
print(st[bg_cov<90, .(state, bg_cov, pop_cov, n_bg)])
cat("\nInterpretation: whole low-coverage states = regional gap to address;",
    "\nscattered BGs with pop_cov~100 = benign (unpopulated polygons).\n")
