# ============================================================
# POLLUTION JOIN DIAGNOSTICS + NATIONAL MAPS
# ------------------------------------------------------------
# Confirms every join worked and visualizes each pollutant:
#   1. Per-pollutant, per-year block-group COVERAGE (% non-missing).
#   2. Value sanity spot-checks at known cities.
#   3. FLOW join: share of MIGRATE origins/destinations matched to the
#      pollution panel (and to the race panel) for a reference year.
#   4. County-aggregated national choropleth maps of every pollutant.
#
# Output: pm25_avoidance_lag3_olag1/diagnostics_coverage.csv,
#         .../figures/diag_coverage_heatmap.png,
#         .../figures/maps/map_<pollutant>.png
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr); library(sf); library(tigris); library(ggplot2) })
options(tigris_use_cache=TRUE); sf::sf_use_s2(FALSE)
base_dir <- "/Users/priyanka/MIGRATE_dir"; mig_dir <- base_dir
out_dir <- file.path(base_dir,"pm25_avoidance_lag3_olag1"); fig_dir <- file.path(out_dir,"figures"); map_dir <- file.path(fig_dir,"maps")
dir.create(map_dir, showWarnings=FALSE, recursive=TRUE)
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")
MAPYR <- 2016
PLAB <- c(pm25="PM2.5 (ug/m3)", no2="NO2", wildfire_pm25="Wildfire PM2.5 (ug/m3)", ozone="Ozone WSMDA8 (ppb)",
          so4="Sulfate", no3="Nitrate", bc="Black carbon", om="Organic matter")

panel <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character")); panel[, geoid_bg:=pad_bg(geoid_bg)]
comp_file <- file.path(base_dir,"pollution_bg","bg_pm25_components_panel.csv")
if (file.exists(comp_file)) { cp<-fread(comp_file,colClasses=c(geoid_bg="character")); cp[,geoid_bg:=pad_bg(geoid_bg)]
  cc<-intersect(c("so4","no3","bc","om"),names(cp)); panel<-merge(panel,cp[,c("geoid_bg","year",cc),with=FALSE],by=c("geoid_bg","year"),all.x=TRUE) }
ndup <- nrow(panel) - nrow(unique(panel, by=c("geoid_bg","year")))
if (ndup>0) { message("Note: ", ndup, " duplicate BG-year rows in panel (deduplicating)."); panel <- unique(panel, by=c("geoid_bg","year")) }
polls <- intersect(names(PLAB), names(panel))
message("Pollutants in panel: ", paste(polls, collapse=", "))

# ---- 1. coverage per pollutant per year ----
cov <- panel[, c(list(n_bg=.N), lapply(.SD, function(x) round(100*mean(!is.na(x)),1))), by=year, .SDcols=polls][order(year)]
fwrite(cov, file.path(out_dir,"diagnostics_coverage.csv"))
cat("\n=== Block-group coverage (% non-missing) by pollutant x year ===\n"); print(cov)

# ---- 2. value sanity at known locations (nearest BG by centroid not available here; report national quantiles) ----
cat("\n=== Value ranges (pooled; sanity check units/magnitudes) ===\n")
for (p in polls) { v<-panel[[p]]; cat(sprintf("  %-14s mean %.2f  p05 %.2f  p95 %.2f  min %.2f  max %.2f\n",
    p, mean(v,na.rm=TRUE), quantile(v,.05,na.rm=TRUE), quantile(v,.95,na.rm=TRUE), min(v,na.rm=TRUE), max(v,na.rm=TRUE))) }

# ---- 3. FLOW join diagnostic (reference year) ----
bgset <- unique(panel[year==MAPYR & !is.na(pm25)]$geoid_bg)
f <- file.path(mig_dir, paste0(MAPYR,".csv"))
if (file.exists(f)) {
  x <- fread(f, select=c("Origin","Destination","Flow"))
  x[, `:=`(o=pad_bg(Origin), d=pad_bg(Destination))]
  cat(sprintf("\n=== Flow join (%d): origins matched to panel %.1f%%, destinations %.1f%% (flow-weighted %.1f%% / %.1f%%) ===\n",
      MAPYR, 100*mean(x$o %in% bgset), 100*mean(x$d %in% bgset),
      100*x[o %in% bgset, sum(Flow)]/x[,sum(Flow)], 100*x[d %in% bgset, sum(Flow)]/x[,sum(Flow)]))
}
rc <- list.files(file.path(base_dir,"pollution_bg"), sprintf("race_bg_%d.csv$",MAPYR), full.names=TRUE)
if (length(rc)) { r<-fread(rc[1],colClasses=c(geoid_bg="character")); r[,geoid_bg:=pad_bg(geoid_bg)]
  cat(sprintf("=== Race join (%d): panel BGs matched to race panel %.1f%% ===\n", MAPYR, 100*mean(bgset %in% r$geoid_bg))) }

# ---- coverage heatmap ----
cl <- melt(cov[, c("year",polls), with=FALSE], id.vars="year", variable.name="pollutant", value.name="pct")
ph <- ggplot(cl, aes(factor(year), pollutant, fill=pct)) + geom_tile(color="white") +
  geom_text(aes(label=round(pct)), size=2.6) +
  scale_fill_gradient(low="#fdd0a2", high="#238b45", limits=c(0,100), name="% BG covered") +
  labs(title="Block-group coverage by pollutant and year", x=NULL, y=NULL) + theme_minimal(base_size=10)
ggsave(file.path(fig_dir,"diag_coverage_heatmap.png"), ph, width=9, height=3.6, dpi=200)

# ---- 4. national county maps of each pollutant (reference year) ----
message("Building county map base ...")
co <- tigris::counties(cb=TRUE, year=2019)
co <- co[!co$STATEFP %in% c("02","15","60","66","69","72","78"), ]      # CONUS
co <- st_transform(st_make_valid(co), 5070); co$county <- str_pad(co$GEOID,5,pad="0")
pm <- panel[year==MAPYR]; pm[, county := substr(geoid_bg,1,5)]
cmean <- pm[, lapply(.SD, function(x) mean(x,na.rm=TRUE)), by=county, .SDcols=polls]
com <- merge(co, cmean, by="county", all.x=TRUE)
for (p in polls) {
  gg <- ggplot(com) + geom_sf(aes(fill=.data[[p]]), color=NA) +
    scale_fill_viridis_c(option="magma", direction=-1, na.value="grey90", name=unname(PLAB[p])) +
    labs(title=paste0(PLAB[p], " by county, ", MAPYR), subtitle="County means of block-group extractions (join diagnostic)") +
    theme_void(base_size=11) + theme(plot.title=element_text(face="bold"))
  ggsave(file.path(map_dir, paste0("map_",p,".png")), gg, width=8, height=5, dpi=200)
}
message("\nDONE. Diagnostics + ", length(polls), " pollutant maps in ", map_dir)
