# ============================================================
# BLOCK-GROUP AREA vs EXPOSURE-GRID RESOLUTION
# ------------------------------------------------------------
# Computes the land area of every block group in the study (those
# present in the pollution panel, CONUS) and compares it with the
# native resolution of the exposure products:
#   ~1 km grids  (PM2.5, components, NO2, ozone)  -> 1 km2 cells
#   10 km grid   (wildfire smoke)                 -> 100 km2 cells
# Reports counts and shares of BGs (and of population) with area
# < 1 km2, 1-10 km2, 10-100 km2, >= 100 km2, plus general summary
# statistics. Uses Census ALAND (land area), so results do not
# depend on projection choices.
#
# Output: pm25_avoidance_lag3_olag1/blockgroup_area_summary.csv
#         pm25_avoidance_lag3_olag1/blockgroup_area_stats.csv
#         (+ console summary)
# Runtime: fast; tigris boundaries already cached from the map script.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr); library(tigris); library(sf) })
options(tigris_use_cache = TRUE)

base_dir <- Sys.getenv("MIGRATE_DIR", unset="/Users/priyanka/MIGRATE_dir")
poll_dir <- file.path(base_dir, "pollution_bg")
gent_dir <- file.path(base_dir, "acs_gentrification")
out_dir  <- file.path(base_dir, "pm25_avoidance_lag3_olag1")
pad_bg <- function(x) str_pad(as.character(x), 12, pad = "0")
noncontig <- c("02","15","60","66","69","72","78")

# ---- study block groups: those in the pollution panel ----
panel <- fread(file.path(poll_dir, "bg_pollution_panel.csv"),
               select = c("geoid_bg","year","pm25"), colClasses = c(geoid_bg = "character"))
panel[, geoid_bg := pad_bg(geoid_bg)]
study_bg <- unique(panel[!is.na(pm25)]$geoid_bg)
message("Study block groups (in pollution panel): ", length(study_bg))

# ---- land area from Census boundary files (ALAND, m2) ----
fc <- unique(tigris::fips_codes$state_code)
fc <- fc[fc <= "56" & !fc %in% noncontig]
bg <- rbindlist(lapply(fc, function(s) {
  b <- block_groups(state = s, cb = TRUE, year = 2019)
  data.table(geoid_bg = b$GEOID, aland_m2 = as.numeric(b$ALAND))
}))
bg <- bg[geoid_bg %in% study_bg]
bg[, area_km2 := aland_m2 / 1e6]
message("Matched to boundary files: ", nrow(bg), " (",
        round(100 * nrow(bg) / length(study_bg), 1), "% of study BGs)")

# ---- population weights (2019 ACS) ----
pop <- fread(file.path(gent_dir, "acs_bg_2019.csv"),
             select = c("geoid_bg","total_pop"), colClasses = c(geoid_bg = "character"))
pop[, geoid_bg := pad_bg(geoid_bg)]
bg <- merge(bg, unique(pop, by = "geoid_bg"), by = "geoid_bg", all.x = TRUE)
bg[is.na(total_pop), total_pop := 0]

# ---- summary statistics ----
qs <- c(.01,.05,.10,.25,.50,.75,.90,.95,.99)
stats <- data.table(
  statistic = c("n_block_groups","mean_km2","sd_km2",
                paste0("p", 100*qs, "_km2"), "min_km2","max_km2",
                "pop_weighted_mean_km2","pop_weighted_median_km2"),
  value = c(nrow(bg), mean(bg$area_km2), sd(bg$area_km2),
            quantile(bg$area_km2, qs), min(bg$area_km2), max(bg$area_km2),
            weighted.mean(bg$area_km2, bg$total_pop),
            { o <- bg[order(area_km2)]; cw <- cumsum(o$total_pop)/sum(o$total_pop)
              o$area_km2[which(cw >= .5)[1]] }))
stats[, value := round(value, 3)]
fwrite(stats, file.path(out_dir, "blockgroup_area_stats.csv"))
cat("\n=== Block-group land area: summary statistics (km2) ===\n"); print(stats)

# ---- comparison with exposure-grid cell sizes ----
bg[, size_class := cut(area_km2, c(-Inf, 1, 10, 100, Inf),
    labels = c("< 1 km2 (smaller than a 1 km grid cell)",
               "1-10 km2",
               "10-100 km2 (smaller than the 10 km wildfire cell)",
               ">= 100 km2 (larger than the 10 km wildfire cell)"))]
cls <- bg[, .(n_bg = .N,
              pct_bg = round(100 * .N / nrow(bg), 1),
              population = sum(total_pop),
              pct_pop = round(100 * sum(total_pop) / sum(bg$total_pop), 1)),
          by = size_class][order(size_class)]
fwrite(cls, file.path(out_dir, "blockgroup_area_summary.csv"))
cat("\n=== Block-group area vs exposure-grid resolution ===\n"); print(cls)

# headline numbers for the text
cat(sprintf("\nBGs smaller than 1 km2:   %s (%.1f%% of BGs; %.1f%% of population)\n",
    format(bg[area_km2 < 1, .N], big.mark = ","),
    100 * bg[area_km2 < 1, .N] / nrow(bg),
    100 * bg[area_km2 < 1, sum(total_pop)] / sum(bg$total_pop)))
cat(sprintf("BGs smaller than 10 km2:  %s (%.1f%% of BGs; %.1f%% of population)\n",
    format(bg[area_km2 < 10, .N], big.mark = ","),
    100 * bg[area_km2 < 10, .N] / nrow(bg),
    100 * bg[area_km2 < 10, sum(total_pop)] / sum(bg$total_pop)))
cat(sprintf("BGs smaller than 100 km2 (one 10 km wildfire cell): %s (%.1f%%; %.1f%% of population)\n",
    format(bg[area_km2 < 100, .N], big.mark = ","),
    100 * bg[area_km2 < 100, .N] / nrow(bg),
    100 * bg[area_km2 < 100, sum(total_pop)] / sum(bg$total_pop)))
message("\nDONE. Wrote blockgroup_area_stats.csv and blockgroup_area_summary.csv to ", out_dir)
