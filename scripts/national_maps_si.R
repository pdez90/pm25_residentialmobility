# ============================================================
# SUPPLEMENTARY NATIONAL MAPS
# Choropleth maps of every quantity used in the paper:
#   (A) the four gentrification measures - gentrifying rate among
#       eligible BGs, and (for the three continuous measures) the
#       mean score percentile;
#   (B) origin opportunity percentile;
#   (C) every destination-context variable in bg_context_table.csv
#       (opportunity, social capital, health, environment, ADI,
#        tree canopy, park cover, amenity access, walkability,
#        segregation - whatever is present in the table).
#
# DEFAULT = county level (population-weighted BG means; fast,
# legible, standard SI format). Set map_level <- "bg" for full
# block-group detail (slower: ~240k polygons per map).
#
# CONUS only (Albers EPSG:5070). Reads the same reference year as
# the static context table so maps and analysis are consistent.
# Self-contained: run directly or via master_run_all.R.
# ============================================================
library(data.table)
library(sf)
library(stringr)
library(ggplot2)
if (!requireNamespace("tigris", quietly = TRUE)) install.packages("tigris")

# ----------------------------
# Paths & switches
# ----------------------------
base_dir <- "/Users/priyanka/MIGRATE_dir"
gent_dir <- file.path(base_dir, "acs_gentrification")
change_lag <- 3
origin_lag <- 1
map_level <- "bg"              # "bg" (default, block-group) or "county"
map_year <- 2016               # reference year for measures/opportunity
lag_tag <- paste0("lag", change_lag)
olag_tag <- paste0("olag", origin_lag)
tag <- paste0(lag_tag, "_", olag_tag)
ctx_tab <- file.path(base_dir, paste0("migration_destination_context_", tag),
                     "tables", "bg_context_table.csv")
fig_dir <- file.path(base_dir, paste0("migration_destination_context_", tag),
                     "figures", "national_maps_si")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
options(tigris_use_cache = TRUE)
noncontig <- c("02", "15", "60", "66", "69", "72", "78")  # AK/HI/territories
pad_bg <- function(x) str_pad(str_remove_all(as.character(x), "\\D"), 12, pad = "0")

# ----------------------------
# Labels (measures + context). Missing -> raw column name.
# ----------------------------
labs <- c(
  # measures / opportunity
  rate_transition = "Gentrifying rate: neighborhood transition",
  rate_expectations = "Gentrifying rate: market expectations",
  rate_classic = "Gentrifying rate: classic Freeman/Bates",
  rate_capital = "Gentrifying rate: capital reinvestment",
  score_transition = "Transition score (mean percentile)",
  score_expectations = "Market-expectation gap (mean percentile)",
  score_capital = "Capital-reinvestment score (mean percentile)",
  eligible_share = "Share of block groups gentrification-eligible",
  opportunity_pctile = "Origin opportunity (mean percentile)",
  # context
  oa_kfr_p25 = "Upward mobility (kfr p25)",
  oa_jail_p25 = "Incarceration (p25)",
  sc_ec = "Economic connectedness",
  sc_support_ratio = "Support ratio",
  sc_volunteering = "Volunteering",
  places_mhlth = "Frequent mental distress (%)",
  places_obesity = "Obesity (%)",
  places_casthma = "Asthma (%)",
  places_phlth = "Poor physical health (%)",
  places_diabetes = "Diabetes (%)",
  places_bphigh = "High blood pressure (%)",
  places_health_burden = "Health burden (composite z)",
  nri_hazard_burden = "Climate hazard burden (composite z)",
  ej_environmental_burden = "Environmental burden (composite z)",
  ej_pm25 = "PM2.5",
  ej_ozone = "Ozone", ej_diesel = "Diesel PM", ej_traffic = "Traffic proximity",
  adi_natrank = "Area Deprivation Index (national rank)",
  tcc_canopy = "Tree canopy cover (%)",
  padus_park_pct = "Park cover (PAD-US-AR, %)",
  amen_park_dist = "Distance to parks (visitor-weighted)",
  amen_grocery_dist = "Distance to grocery",
  amen_healthcare_dist = "Distance to healthcare",
  amen_school_dist = "Distance to schools",
  amen_access_composite = "Amenity access composite (z)",
  walkability = "Walkability (NatWalkInd)",
  travel_distance = "Travel distance (visitor-weighted)",
  seg_experienced = "Experienced income segregation",
  seg_residential = "Residential income segregation"
)
lab_of <- function(v) if (v %in% names(labs)) labs[[v]] else v

# ----------------------------
# 1. Block-group values: measures + opportunity + context
# ----------------------------
combined <- file.path(
  gent_dir, paste0("blockgroup_gentrification_2009_2019_", lag_tag, "_with_HMDA.csv"))
stopifnot(file.exists(combined))
ghdr <- names(fread(combined, nrows = 0))
cap_flag <- if ("gentrifying_capital_top10_screened" %in% ghdr)
  "gentrifying_capital_top10_screened" else "gentrifying_capital_top10"
meas <- data.table(
  measure = c("transition", "expectations", "classic", "capital"),
  flag = c("gentrifying_transition_top10", "gentrifying_expectation_top10",
           "gentrifying_classic_binary", cap_flag),
  score = c("gent_score_transition_pctile", "expectation_gap_pctile",
            NA, "capital_reinvestment_score_pctile"))
gsel <- unique(c("year", "geoid_bg", "statefp", "countyfp", "total_pop",
                 "med_hh_income", "med_home_value", "bach_pct", "poverty_pct",
                 "eligible_gentrifiable", meas$flag,
                 na.omit(meas$score)))
gsel <- intersect(gsel, ghdr)
g <- fread(combined, select = gsel, colClasses = c(geoid_bg = "character"))
g <- g[year == map_year]
g[, geoid_bg := pad_bg(geoid_bg)]
# opportunity percentile (same construction as the models)
g[, opportunity_pctile := {
  zi <- scale(med_hh_income); zb <- scale(bach_pct)
  zv <- scale(med_home_value); zp <- scale(poverty_pct)
  s <- as.numeric(zi + zb + zv - zp)
  frank(s, ties.method = "average", na.last = "keep") / sum(!is.na(s))
}]
g[, county := paste0(str_pad(statefp, 2, pad = "0"),
                     str_pad(countyfp, 3, pad = "0"))]

ctx <- NULL
if (file.exists(ctx_tab)) {
  ctx <- fread(ctx_tab, colClasses = c(geoid_bg = "character"))
  ctx[, geoid_bg := pad_bg(geoid_bg)]
  drop <- intersect(c("geoid_tract"), names(ctx))
  if (length(drop)) ctx[, (drop) := NULL]
  ctx_vars <- setdiff(names(ctx), "geoid_bg")
  ctx_vars <- ctx_vars[sapply(ctx_vars, function(v) is.numeric(ctx[[v]]))]
  message("Context variables to map (", length(ctx_vars), "): ",
          paste(ctx_vars, collapse = ", "))
} else {
  ctx_vars <- character(0)
  message("bg_context_table.csv not found - context maps skipped. ",
          "(Run destination_context_extension.R first.)")
}

bg <- g[, .(geoid_bg, county, total_pop, opportunity_pctile)]
for (i in seq_len(nrow(meas))) {
  fl <- meas$flag[i]; sc <- meas$score[i]; m <- meas$measure[i]
  if (fl %in% names(g)) bg[, (paste0("flag_", m)) := g[[fl]]]
  if (!is.na(sc) && sc %in% names(g)) bg[, (paste0("score_", m)) := g[[sc]]]
}
bg[, eligible := as.integer(g$eligible_gentrifiable == TRUE)]
if (!is.null(ctx)) bg <- merge(bg, ctx, by = "geoid_bg", all.x = TRUE)

# ----------------------------
# 2. Aggregate to mapping unit
# ----------------------------
wmean <- function(x, w) {
  ok <- !is.na(x) & !is.na(w); if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}
measure_names <- meas$measure
score_names <- meas[!is.na(score), measure]
if (map_level == "county") {
  agg <- bg[, {
    out <- list(
      eligible_share = mean(eligible, na.rm = TRUE),
      opportunity_pctile = wmean(opportunity_pctile, total_pop))
    for (m in measure_names) {
      fc <- paste0("flag_", m)
      out[[paste0("rate_", m)]] <-
        if (fc %in% names(.SD)) mean(.SD[[fc]], na.rm = TRUE) else NA_real_
    }
    for (m in score_names) {
      scn <- paste0("score_", m)
      out[[scn]] <-
        if (scn %in% names(.SD)) wmean(.SD[[scn]], total_pop) else NA_real_
    }
    for (v in ctx_vars) out[[v]] <- wmean(.SD[[v]], total_pop)
    out
  }, by = county]
  unit_key <- "county"
} else {
  agg <- copy(bg)
  for (m in measure_names) {
    fc <- paste0("flag_", m)
    if (fc %in% names(agg)) setnames(agg, fc, paste0("rate_", m))
  }
  setnames(agg, "eligible", "eligible_share")
  unit_key <- "geoid_bg"
}

# ----------------------------
# 3. Geometry
# ----------------------------
if (map_level == "county") {
  geo <- tigris::counties(cb = TRUE, year = 2019, progress_bar = FALSE)
  geo <- geo[!geo$STATEFP %in% noncontig, ]
  geo$county <- paste0(str_pad(geo$STATEFP, 2, pad = "0"),
                       str_pad(geo$COUNTYFP, 3, pad = "0"))
  geo <- geo[, "county"]
} else {
  bg_cache <- file.path(base_dir, "bg_polygons_2010vintage.gpkg")
  if (file.exists(bg_cache)) {
    geo <- st_read(bg_cache, quiet = TRUE)
  } else {
    message("Building BG polygons via tigris (one-time; cached)")
    geo <- do.call(rbind, lapply(setdiff(c(state.abb, "DC"),
                                         c("AK", "HI")), function(st) {
      x <- tigris::block_groups(state = st, year = 2019, cb = TRUE,
                                progress_bar = FALSE)
      x[, "GEOID"]
    }))
    st_write(geo, bg_cache, quiet = TRUE)
  }
  geo$geoid_bg <- pad_bg(geo$GEOID)
  geo <- geo[!str_sub(geo$geoid_bg, 1, 2) %in% noncontig, "geoid_bg"]
}
geo <- st_transform(geo, 5070)                       # CONUS Albers
geo <- st_simplify(geo, dTolerance = if (map_level == "county") 500 else 300)
mdt <- merge(geo, agg, by = unit_key, all.x = TRUE)

# ----------------------------
# 4. Draw + save one map per quantity
# ----------------------------
theme_map <- theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        legend.position = "right",
        legend.key.height = unit(1.1, "lines"))
draw <- function(var, kind = c("quantile", "rate")) {
  kind <- match.arg(kind)
  if (!var %in% names(mdt)) return(invisible())
  x <- mdt[[var]]
  if (all(is.na(x))) return(invisible())
  d <- mdt
  if (kind == "rate") {
    p <- ggplot(d) +
      geom_sf(aes(fill = .data[[var]]), color = NA) +
      scale_fill_viridis_c(option = "magma", direction = -1,
                           labels = scales::percent_format(accuracy = 1),
                           na.value = "grey90", name = NULL) +
      labs(title = lab_of(var)) + theme_map
  } else {
    # decile bins on non-missing values -> readable legend
    br <- unique(quantile(x, probs = seq(0, 1, 0.1), na.rm = TRUE))
    if (length(br) < 3) return(invisible())
    d$bin <- cut(x, breaks = br, include.lowest = TRUE, dig.lab = 3)
    p <- ggplot(d) +
      geom_sf(aes(fill = bin), color = NA) +
      scale_fill_viridis_d(option = "viridis", na.value = "grey90",
                           name = "Decile", guide = guide_legend(reverse = TRUE)) +
      labs(title = lab_of(var),
           subtitle = "National deciles (grey = no data)") +
      theme_map + theme(plot.subtitle = element_text(hjust = 0.5, size = 9))
  }
  ggsave(file.path(fig_dir, paste0("map_", var, "_", map_level, "_", tag, ".png")),
         p, width = 9, height = 6, dpi = 220)
  message("  wrote map: ", var)
}

# measures (rate) + eligible share
for (m in measure_names) draw(paste0("rate_", m), "rate")
draw("eligible_share", "rate")
# opportunity + measure score percentiles (quantile)
draw("opportunity_pctile", "quantile")
for (m in score_names) draw(paste0("score_", m), "quantile")
# every context variable (quantile)
for (v in ctx_vars) draw(v, "quantile")

# ----------------------------
# 5. Combined 4-panel measures figure (gentrifying rate)
# ----------------------------
rate_vars <- paste0("rate_", measure_names)
rate_vars <- rate_vars[rate_vars %in% names(mdt)]
if (length(rate_vars) >= 2) {
  # sf keeps its geometry column sticky, so select only the data
  # column (its name varies: "geometry" from tigris, "geom" from a
  # GeoPackage cache) and let sf carry geometry along.
  long <- do.call(rbind, lapply(rate_vars, function(v) {
    dd <- mdt[, v]
    names(dd)[names(dd) == v] <- "rate"
    dd$measure <- lab_of(v)
    dd
  }))
  pm <- ggplot(long) +
    geom_sf(aes(fill = rate), color = NA) +
    scale_fill_viridis_c(option = "magma", direction = -1,
                         labels = scales::percent_format(accuracy = 1),
                         na.value = "grey90", name = "Gentrifying\nrate") +
    facet_wrap(~ measure, ncol = 2) +
    labs(title = "Gentrification prevalence by measure (share of eligible block groups)") +
    theme_map + theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(fig_dir, paste0("map_measures_panel_", map_level, "_", tag, ".png")),
         pm, width = 12, height = 8, dpi = 220)
  message("  wrote combined measures panel")
}

message("\nNational SI maps written to: ", fig_dir,
        "\nmap_level = ", map_level, " | reference year = ", map_year)
