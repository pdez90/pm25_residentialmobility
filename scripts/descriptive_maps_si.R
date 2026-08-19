# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# DESCRIPTIVE MAPS FOR THE SI: start vs end year
# ------------------------------------------------------------
# National maps of every quantity the paper uses, at the start and
# end of its data window, so the reader can see the raw surfaces:
#   1. All 8 pollutants, 2011 vs 2019, at raw BLOCK-GROUP resolution
#      (~240k polygons per panel; the slowest step, ~2-5 min per map).
#   2. ACS covariates used in Eq. 2 and affordability, at BLOCK-GROUP
#      level: median household income, gross rent, home value,
#      % bachelor's, % poverty (2013 vs 2019; block-group ACS starts 2013).
#   3. Race/ethnicity at BLOCK-GROUP level: % Black, % Hispanic, % white,
#      and the majority category used in the models, 2011 vs 2019.
#   4. Gentrification: county share of eligible block groups, and share
#      of eligible BGs flagged gentrifying during 2011-2019 (4 measures).
#   5. Segregation: CBSA Black-white dissimilarity, 2011 vs 2019.
#   6. Migration-PUMA maps (ACS replication geography), 2012 vs 2019:
#      in-mover volume, % Black and % Hispanic among in-movers, and
#      MIGPUMA mean PM2.5.
#
# Sankeys are NOT re-made here: pollution_sankey.R already wrote them.
#
# Output: pm25_avoidance_lag3_olag1/figures/descriptive_maps/*.png
# Runtime: ~10-20 min (tigris downloads cached; IPUMS read is the slow step).
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(stringr); library(sf); library(tigris)
  library(ggplot2); library(scales)
})
options(tigris_use_cache = TRUE); sf::sf_use_s2(FALSE)

base_dir <- migrate_data_dir()
poll_dir <- file.path(base_dir, "pollution_bg")
gent_dir <- file.path(base_dir, "acs_gentrification")
ipums_dir <- file.path(base_dir, "ipums")
out_dir  <- file.path(base_dir, "pm25_avoidance_lag3_olag1")
fig_dir  <- file.path(out_dir, "figures", "descriptive_maps")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

pad_bg  <- function(x) str_pad(as.character(x), 12, pad = "0")
YRS_POLL <- c(2011, 2019); YRS_ACS <- c(2013, 2019); YRS_PUMA <- c(2012, 2019)
noncontig <- c("02","15","60","66","69","72","78")   # CONUS only
ALBERS <- 5070

PLAB <- c(pm25 = "PM2.5 (µg/m³)", no2 = "NO2 (ppb)",
          ozone = "Warm-season ozone (ppb)", wildfire_pm25 = "Wildfire smoke PM2.5 (µg/m³)",
          so4 = "Sulfate (µg/m³)", no3 = "Nitrate (µg/m³)",
          bc = "Black carbon (µg/m³)", om = "Organic matter (µg/m³)")

# ---------- shared county geometry ----------
cty <- counties(cb = TRUE, year = 2019)
cty <- cty[!cty$STATEFP %in% noncontig, c("GEOID")]
cty <- st_transform(cty, ALBERS)

theme_map <- theme_void(base_size = 11) +
  theme(legend.position = "bottom",
        legend.key.width = unit(1.4, "cm"), legend.key.height = unit(.3, "cm"),
        plot.title = element_text(face = "bold", hjust = .5),
        strip.text = element_text(face = "bold", size = 11))

# draw a 2-panel (year-faceted) county map with a shared, outlier-trimmed scale
map2yr <- function(dat, valcol, legend, title, file, trans = "identity",
                   pal = "viridis", digits = 1) {
  lim <- quantile(dat[[valcol]], c(.02, .98), na.rm = TRUE)
  dat$v <- pmin(pmax(dat[[valcol]], lim[1]), lim[2])
  m <- merge(cty, dat[, c("GEOID", "year", "v")], by = "GEOID")
  gg <- ggplot(m) + geom_sf(aes(fill = v), color = NA) +
    facet_wrap(~year, ncol = 2) +
    scale_fill_viridis_c(option = pal, name = legend, trans = trans,
                         labels = label_number(accuracy = 10^-digits),
                         na.value = "grey85") +
    labs(title = title) + theme_map
  ggsave(file.path(fig_dir, file), gg, width = 12, height = 4.8, dpi = 200)
  message("wrote ", file)
}

# ---------- block-group -> county aggregation helper ----------
# pop-weighted county mean of a BG variable, per year
to_county <- function(dt, valcol, wcol = "w") {
  d <- dt[!is.na(get(valcol)) & !is.na(get(wcol)) & get(wcol) > 0]
  out <- d[, .(v = weighted.mean(get(valcol), get(wcol))),
           by = .(GEOID = substr(geoid_bg, 1, 5), year)]
  setnames(out, "v", valcol); out[]
}

# ============================================================
# 1. POLLUTANTS, 2011 vs 2019
# ============================================================
panel <- fread(file.path(poll_dir, "bg_pollution_panel.csv"),
               colClasses = c(geoid_bg = "character"))
panel[, geoid_bg := pad_bg(geoid_bg)]
comp_file <- file.path(poll_dir, "bg_pm25_components_panel.csv")
if (file.exists(comp_file)) {
  cp <- fread(comp_file, colClasses = c(geoid_bg = "character"))
  cp[, geoid_bg := pad_bg(geoid_bg)]
  cc <- intersect(c("so4","no3","bc","om"), names(cp))
  panel <- merge(panel, cp[, c("geoid_bg","year",cc), with = FALSE],
                 by = c("geoid_bg","year"), all.x = TRUE)
}
panel <- unique(panel[year %in% YRS_POLL], by = c("geoid_bg","year"))
polls <- intersect(names(PLAB), names(panel))

# population weights (ACS total_pop; nearest available year for 2011)
gfile <- file.path(gent_dir, "blockgroup_gentrification_2009_2019_lag3_with_HMDA.csv")
gpop <- fread(gfile, select = c("year","geoid_bg","total_pop"),
              colClasses = c(geoid_bg = "character"))
gpop[, geoid_bg := pad_bg(geoid_bg)]
gpop <- unique(gpop, by = c("geoid_bg","year"))
wt <- function(yr) { yy <- gpop[year == ifelse(yr < 2013, 2013, yr)]
  setnames(yy[, .(geoid_bg, w = total_pop)], c("geoid_bg","w")) }
for (yr in YRS_POLL) panel[year == yr, w := wt(yr)[panel[year == yr], on = "geoid_bg", x.w]]

# block-group geometry: cartographic-boundary files, 2010-vintage BG
# definitions (matching the panel geoids). ~49 state downloads, cached.
message("loading block-group geometry (slow on first run; cached by tigris)")
fc <- unique(tigris::fips_codes$state_code)
fc <- fc[fc <= "56" & !fc %in% noncontig]
bg_sf <- do.call(rbind, lapply(fc, function(s)
  block_groups(state = s, cb = TRUE, year = 2019)[, "GEOID"]))
bg_sf <- st_transform(bg_sf, ALBERS)
message("block groups loaded: ", nrow(bg_sf))

# generic two-year block-group map (unit of analysis in the paper)
map2yr_bg <- function(dat, valcol, legend, title, file, trans = "identity") {
  dat <- dat[!is.na(get(valcol))]
  lim <- quantile(dat[[valcol]], c(.02, .98))
  dat[, v := pmin(pmax(get(valcol), lim[1]), lim[2])]
  m <- merge(bg_sf, dat[, .(GEOID = geoid_bg, year, v)], by = "GEOID")
  gg <- ggplot(m) + geom_sf(aes(fill = v), color = NA) +
    facet_wrap(~year, ncol = 2) +
    scale_fill_viridis_c(name = legend, trans = trans, na.value = "grey85") +
    labs(title = title) + theme_map
  ggsave(file.path(fig_dir, file), gg, width = 12, height = 4.8, dpi = 200)
  rm(m); gc(verbose = FALSE)
  message("wrote ", file, " (block-group level)")
}

for (p in polls)
  map2yr_bg(panel, p, PLAB[[p]],
            paste0(sub(" \\(.*", "", PLAB[[p]]), ", block-group annual mean"),
            paste0("map_", p, "_2011_2019.png"))

# ============================================================
# 2. ACS COVARIATES, 2013 vs 2019
# ============================================================
acs <- rbindlist(lapply(YRS_ACS, function(yr)
  fread(file.path(gent_dir, paste0("acs_bg_", yr, ".csv")),
        colClasses = c(geoid_bg = "character"))), fill = TRUE)
acs[, geoid_bg := pad_bg(geoid_bg)]
acs[, `:=`(bach_pct = 100 * bach_plus / pmax(total_25plus, 1),
           poverty_pct = 100 * poverty_num / pmax(poverty_den, 1),
           w = total_pop)]
ACSV <- c(med_hh_income = "Median household income ($)",
          med_gross_rent = "Median gross rent ($/month)",
          med_home_value = "Median home value ($)",
          bach_pct = "Bachelor's degree or higher (%)",
          poverty_pct = "Population below poverty (%)")
for (v in names(ACSV))
  map2yr_bg(acs, v, ACSV[[v]], paste0(ACSV[[v]], ", block-group value"),
            paste0("map_acs_", v, "_2013_2019.png"),
            trans = if (v %in% c("med_hh_income","med_home_value")) "log10" else "identity")

# ============================================================
# 3. RACE / ETHNICITY, 2011 vs 2019
# ============================================================
rc <- rbindlist(lapply(YRS_POLL, function(yr) {
  r <- fread(file.path(poll_dir, paste0("race_bg_", yr, ".csv")),
             colClasses = c(geoid_bg = "character")); r[, year := yr]; r }), fill = TRUE)
rc[, geoid_bg := pad_bg(geoid_bg)]
for (yr in YRS_POLL) rc[year == yr, w := wt(yr)[rc[year == yr], on = "geoid_bg", x.w]]
for (v in c("sh_black","sh_hisp","sh_white")) rc[, (v) := 100 * get(v)]
RLAB <- c(sh_black = "Black population share (%)",
          sh_hisp = "Hispanic population share (%)",
          sh_white = "White population share (%)")
for (v in names(RLAB))
  map2yr_bg(rc, v, RLAB[[v]], paste0(RLAB[[v]], ", block-group value"),
            paste0("map_race_", v, "_2011_2019.png"))
# categorical majority map at BLOCK-GROUP level (the paper's origin categories)
rc[, maj := fifelse(sh_black > 50, "Majority Black",
             fifelse(sh_hisp > 50, "Majority Hispanic",
              fifelse(sh_white > 50, "Majority white", "No majority")))]
m <- merge(bg_sf, rc[!is.na(maj), .(GEOID = geoid_bg, year, maj)], by = "GEOID")
gg <- ggplot(m) + geom_sf(aes(fill = maj), color = NA) + facet_wrap(~year, ncol = 2) +
  scale_fill_manual(name = "Block-group racial/ethnic majority",
    values = c("Majority Black" = "#7b3294", "Majority Hispanic" = "#e66101",
               "Majority white" = "#80cdc1", "No majority" = "#f7e5a4")) +
  labs(title = "Neighborhood racial/ethnic majority (block group)") + theme_map
ggsave(file.path(fig_dir, "map_race_majority_2011_2019.png"), gg, width = 12, height = 4.8, dpi = 200)
rm(m); gc(verbose = FALSE)
message("wrote map_race_majority_2011_2019.png (block-group level)")

# ============================================================
# 4. GENTRIFICATION (window-based: eligibility + ever-gentrifying 2011-2019)
# ============================================================
ghdr <- names(fread(gfile, nrows = 0))
flags <- intersect(c("gentrifying_transition_top10","gentrifying_expectation_top10",
                     "gentrifying_classic_binary",
                     if ("gentrifying_capital_top10_screened" %in% ghdr)
                       "gentrifying_capital_top10_screened" else "gentrifying_capital_top10"), ghdr)
g <- fread(gfile, select = c("year","geoid_bg","eligible_gentrifiable", flags),
           colClasses = c(geoid_bg = "character"))
g <- g[year %in% 2011:2019]; g[, geoid_bg := pad_bg(geoid_bg)]
ge <- g[, c(.(eligible = as.integer(any(eligible_gentrifiable == 1, na.rm = TRUE))),
            lapply(.SD, function(x) as.integer(any(x == 1, na.rm = TRUE)))),
        by = geoid_bg, .SDcols = flags]
ge[, GEOID := substr(geoid_bg, 1, 5)]
gc <- ge[, c(.(elig_share = 100 * mean(eligible)),
             lapply(.SD, function(x) 100 * mean(x[eligible == 1]))),
         by = GEOID, .SDcols = flags]
gl <- melt(gc, id.vars = "GEOID")
GLAB <- c(elig_share = "Eligible for gentrification (% of BGs)",
          gentrifying_transition_top10 = "Gentrifying: transition (% of eligible)",
          gentrifying_expectation_top10 = "Gentrifying: market expectations (% of eligible)",
          gentrifying_classic_binary = "Gentrifying: classic Freeman/Bates (% of eligible)",
          gentrifying_capital_top10 = "Gentrifying: capital reinvestment (% of eligible)",
          gentrifying_capital_top10_screened = "Gentrifying: capital reinvestment (% of eligible)")
gl[, panel := factor(GLAB[as.character(variable)], levels = unique(GLAB))]
m <- merge(cty, gl, by = "GEOID")
gg <- ggplot(m) + geom_sf(aes(fill = pmin(value, 100)), color = NA) +
  facet_wrap(~panel, ncol = 2) +
  scale_fill_viridis_c(name = "% of block groups", limits = c(0, 100)) +
  labs(title = "Gentrification eligibility and gentrifying block groups, 2011–2019") + theme_map
ggsave(file.path(fig_dir, "map_gentrification_2011_2019.png"), gg, width = 12, height = 9, dpi = 200)
message("wrote map_gentrification_2011_2019.png")

# ============================================================
# 5. SEGREGATION: CBSA Black-white dissimilarity, 2011 vs 2019
# ============================================================
xw <- fread(file.path(poll_dir, "county_cbsa_xwalk.csv"),
            colClasses = "character")
xw <- xw[!grepl("^nonmetro", cbsa)]
seg <- rc[!is.na(w) & w > 0,
          .(black = sum(sh_black / 100 * w), white = sum(sh_white / 100 * w)),
          by = .(GEOID_bg = geoid_bg, county = substr(geoid_bg, 1, 5), year)]
seg <- merge(seg, xw, by = "county")
dis <- seg[, .(black = black / sum(black), white = white / sum(white)),
           by = .(cbsa, year, GEOID_bg)][
       , .(D = 0.5 * sum(abs(black - white), na.rm = TRUE)), by = .(cbsa, year)]
cb <- core_based_statistical_areas(cb = TRUE, year = 2019)
cb <- st_transform(cb[, c("GEOID")], ALBERS)
st <- states(cb = TRUE, year = 2019); st <- st_transform(st[!st$STATEFP %in% noncontig, "STATEFP"], ALBERS)
m <- merge(cb, dis, by.x = "GEOID", by.y = "cbsa")
m <- st_intersection(m, st_union(st))
gg <- ggplot() + geom_sf(data = st, fill = "grey92", color = "white", linewidth = .2) +
  geom_sf(data = m, aes(fill = D), color = NA) + facet_wrap(~year, ncol = 2) +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       name = "Black–white dissimilarity index", limits = c(0, 1)) +
  labs(title = "Metropolitan Black–white residential segregation (Duncan dissimilarity)") + theme_map
ggsave(file.path(fig_dir, "map_segregation_dissimilarity_2011_2019.png"), gg,
       width = 12, height = 4.8, dpi = 200)
message("wrote map_segregation_dissimilarity_2011_2019.png")

# ============================================================
# 6. MIGRATION-PUMA MAPS (ACS replication geography), 2012 vs 2019
# ============================================================
# MIGPUMA polygons: dissolve 2010-definition PUMAs via the same crosswalk
# used in acs_micro_validation.R
pu <- pumas(cb = TRUE, year = 2019)   # cb_2019 national file, 2010 PUMA definitions
pu <- pu[!pu$STATEFP10 %in% noncontig, c("STATEFP10","PUMACE10")]
pu <- st_transform(pu, ALBERS)
mx <- fread(file.path(ipums_dir, "puma_migpuma_xwalk.csv"), colClasses = "character")
pu <- merge(pu, mx, by.x = c("STATEFP10","PUMACE10"), by.y = c("statefip","puma"))
pu$mig_id <- paste0(pu$STATEFP10, "_", pu$migpuma)
gl <- split(st_geometry(pu), pu$mig_id)
mp <- st_sf(mig_id = names(gl),
            geometry = st_sfc(do.call(c, lapply(gl, st_union)), crs = st_crs(pu)))

# movers from IPUMS (same filters as acs_micro_validation.R)
ipf <- file.path(ipums_dir, "usa_acs_migration.csv")
ac <- fread(ipf, select = c("YEAR","PERWT","STATEFIP","PUMA","MIGRATE1",
                            "MIGPLAC1","MIGPUMA1","RACE","HISPAN"))
ac <- ac[YEAR %in% YRS_PUMA & MIGRATE1 %in% c(2,3) & MIGPLAC1 >= 1 & MIGPLAC1 <= 56 & MIGPUMA1 > 0]
ac[, race3 := fifelse(HISPAN > 0 & HISPAN < 9, "Hispanic",
              fifelse(RACE == 2, "Black", fifelse(RACE == 1, "White", "Other")))]
# destination MIGPUMA of the current residence
ac[, `:=`(statefip = str_pad(STATEFIP, 2, pad = "0"), puma = str_pad(PUMA, 5, pad = "0"))]
ac <- merge(ac, mx, by = c("statefip","puma"))
ac[, mig_id := paste0(statefip, "_", migpuma)]
pm <- ac[, .(inflow = sum(PERWT),
             sh_black = 100 * sum(PERWT[race3 == "Black"]) / sum(PERWT),
             sh_hisp  = 100 * sum(PERWT[race3 == "Hispanic"]) / sum(PERWT)),
         by = .(mig_id, year = YEAR)]

# MIGPUMA mean PM2.5 (pop-weighted BG centroids -> PUMA -> MIGPUMA)
cen <- fread(file.path(poll_dir, "bg_centroids.csv"), colClasses = c(geoid_bg = "character"))
cen[, geoid_bg := pad_bg(geoid_bg)]
cs <- st_as_sf(cen, coords = c("lon","lat"), crs = 4326); cs <- st_transform(cs, ALBERS)
cs$puma_i <- as.integer(st_intersects(cs, pu, sparse = TRUE))
cen[, mig_id := ifelse(is.na(cs$puma_i), NA, paste0(pu$STATEFP10, "_", pu$migpuma)[cs$puma_i])]
pp <- merge(panel[year %in% YRS_POLL, .(geoid_bg, year, pm25, w)],
            cen[, .(geoid_bg, mig_id)], by = "geoid_bg")
pp <- pp[!is.na(mig_id) & !is.na(w) & w > 0,
         .(pm25 = weighted.mean(pm25, w, na.rm = TRUE)), by = .(mig_id, year)]

pmap <- function(dat, valcol, legend, title, file, trans = "identity") {
  lim <- quantile(dat[[valcol]], c(.02, .98), na.rm = TRUE)
  dat$v <- pmin(pmax(dat[[valcol]], lim[1]), lim[2])
  m <- merge(mp, dat[, c("mig_id","year","v")], by = "mig_id")
  gg <- ggplot(m) + geom_sf(aes(fill = v), color = NA) + facet_wrap(~year, ncol = 2) +
    scale_fill_viridis_c(name = legend, trans = trans, na.value = "grey85") +
    labs(title = title) + theme_map
  ggsave(file.path(fig_dir, file), gg, width = 12, height = 4.8, dpi = 200)
  message("wrote ", file)
}
pmap(pm, "inflow", "In-movers (weighted persons)", "ACS in-movers by migration PUMA",
     "map_migpuma_inflow_2012_2019.png", trans = "log10")
pmap(pm, "sh_black", "% of in-movers", "Black share of in-movers, by migration PUMA",
     "map_migpuma_black_inmovers_2012_2019.png")
pmap(pm, "sh_hisp", "% of in-movers", "Hispanic share of in-movers, by migration PUMA",
     "map_migpuma_hispanic_inmovers_2012_2019.png")
pmap(pp, "pm25", "PM2.5 (µg/m³)", "Mean PM2.5 by migration PUMA",
     "map_migpuma_pm25_2012_2019.png")

message("\nDONE. All descriptive maps in ", fig_dir)
