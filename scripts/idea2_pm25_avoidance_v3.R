# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# IDEA 2 v3 (PM2.5): FULLY-REPRODUCIBLE FIGURES + AMBITIOUS SPECS
# ------------------------------------------------------------
# Everything the paper needs, generated from code:
#   A. Disparity trends, ABSOLUTE (levels) and RELATIVE (ratios), by
#      race, 2011-2019. Speaks to the "narrowing vs widening" debate:
#      Colmer et al. 2024 (absolute Black-White gap narrowing) vs
#      Jbaily et al. 2022 / Kerr et al. 2023 (relative / health-burden
#      disparities persisting or widening).
#   B. Shift-share decomposition of each group's exposure change AND of
#      the Black-White / Hispanic-White GAP into in-place cleanup and
#      population sorting.
#   C. Conditional avoidance models with THREE fixed-effects specs:
#      CBSA (metro) = MAIN; county and national = SI robustness.
#   D. Mean-reversion scatter.  E. Raw avoidance.  F. Distributional HIA.
# All figures are written by this script (ggplot2).
#
# WHY CBSA (not county) FE IS THE MAIN SPEC:
#   The metropolitan area is the realistic residential choice set, so
#   CBSA FE compares movers facing the same metro housing market while
#   absorbing metro-level confounds. COUNTY FE is stricter: it discards
#   the cross-county avoidance margin (moving to a cleaner county/suburb),
#   where much real PM2.5 reduction happens, so it is a conservative
#   lower bound that changes the estimand to WITHIN-county avoidance.
#   NATIONAL (year FE only) keeps all variation. Reporting all three
#   brackets the result; the privilege is robust across them.
#
# Output: $MIGRATE_DATA_DIR/pm25_avoidance_lag3_olag1/figures/
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(stringr); library(fixest); library(splines)
  library(ggplot2); library(sf); library(tigris)
})
for (p in c("ggrepel")) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
options(tigris_use_cache = TRUE); sf::sf_use_s2(FALSE)

base_dir <- migrate_data_dir()
tag <- "lag3_olag1"; migrate_years <- 2011:2019
out_dir <- file.path(base_dir, paste0("pm25_avoidance_", tag)); fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
# NOTE: PM2.5 mortality CRF and averted-death conversion removed; the design
# (summed multi-year move-events, adjusted group contrast, no exposure duration
# or age-specific rates) does not support attributable-mortality claims. We report
# only an exposure-reallocation signature (person-microgram) as context.
pad_bg <- function(x) str_pad(as.character(x), 12, pad = "0")
wmean  <- function(x,w){ok<-!is.na(x)&!is.na(w); if(!any(ok)) return(NA_real_); sum(x[ok]*w[ok])/sum(w[ok])}
theme_pub <- theme_minimal(base_size = 11) + theme(plot.title = element_text(face="bold"))
racecol <- c(Black="#1b7837", Hispanic="#762a83", White="#2166ac",
             `Majority Black`="#1b7837", `Majority Hispanic`="#762a83", `Majority white (NH)`="#2166ac", `No majority`="#999999")

# ---- build poll / acs / oy if not in memory ----
if (!exists("oy") || !exists("acs") || !exists("poll"))
  source(file.path(base_dir, "scripts", "idea2_pm25_avoidance.R"))

# ---- attach race shares to acs (from cached race pulls) ----
if (!all(c("sh_black","sh_hisp","sh_white") %in% names(acs))) {
  rcs <- list.files(file.path(base_dir,"pollution_bg"), "^race_bg_.*csv$", full.names=TRUE)
  rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}"))
        r<-fread(f, colClasses=c(geoid_bg="character")); r[,year:=y]; r }), fill=TRUE)
  rc[, geoid_bg := pad_bg(geoid_bg)]
  acs <- merge(acs, rc[, .(year, geoid_bg, sh_black, sh_hisp, sh_white)], by=c("year","geoid_bg"), all.x=TRUE)
}

# ============================================================
# County -> CBSA crosswalk (tigris), cached
# ============================================================
xw_file <- file.path(base_dir, "pollution_bg", "county_cbsa_xwalk.csv")
if (file.exists(xw_file)) {
  cxw <- fread(xw_file, colClasses = c(county="character", cbsa="character"))
} else {
  cb <- st_transform(tigris::core_based_statistical_areas(cb = TRUE, year = 2019), 5070)
  co <- st_transform(tigris::counties(cb = TRUE, year = 2019), 5070)
  cc <- st_join(st_centroid(co), cb["CBSAFP"])
  cxw <- data.table(county = str_pad(cc$GEOID, 5, pad="0"), cbsa = cc$CBSAFP)
  cxw[is.na(cbsa), cbsa := paste0("nonmetro_", county)]
  fwrite(cxw, xw_file)
}
oy[, o_county := str_pad(o_state_county, 5, pad="0")]
oy <- merge(oy, cxw, by.x="o_county", by.y="county", all.x=TRUE)
setnames(oy, "cbsa", "o_cbsa"); oy[is.na(o_cbsa), o_cbsa := paste0("nonmetro_", o_county)]

# ============================================================
# A. Disparity trends: ABSOLUTE levels and RELATIVE ratios
# ============================================================
ts <- acs[!is.na(pm25), .(
  Black = wmean(pm25, total_pop*sh_black),
  Hispanic = wmean(pm25, total_pop*sh_hisp),
  White = wmean(pm25, total_pop*sh_white)), by=year][order(year)][year %in% migrate_years]
ts[, `:=`(gap_BW = Black-White, gap_HW = Hispanic-White,
          ratio_BW = Black/White, ratio_HW = Hispanic/White)]
fwrite(ts, file.path(out_dir, "disparity_trends.csv"))
cat("\n=== Disparity trends (absolute gap and relative ratio) ===\n"); print(ts)

lev <- melt(ts[, .(year, Black, Hispanic, White)], id.vars="year", variable.name="group", value.name="pm25")
lev[, panel := "Absolute level (ug/m3)"]
rat <- melt(ts[, .(year, Black=ratio_BW, Hispanic=ratio_HW)], id.vars="year", variable.name="group", value.name="pm25")
rat[, `:=`(panel="Relative to white (ratio)")]
p1 <- ggplot(rbind(lev, rat, fill=TRUE), aes(year, pm25, color=group)) +
  geom_hline(data=data.table(panel="Relative to white (ratio)", y=1), aes(yintercept=y), linetype=3, color="grey60") +
  geom_line(linewidth=1) + geom_point(size=1.5) +
  facet_wrap(~panel, scales="free_y") + scale_color_manual(values=racecol) +
  labs(title="Absolute PM2.5 gap narrowed, but the relative disparity did not",
       subtitle="Left: population-weighted mean PM2.5 fell for all groups. Right: minority/white exposure RATIO held or rose.",
       x=NULL, y=NULL, color=NULL) + theme_pub
ggsave(file.path(fig_dir,"fig1_disparity_absolute_vs_relative.png"), p1, width=10, height=4.3, dpi=200)

# ============================================================
# B. Decomposition: each group + the GAP (cleanup vs sorting)
# ============================================================
decomp_group <- function(share) {
  d <- acs[year %in% range(migrate_years), .(geoid_bg, year, pm25, pop = total_pop*get(share))]
  w <- dcast(d, geoid_bg ~ year, value.var=c("pm25","pop")); setnames(w, names(w), sub("_","",names(w)))
  t0 <- min(migrate_years); t1 <- max(migrate_years); w <- w[stats::complete.cases(w)]
  P0<-sum(w[[paste0("pop",t0)]]); P1<-sum(w[[paste0("pop",t1)]])
  data.table(
    cleanup = sum((w[[paste0("pop",t0)]]/P0)*(w[[paste0("pm25",t1)]]-w[[paste0("pm25",t0)]])),
    sorting = sum((w[[paste0("pop",t1)]]/P1 - w[[paste0("pop",t0)]]/P0)*w[[paste0("pm25",t1)]]))
}
dg <- rbindlist(list(cbind(group="Black",decomp_group("sh_black")),
                     cbind(group="Hispanic",decomp_group("sh_hisp")),
                     cbind(group="White",decomp_group("sh_white"))), fill=TRUE)
dg[, total := cleanup+sorting]; fwrite(dg, file.path(out_dir,"exposure_decomposition.csv"))
# gap decomposition
gapdec <- rbindlist(list(
  data.table(gap="Black - White", cleanup=dg[group=="Black"]$cleanup-dg[group=="White"]$cleanup,
             sorting=dg[group=="Black"]$sorting-dg[group=="White"]$sorting),
  data.table(gap="Hispanic - White", cleanup=dg[group=="Hispanic"]$cleanup-dg[group=="White"]$cleanup,
             sorting=dg[group=="Hispanic"]$sorting-dg[group=="White"]$sorting)))
fwrite(gapdec, file.path(out_dir,"gap_decomposition.csv"))
cat("\n=== Change in the exposure GAP: cleanup vs sorting ===\n"); print(gapdec)

dgl <- melt(dg[, .(group, cleanup, sorting)], id.vars="group", variable.name="channel", value.name="dpm")
p2 <- ggplot(dgl, aes(group, dpm, fill=channel)) + geom_col(position="dodge") + geom_hline(yintercept=0, linewidth=.4) +
  scale_fill_manual(values=c(cleanup="#4d4d4d", sorting="#d6604d"), labels=c("In-place cleanup","Migration (sorting)")) +
  labs(title="In-place cleanup, not migration, drove the PM2.5 decline (2011-2019)",
       x=NULL, y=expression(Delta~PM[2.5]~", "~mu*g/m^3), fill=NULL) + theme_pub
ggsave(file.path(fig_dir,"fig2_decomposition.png"), p2, width=7.5, height=4.3, dpi=200)

# ============================================================
# C. Conditional avoidance: CBSA (main) + county + national (SI)
# ============================================================
run_cond <- function(dim_col, ref, fe) {
  d <- oy[!is.na(get(dim_col)) & out_flow>0 & !is.na(o_med_hh_income) & o_med_hh_income>0 & !is.na(orig_pm25)]
  d[, grp := relevel(factor(get(dim_col)), ref=ref)]
  fe_str <- switch(fe, cbsa="o_cbsa + year", county="o_county + year", none="year")
  clu <- if (fe=="none") ~o_county else as.formula(paste0("~", sub(" \\+ year","",fe_str)))
  f <- as.formula(paste0("re_pm25 ~ grp + ns(orig_pm25,4) + log(o_med_hh_income) + o_poverty_pct + o_bach_pct | ", fe_str))
  fit <- tryCatch(feols(f, data=d, weights=~out_flow, cluster=clu), error=function(e) NULL); if (is.null(fit)) return(NULL)
  ct <- as.data.table(coeftable(fit), keep.rownames="term")[grepl("^grp",term)]
  ct[, `:=`(dimension=dim_col, fe=fe)][]
}
specs <- CJ(fe=c("cbsa","county","none"))
cond <- rbindlist(lapply(specs$fe, function(fe) rbindlist(list(
  if ("o_race_majority" %in% names(oy)) run_cond("o_race_majority","Majority white (NH)",fe),
  run_cond("o_inc_tercile","High income",fe),
  run_cond("o_gent_group","Eligible, not gentrifying",fe)), fill=TRUE)), fill=TRUE)
fwrite(cond, file.path(out_dir,"conditional_models_all_fe.csv"))
cat("\n=== Conditional avoidance across FE specs (+ = reached dirtier air) ===\n"); print(cond[order(dimension,term,fe)])

cond[, `:=`(lo=Estimate-1.96*`Std. Error`, hi=Estimate+1.96*`Std. Error`, name=sub("^grp","",term))]
felab <- c(cbsa="CBSA FE (main)", county="County FE (SI)", none="National, year FE (SI)")
cond[, fe_lab := factor(felab[fe], levels=unname(felab))]
# MAIN figure: CBSA spec
p3 <- ggplot(cond[fe=="cbsa"], aes(Estimate, reorder(name, Estimate))) +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_point(aes(color=Estimate>0), size=2.8) + geom_errorbarh(aes(xmin=lo,xmax=hi), height=.2) +
  scale_color_manual(values=c(`TRUE`="#b2182b",`FALSE`="#2166ac"), guide="none") +
  labs(title="Conditional pollution-avoidance gap (net of mean reversion)",
       subtitle="Main spec: origin-PM2.5 spline + income, poverty, education; CBSA + year FE. + = reached dirtier air.",
       x=expression(Adjusted~Delta~destination~PM[2.5]~vs~reference~", "~mu*g/m^3), y=NULL) + theme_pub
ggsave(file.path(fig_dir,"fig3_conditional_main_cbsa.png"), p3, width=8, height=4.6, dpi=200)
# SI figure: all three specs
p3si <- ggplot(cond[dimension=="o_race_majority"], aes(Estimate, reorder(name,Estimate), color=fe_lab)) +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_point(position=position_dodge(.5), size=2.4) +
  geom_errorbarh(aes(xmin=lo,xmax=hi), height=.2, position=position_dodge(.5)) +
  scale_color_manual(values=c("#1b7837","#e08214","#756bb1")) +
  labs(title="Robustness: conditional race gap across fixed-effects specifications",
       x=expression(Adjusted~Delta~destination~PM[2.5]~", "~mu*g/m^3), y=NULL, color=NULL) + theme_pub
ggsave(file.path(fig_dir,"figS_conditional_all_fe.png"), p3si, width=8, height=4, dpi=200)

# ============================================================
# D-F. mean reversion, raw avoidance, distributional HIA
# ============================================================
sub_cols <- intersect(c("o_race_majority","o_inc_tercile","o_high_poverty","o_gent_group"), names(oy))
avoid <- rbindlist(lapply(sub_cols, function(g) {
  s <- oy[!is.na(get(g)), .(avoidance = -wmean(re_pm25,out_flow), orig_pm25=wmean(orig_pm25,out_flow),
                            delta=wmean(re_pm25,out_flow), flow=sum(out_flow,na.rm=TRUE)), by=g]
  setnames(s,g,"subgroup"); s[, dimension:=sub("^o_","",g)][] }), fill=TRUE)
avoid[, exposure_reallocated := delta * flow]   # person-microgram signature (- = net reduction); NOT averted deaths
fwrite(avoid, file.path(out_dir,"avoidance_signature_by_subgroup.csv"))

# D. mean reversion
mr <- avoid[dimension %in% c("race_majority","inc_tercile")]
p4 <- ggplot(mr, aes(orig_pm25, avoidance)) +
  geom_smooth(method="lm", se=FALSE, color="grey60", linetype=2) +
  geom_point(aes(size=flow, color=subgroup)) +
  ggrepel::geom_text_repel(aes(label=subgroup), size=3) +
  scale_color_manual(values=racecol, guide="none") + scale_size(guide="none") +
  geom_hline(yintercept=0, linewidth=.4) +
  labs(title="Avoidance tracks origin dirtiness (mean reversion)",
       subtitle="Points below the line under-capture avoidance for their baseline; marker size proportional to flow",
       x=expression(Origin~PM[2.5]~", "~mu*g/m^3~"(flow-weighted)"),
       y=expression(PM[2.5]~reduction~by~moving~", "~mu*g/m^3)) + theme_pub
tryCatch(ggsave(file.path(fig_dir,"fig4_mean_reversion.png"), p4, width=7, height=5, dpi=200),
         error=function(e) message("fig4 needs ggrepel: install.packages('ggrepel')"))

# E. raw avoidance
rawf <- avoid[dimension=="race_majority"]
p5 <- ggplot(rawf, aes(reorder(subgroup, avoidance), avoidance, fill=subgroup)) +
  geom_col() + scale_fill_manual(values=racecol, guide="none") + geom_hline(yintercept=0, linewidth=.4) +
  labs(title="Raw pollution avoidance by origin race (uncorrected)", x=NULL,
       y=expression(PM[2.5]~reduction~by~moving~", "~mu*g/m^3)) + coord_flip() + theme_pub
ggsave(file.path(fig_dir,"fig5_raw_avoidance.png"), p5, width=7, height=3.6, dpi=200)

# F. exposure-reallocation signature (context, NOT mortality)
sig <- copy(avoid); sig[, name := subgroup]
p6 <- ggplot(sig, aes(exposure_reallocated, reorder(name, exposure_reallocated), fill=exposure_reallocated<0)) +
  geom_col() + scale_fill_manual(values=c(`TRUE`="#2166ac",`FALSE`="#b2182b"), guide="none") +
  geom_vline(xintercept=0, linewidth=.5) + facet_grid(dimension~., scales="free_y", space="free_y") +
  labs(title="Exposure-reallocation signature by origin group (environmental-health context, not mortality)",
       x=expression("Exposure reallocated (person-"*mu*"g/m"^3*"; negative = net reduction)"), y=NULL) + theme_pub
ggsave(file.path(fig_dir,"fig6_exposure_signature.png"), p6, width=8, height=6, dpi=200)

cat("\nDONE. All figures written to ", fig_dir, "\n", sep="")
