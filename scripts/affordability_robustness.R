# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# METHODOLOGICAL HARDENING: affordability, selection, measurement
# error, pollutant collinearity, and mean reversion
# ------------------------------------------------------------
# Addresses the reviewer's methodological concerns with public data:
#   A. AFFORDABILITY decomposition. Does the racial avoidance gap
#      operate through reaching cheaper (dirtier) destinations? Gelbach
#      decomposition: how much of the Black/Hispanic gap is explained by
#      destination housing cost (ACS median home value / gross rent).
#   B. PROBABILITY of reaching a cleaner neighborhood as a function of
#      origin race, destination housing cost, and move distance.
#   C. SELECTION into moving. Re-estimate the conditional gap on
#      within-CBSA moves and on short-distance moves only, and adding
#      move distance as a control.
#   D. MEASUREMENT ERROR. Re-estimate weighting by inverse PM2.5 variance
#      (WashU AbsoluteUncertainty), down-weighting high-uncertainty BGs.
#   E. POLLUTANT COLLINEARITY. Correlation matrix of the pollutants
#      (they are expected to correlate: all mark combustion infrastructure).
#   F. MEAN REVERSION. A figure and a simulation showing why the RAW
#      racial gap is baseline-driven and reverses under adjustment.
#
# Housing cost and race come from public ACS; PM2.5 uncertainty from the
# V6.NA AbsoluteUncertainty product. Estimand note: contrasts are at the
# level of neighborhood context (streams leaving majority-X neighborhoods),
# not individual race.
#
# Output: pm25_avoidance_lag3_olag1/  (CSVs + figures)
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr); library(fixest); library(splines); library(ggplot2); library(sf) })
sf::sf_use_s2(FALSE)
base_dir <- migrate_data_dir(); mig_dir <- base_dir
gent_dir <- file.path(base_dir,"acs_gentrification"); change_lag <- 3
out_dir <- file.path(base_dir,"pm25_avoidance_lag3_olag1"); fig_dir <- file.path(out_dir,"figures")
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
migrate_years <- 2011:2019; SHORT_KM <- 50
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")
wmean  <- function(x,w){ok<-!is.na(x)&!is.na(w); if(!any(ok)) return(NA_real_); sum(x[ok]*w[ok])/sum(w[ok])}
hav <- function(lo1,la1,lo2,la2){ R<-6371; dlo<-(lo2-lo1)*pi/180; dla<-(la2-la1)*pi/180
  a<-sin(dla/2)^2+cos(la1*pi/180)*cos(la2*pi/180)*sin(dlo/2)^2; 2*R*asin(pmin(1,sqrt(a))) }
bcoef <- function(fit,term="grpMajority Black") tryCatch(coeftable(fit)[term,"Estimate"], error=function(e) NA_real_)

# ---- panel + uncertainty ----
poll <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character")); poll[, geoid_bg:=pad_bg(geoid_bg)]
poll <- unique(poll, by=c("geoid_bg","year"))
unc_file <- file.path(base_dir,"pollution_bg","bg_pm25_uncertainty_panel.csv")
if (file.exists(unc_file)) { un<-fread(unc_file,colClasses=c(geoid_bg="character")); un[,geoid_bg:=pad_bg(geoid_bg)]
  poll <- merge(poll, unique(un,by=c("geoid_bg","year")), by=c("geoid_bg","year"), all.x=TRUE) } else message("No PM2.5 uncertainty panel; skipping measurement-error robustness.")

# ---- composition + housing + gentrification ----
hdr <- names(fread(file.path(gent_dir, paste0("blockgroup_gentrification_2009_2019_lag",change_lag,"_with_HMDA.csv")), nrows=0))
cf  <- file.path(gent_dir, paste0("blockgroup_gentrification_2009_2019_lag",change_lag,"_with_HMDA.csv"))
acs <- fread(cf, colClasses=c(geoid_bg="character", state_county="character"),
  select=unique(c("year","geoid_bg","state_county","med_hh_income","poverty_pct","bach_pct","total_pop",
                  "med_home_value","med_gross_rent")))
acs[, geoid_bg:=pad_bg(geoid_bg)]
rcs <- list.files(file.path(base_dir,"pollution_bg"), "^race_bg_.*csv$", full.names=TRUE)
rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}")); r<-fread(f,colClasses=c(geoid_bg="character")); r[,year:=y]; r }), fill=TRUE)
rc[, geoid_bg:=pad_bg(geoid_bg)]
acs <- merge(acs, rc[,.(year,geoid_bg,sh_black,sh_hisp,sh_white)], by=c("year","geoid_bg"), all.x=TRUE)
acs[, race_majority := fifelse(sh_black>0.5,"Majority Black", fifelse(sh_hisp>0.5,"Majority Hispanic", fifelse(sh_white>0.5,"Majority white (NH)","No majority")))]
acs <- merge(acs, poll, by=c("geoid_bg","year"), all.x=TRUE)
acs[, county := substr(geoid_bg,1,5)]
xw <- file.path(base_dir,"pollution_bg","county_cbsa_xwalk.csv")
if (file.exists(xw)) { cxw<-fread(xw,colClasses=c(county="character",cbsa="character")); acs<-merge(acs,unique(cxw,by="county"),by="county",all.x=TRUE); acs[is.na(cbsa),cbsa:=paste0("nonmetro_",county)] } else acs[, cbsa:=county]
acs <- unique(acs, by=c("year","geoid_bg"))

# ---- BG centroids for move distance ----
bg_cache <- file.path(base_dir,"pollution_bg","cache_multi","bg_national.rds")
if (!file.exists(bg_cache)) bg_cache <- file.path(base_dir,"pollution_bg","cache","bg_national.rds")
cent <- data.table(); if (file.exists(bg_cache)) {
  bg <- readRDS(bg_cache); cc <- sf::st_coordinates(sf::st_centroid(bg)); cent <- data.table(geoid_bg=pad_bg(bg$geoid_bg), lon=cc[,1], lat=cc[,2])
} else message("BG centroid cache not found; distance-based robustness will be skipped.")

# ---- origin / dest lookups ----
oa <- acs[, .(cy=year, origin_bg=geoid_bg, o_race=race_majority, o_inc0=med_hh_income, o_pov=poverty_pct,
              o_bach=bach_pct, o_county=county, o_cbsa=cbsa, o_pm25=pm25,
              o_pm25_se=if("pm25_se"%in%names(acs)) pm25_se else NA_real_)]
if (nrow(cent)) oa <- merge(oa, cent[,.(origin_bg=geoid_bg,o_lon=lon,o_lat=lat)], by="origin_bg", all.x=TRUE)
setkey(oa, cy, origin_bg)
da <- acs[, .(cy=year, dest_bg=geoid_bg, d_pm25=pm25, d_cbsa=cbsa, d_homeval=med_home_value, d_rent=med_gross_rent,
              d_pm25_se=if("pm25_se"%in%names(acs)) pm25_se else NA_real_)]
if (nrow(cent)) da <- merge(da, cent[,.(dest_bg=geoid_bg,d_lon=lon,d_lat=lat)], by="dest_bg", all.x=TRUE)
setkey(da, cy, dest_bg)

# ---- MIGRATE pass: origin-BG-year aggregates (all / within-CBSA / short) ----
acc <- list()
for (yr in migrate_years) {
  x <- fread(file.path(mig_dir, paste0(yr,".csv")), select=c("Origin","Destination","Flow"))
  x <- x[, .(origin_bg=pad_bg(Origin), dest_bg=pad_bg(Destination), flow=as.numeric(Flow))][flow>0 & origin_bg!=dest_bg]
  x[, cy:=yr]; setkey(x,cy,origin_bg); x<-oa[x,on=.(cy,origin_bg)]; setkey(x,cy,dest_bg); x<-da[x,on=.(cy,dest_bg)]
  x[, re := d_pm25 - o_pm25]; x[, cleaner := as.integer(re < 0)]
  x[, within := as.integer(o_cbsa==d_cbsa)]
  x[, dist := if ("o_lon"%in%names(x)) hav(o_lon,o_lat,d_lon,d_lat) else NA_real_]
  acc[[as.character(yr)]] <- x[, .(
    year=yr, out_flow=sum(flow,na.rm=TRUE),
    re_pm25 = wmean(re,flow), orig_pm25 = wmean(o_pm25,flow),
    dest_homeval = wmean(d_homeval,flow), dest_rent = wmean(d_rent,flow),
    mean_dist = wmean(dist,flow), share_cleaner = wmean(cleaner,flow),
    o_se = data.table::first(o_pm25_se), dest_se = wmean(d_pm25_se,flow),
    out_flow_wc = sum(flow[within==1],na.rm=TRUE), re_wc = wmean(re[within==1], flow[within==1]),
    out_flow_sd = sum(flow[!is.na(dist)&dist<=SHORT_KM],na.rm=TRUE), re_sd = wmean(re[!is.na(dist)&dist<=SHORT_KM], flow[!is.na(dist)&dist<=SHORT_KM])
  ), by=.(origin_bg, o_race, o_inc0, o_pov, o_bach, o_county, o_cbsa)]
}
oy <- rbindlist(acc, fill=TRUE)
oy <- oy[!is.na(o_race) & out_flow>0 & !is.na(o_inc0) & o_inc0>0 & !is.na(orig_pm25)]
oy[, grp := relevel(factor(o_race), ref="Majority white (NH)")]
fwrite(oy, file.path(out_dir,"origin_year_affordability.csv"))
ctrl <- "+ ns(orig_pm25,4) + log(o_inc0) + o_pov + o_bach"

# ============================================================
# A. AFFORDABILITY: sequential covariate attenuation on a COMMON sample.
#    This is a Gelbach-STYLE descriptive diagnostic, NOT causal mediation:
#    destination housing cost is a post-move, chosen quantity. All models use
#    the SAME sample (dest_homeval>0 & dest_rent>0) so the change in the racial
#    coefficient reflects the added covariates, not a shifting sample.
# ============================================================
m0 <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," | o_cbsa + year")), oy, weights=~out_flow, cluster=~o_cbsa)  # full-sample baseline (reused below)
oy_aff <- oy[dest_homeval>0 & dest_rent>0]                                                                        # common sample for the attenuation
a0 <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," | o_cbsa + year")), oy_aff, weights=~out_flow, cluster=~o_cbsa)
a1 <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," + log(dest_homeval) | o_cbsa + year")), oy_aff, weights=~out_flow, cluster=~o_cbsa)
a2 <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," + log(dest_homeval) + log(dest_rent) | o_cbsa + year")), oy_aff, weights=~out_flow, cluster=~o_cbsa)
bse <- function(fit,t) tryCatch(coeftable(fit)[t,"Std. Error"], error=function(e) NA_real_)
afford <- rbindlist(lapply(c("grpMajority Black","grpMajority Hispanic"), function(t){
  b0<-bcoef(a0,t); b1<-bcoef(a1,t); b2<-bcoef(a2,t)
  data.table(group=sub("grpMajority ","",t), gap_base=b0, se_base=bse(a0,t),
             gap_plus_homeval=b1, se_homeval=bse(a1,t), gap_plus_rent=b2, se_rent=bse(a2,t),
             pct_attenuation_homeval=round(100*(b0-b1)/b0,1), pct_attenuation_housing=round(100*(b0-b2)/b0,1)) }), fill=TRUE)
fwrite(afford, file.path(out_dir,"affordability_decomposition.csv"))
cat("\n=== AFFORDABILITY (sequential attenuation on a common sample; descriptive, not mediation) ===\n"); print(afford)
cat("  common sample:", nrow(oy_aff), "of", nrow(oy), "origin-years\n")

# ============================================================
# A2. GELBACH DECOMPOSITION (order-invariant, exact).
#    The sequential differences above give the TOTAL attenuation but split it
#    in a path-dependent way. Gelbach (2016, J. Labor Econ. 34:509) attributes
#    that total to each covariate in an ORDER-INVARIANT way via the identity
#         b_base - b_full  =  sum_k  delta_k * gamma_k
#    gamma_k = covariate k's coefficient in the FULL model; delta_k = the
#    coefficient on the race dummy when covariate k is regressed on the BASE
#    design (race + controls | FE). All regressions use the SAME common sample,
#    controls, FE, and flow weights, so by Frisch-Waugh-Lovell the contributions
#    sum EXACTLY to the total attenuation (self-checked below). Still DESCRIPTIVE:
#    destination housing cost is a post-move, chosen quantity, not a treatment.
# ============================================================
X3   <- c("log(dest_homeval)","log(dest_rent)")                                  # covariates being attributed
tmap <- c(`grpMajority Black`="Black", `grpMajority Hispanic`="Hispanic", `grpNo majority`="No majority")
gam  <- sapply(X3, function(v) tryCatch(coeftable(a2)[v,"Estimate"], error=function(e) NA_real_))  # gamma_k (full model)
aux  <- lapply(X3, function(v)                                                    # auxiliary regressions: X3_k on the base design
          feols(as.formula(paste0(v," ~ grp",ctrl," | o_cbsa + year")), oy_aff, weights=~out_flow, cluster=~o_cbsa))
names(aux) <- X3
gelb <- rbindlist(lapply(names(tmap), function(t){
  b0 <- bcoef(a0,t); b2 <- bcoef(a2,t)                                            # base vs full gap on the COMMON sample
  delt  <- sapply(X3, function(v) bcoef(aux[[v]], t))                             # delta_k for this race term
  contr <- delt * gam                                                            # per-covariate contribution (order-invariant)
  data.table(group=tmap[[t]], gap_base=b0, gap_full=b2, total_explained=b0-b2,
             contrib_homeval=contr[["log(dest_homeval)"]], contrib_rent=contr[["log(dest_rent)"]],
             pct_homeval=round(100*contr[["log(dest_homeval)"]]/b0,2),
             pct_rent   =round(100*contr[["log(dest_rent)"]]/b0,2),
             pct_total  =round(100*sum(contr)/b0,2),
             identity_gap=(b0-b2)-sum(contr)) }), fill=TRUE)                      # must be ~0 (FWL)
fwrite(gelb, file.path(out_dir,"gelbach_decomposition.csv"))
cat("\n=== GELBACH DECOMPOSITION (order-invariant; per-covariate shares sum to total attenuation) ===\n"); print(gelb)
# HARD CHECK: the additivity identity is exact under FWL. Fail loudly if it isn't
# (that would signal a sample/weight/spec mismatch between the base, full, and
# auxiliary regressions -- not a real decomposition).
maxgap <- max(abs(gelb$identity_gap), na.rm=TRUE)
cat(sprintf("  additivity check: max|(b_base - b_full) - sum_k delta_k*gamma_k| = %.3e\n", maxgap))
if (!is.finite(maxgap) || maxgap > 1e-5)
  stop(sprintf("Gelbach additivity identity violated (max gap %.3e). The base (a0), full (a2), and auxiliary regressions must share the SAME sample, weights, controls, and fixed effects.", maxgap))

# ============================================================
# B. PROBABILITY of reaching a cleaner neighborhood
# ============================================================
pc0 <- feols(as.formula(paste0("share_cleaner ~ grp",ctrl," | o_cbsa + year")), oy, weights=~out_flow, cluster=~o_cbsa)
pc1 <- feols(as.formula(paste0("share_cleaner ~ grp",ctrl," + log(dest_homeval) + mean_dist | o_cbsa + year")), oy[dest_homeval>0], weights=~out_flow, cluster=~o_cbsa)
prob <- rbindlist(lapply(c("grpMajority Black","grpMajority Hispanic"), function(t){
  data.table(group=sub("grpMajority ","",t), gap_base=bcoef(pc0,t), gap_plus_afford=bcoef(pc1,t),
             pct_explained=round(100*(bcoef(pc0,t)-bcoef(pc1,t))/bcoef(pc0,t),1)) }), fill=TRUE)
fwrite(prob, file.path(out_dir,"prob_cleaner_decomposition.csv"))
cat("\n=== PROBABILITY of reaching cleaner air: racial gap and share explained by affordability+distance ===\n"); print(prob)

# ============================================================
# C. SELECTION robustness: within-CBSA, short-distance, +distance control
# ============================================================
rob <- list(all=m0)
rob$within_cbsa <- feols(as.formula(paste0("re_wc ~ grp",ctrl," | o_cbsa + year")), oy[out_flow_wc>0], weights=~out_flow_wc, cluster=~o_cbsa)
if ("mean_dist" %in% names(oy) && any(!is.na(oy$mean_dist))) {
  rob$short_dist <- feols(as.formula(paste0("re_sd ~ grp",ctrl," | o_cbsa + year")), oy[out_flow_sd>0], weights=~out_flow_sd, cluster=~o_cbsa)
  rob$plus_dist  <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," + mean_dist | o_cbsa + year")), oy, weights=~out_flow, cluster=~o_cbsa)
}
seltab <- rbindlist(lapply(names(rob), function(n){
  ct <- as.data.table(coeftable(rob[[n]]), keep.rownames="term")[grepl("^grpMajority",term)]
  ct[, .(spec=n, group=sub("grpMajority ","",term), est=Estimate, se=`Std. Error`, p=`Pr(>|t|)`)] }), fill=TRUE)
fwrite(seltab, file.path(out_dir,"selection_robustness.csv"))
cat("\n=== SELECTION robustness: conditional gap across move subsamples ===\n"); print(seltab)

# ============================================================
# D. UNCERTAINTY-WEIGHTED SENSITIVITY (NOT a formal measurement-error
#    correction: AbsoluteUncertainty is not an ordinary sampling standard
#    error, origin and destination errors may be spatially correlated, and the
#    variance of a difference is not the sum of the two variances. This simply
#    down-weights high-uncertainty flows as a robustness check.)
# ============================================================
if (all(c("o_se","dest_se") %in% names(oy)) && any(!is.na(oy$o_se))) {
  oy[, ivw := out_flow / pmax(o_se^2 + dest_se^2, 1e-6)]
  me <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," | o_cbsa + year")), oy[!is.na(ivw)], weights=~ivw, cluster=~o_cbsa)
  metab <- data.table(spec=c("flow-weighted (main)","uncertainty-weighted (sensitivity)"),
                      Black=c(bcoef(m0),bcoef(me)),
                      Black_se=c(bse(m0,"grpMajority Black"),bse(me,"grpMajority Black")),
                      Hispanic=c(bcoef(m0,"grpMajority Hispanic"),bcoef(me,"grpMajority Hispanic")),
                      Hispanic_se=c(bse(m0,"grpMajority Hispanic"),bse(me,"grpMajority Hispanic")))
  fwrite(metab, file.path(out_dir,"measurement_error_robustness.csv"))
  cat("\n=== UNCERTAINTY-WEIGHTED SENSITIVITY: conditional gap, flow- vs inverse-uncertainty weighting ===\n"); print(metab)
}

# ============================================================
# E. POLLUTANT COLLINEARITY
# ============================================================
pc <- intersect(c("pm25","no2","wildfire_pm25","ozone"), names(poll))
comp <- file.path(base_dir,"pollution_bg","bg_pm25_components_panel.csv")
if (file.exists(comp)) { cp<-fread(comp,colClasses=c(geoid_bg="character")); cp[,geoid_bg:=pad_bg(geoid_bg)]
  cc2<-intersect(c("so4","no3","bc","om"),names(cp)); poll<-merge(poll,cp[,c("geoid_bg","year",cc2),with=FALSE],by=c("geoid_bg","year"),all.x=TRUE); pc<-c(pc,cc2) }
cormat <- round(cor(poll[year==2016, ..pc], use="pairwise.complete.obs"), 2)
fwrite(as.data.table(cormat, keep.rownames="pollutant"), file.path(out_dir,"pollutant_correlations.csv"))
cat("\n=== POLLUTANT correlation matrix (2016) ===\n"); print(cormat)

# ============================================================
# F. MEAN REVERSION: figure + simulation
# ============================================================
oy[, opm_bin := cut(orig_pm25, quantile(orig_pm25, 0:10/10, na.rm=TRUE), include.lowest=TRUE)]
mrb <- oy[o_race %in% c("Majority Black","Majority white (NH)","Majority Hispanic"),
          .(avoid = -wmean(re_pm25,out_flow), obin_mid = wmean(orig_pm25,out_flow)), by=.(o_race,opm_bin)]
pmr <- ggplot(mrb[!is.na(opm_bin)], aes(obin_mid, avoid, color=o_race)) + geom_line(linewidth=1) + geom_point(size=1.4) +
  geom_hline(yintercept=0, linewidth=.3) +
  scale_color_manual(values=c(`Majority Black`="#1b7837",`Majority Hispanic`="#762a83",`Majority white (NH)`="#2166ac"), name=NULL) +
  labs(title="Why the raw racial gap largely reflects mean reversion",
       subtitle="Within bins of origin PM2.5, avoidance is nearly identical across origin race; the raw gap comes from different baselines",
       x=expression(Origin~PM[2.5]~", "~mu*g/m^3), y=expression(PM[2.5]~reduction~by~moving~", "~mu*g/m^3)) +
  theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"), legend.position="bottom")
ggsave(file.path(fig_dir,"figS_mean_reversion_explained.png"), pmr, width=8, height=5, dpi=200)

# STYLIZED ILLUSTRATION (not evidence): regression toward a COMMON long-run
# mean makes higher-baseline groups show larger raw reductions, so a raw racial
# gap can appear even when true avoidance is identical. (Reverting to each
# group's OWN mean, as before, would make both groups' expected change equal
# and would NOT reproduce a baseline-driven gap.)
set.seed(1)
common_target <- 8.0
sim <- function(base, n=5e4, true_avoid=0.05, revert=0.30, sd=1.0){
  o <- rnorm(n, base, sd); d <- o - true_avoid - revert*(o - common_target) + rnorm(n,0,sd); data.table(o,d,re=d-o) }
sB <- sim(8.8); sW <- sim(7.9)
cat(sprintf("\n=== MEAN-REVERSION (stylized illustration; common target %.1f, identical true avoidance 0.05) ===\n  Raw avoidance: higher-baseline group %.3f vs lower-baseline group %.3f (spurious gap %.3f)\n  Illustrative only — not evidence that the empirical result is mean reversion.\n",
    common_target, -mean(sB$re), -mean(sW$re), -mean(sB$re)-(-mean(sW$re))))

cat("\nDONE. Robustness outputs + figures in ", out_dir, "\n", sep="")
