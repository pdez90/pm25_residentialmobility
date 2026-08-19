# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# ENVIRONMENTAL OPPORTUNITY, HETEROGENEITY, and AFFORDABILITY+
# ------------------------------------------------------------
# Three reviewer-requested extensions, all on data in hand:
#  (2) MEASURE environmental opportunity. For each origin, the share of
#      clean (below-metro-median PM2.5) neighborhoods in its metro that
#      it could AFFORD (home value <= AFFORD_MULT x origin income). This
#      operationalizes "environmental opportunity" as a measured
#      opportunity SET rather than an inferred concept, tests whether it
#      is smaller for minority origins, and asks whether it explains the
#      avoidance gap.
#  (4) HETEROGENEITY. Does the avoidance gap scale with metropolitan
#      SEGREGATION (a testable prediction of the discrimination/steering
#      account), and does it hold across metros/regions?
#  (5) AFFORDABILITY+. Strengthen the null: add the origin-to-destination
#      CHANGE in housing cost, a rent-to-income measure, and the
#      opportunity-set share to the Gelbach decomposition.
#
# Inputs: origin_year_affordability.csv (stage 7), the pollution panel,
# the gentrification file, cached race, and the CBSA crosswalk.
# Output: pm25_avoidance_lag3_olag1/  (CSVs + figures)
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr); library(fixest); library(splines); library(ggplot2) })
base_dir <- migrate_data_dir()
out_dir <- file.path(base_dir,"pm25_avoidance_lag3_olag1"); fig_dir <- file.path(out_dir,"figures")
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
migrate_years <- 2011:2019; AFFORD_MULT <- 3
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")
wmean  <- function(x,w){ok<-!is.na(x)&!is.na(w); if(!any(ok)) return(NA_real_); sum(x[ok]*w[ok])/sum(w[ok])}
bcoef <- function(fit,t="grpMajority Black") tryCatch(coeftable(fit)[t,"Estimate"], error=function(e) NA_real_)

# ---- origin-year (realized moves) from stage 7 ----
oyf <- file.path(out_dir,"origin_year_affordability.csv")
stopifnot("run affordability_robustness.R (stage 7) first" = file.exists(oyf))
oy <- fread(oyf, colClasses=c(origin_bg="character")); oy[, origin_bg:=pad_bg(origin_bg)]

# ---- rebuild BG-year attributes (pollution + housing + race + cbsa) ----
poll <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character")); poll[, geoid_bg:=pad_bg(geoid_bg)]
poll <- unique(poll[, .(geoid_bg, year, pm25)], by=c("geoid_bg","year"))
cf <- file.path(base_dir,"acs_gentrification", "blockgroup_gentrification_2009_2019_lag3_with_HMDA.csv")
acs <- fread(cf, colClasses=c(geoid_bg="character"), select=c("year","geoid_bg","med_hh_income","med_home_value","med_gross_rent","total_pop"))
acs[, geoid_bg:=pad_bg(geoid_bg)]
rcs <- list.files(file.path(base_dir,"pollution_bg"), "^race_bg_.*csv$", full.names=TRUE)
rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}")); r<-fread(f,colClasses=c(geoid_bg="character")); r[,year:=y]; r }), fill=TRUE)
rc[, geoid_bg:=pad_bg(geoid_bg)]
acs <- merge(acs, rc[,.(year,geoid_bg,sh_black,sh_hisp,sh_white)], by=c("year","geoid_bg"), all.x=TRUE)
acs <- merge(acs, poll, by=c("geoid_bg","year"), all.x=TRUE)
acs[, county := substr(geoid_bg,1,5)]
xw <- fread(file.path(base_dir,"pollution_bg","county_cbsa_xwalk.csv"), colClasses=c(county="character",cbsa="character"))
acs <- merge(acs, unique(xw,by="county"), by="county", all.x=TRUE); acs[is.na(cbsa), cbsa:=paste0("nonmetro_",county)]
acs <- unique(acs, by=c("geoid_bg","year"))

# ============================================================
# (2) OPPORTUNITY SET: share of CLEAN metro neighborhoods AFFORDABLE to origin
# ============================================================
a2 <- acs[!is.na(pm25) & !is.na(med_home_value) & !is.na(med_hh_income) & med_hh_income>0]
a2[, opp_share := {
    mp <- median(pm25, na.rm=TRUE)
    ch <- sort(med_home_value[pm25 < mp])            # clean = below metro-median PM2.5
    if (length(ch)==0) rep(0, .N) else findInterval(med_hh_income*AFFORD_MULT, ch)/.N
  }, by=.(cbsa, year)]
opp <- a2[, .(origin_bg=geoid_bg, year, opp_share)]
oy <- merge(oy, opp, by=c("origin_bg","year"), all.x=TRUE)

# opportunity captured by origin race (flow-weighted)
oppr <- oy[!is.na(o_race), .(opp = wmean(opp_share, out_flow), inc = wmean(o_inc0, out_flow)), by=o_race]
fwrite(oppr, file.path(out_dir,"opportunity_by_race.csv"))
cat("\n=== ENVIRONMENTAL OPPORTUNITY captured (share of clean metro BGs affordable), by origin race ===\n"); print(oppr)

# ============================================================
# (4) SEGREGATION heterogeneity: does the gap scale with CBSA dissimilarity?
# ============================================================
acs[, `:=`(bpop=sh_black*total_pop, wpop=sh_white*total_pop)]
seg <- acs[cbsa %in% unique(oy$o_cbsa), .(D = 0.5*sum(abs(bpop/sum(bpop,na.rm=TRUE) - wpop/sum(wpop,na.rm=TRUE)), na.rm=TRUE)), by=.(cbsa,year)][
  , .(seg=mean(D,na.rm=TRUE)), by=cbsa]
oy <- merge(oy, seg, by.x="o_cbsa", by.y="cbsa", all.x=TRUE)
oy[, grp := relevel(factor(o_race), ref="Majority white (NH)")]
oy[, seg_z := as.numeric(scale(seg))]
ctrl <- "ns(orig_pm25,4) + log(o_inc0) + o_pov + o_bach"
# interaction: Black gap x segregation (year FE only; seg is CBSA-level)
mseg <- feols(as.formula(paste0("re_pm25 ~ grp*seg_z + ",ctrl," | year")), oy[!is.na(seg_z) & out_flow>0], weights=~out_flow, cluster=~o_cbsa)
segtab <- as.data.table(coeftable(mseg), keep.rownames="term")[grepl("grp",term)]
fwrite(segtab, file.path(out_dir,"segregation_interaction.csv"))
cat("\n=== SEGREGATION interaction (does the avoidance gap grow with metro segregation?) ===\n")
print(segtab[grepl("seg_z", term), .(term, est=round(Estimate,3), p=signif(`Pr(>|t|)`,2))])
# gap by segregation quartile
oy[, seg_q := cut(seg, quantile(seg, 0:4/4, na.rm=TRUE), labels=c("Q1 least","Q2","Q3","Q4 most"), include.lowest=TRUE)]
segq <- rbindlist(lapply(levels(oy$seg_q), function(q){
  d <- oy[seg_q==q & out_flow>0]; if(nrow(d)<50) return(NULL)
  fit <- tryCatch(feols(as.formula(paste0("re_pm25 ~ grp + ",ctrl," | o_cbsa + year")), d, weights=~out_flow, cluster=~o_cbsa), error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  data.table(seg_q=q, black_gap=bcoef(fit), se=tryCatch(coeftable(fit)["grpMajority Black","Std. Error"],error=function(e)NA)) }), fill=TRUE)
fwrite(segq, file.path(out_dir,"gap_by_segregation_quartile.csv"))
cat("\n=== Black avoidance gap by metro segregation quartile ===\n"); print(segq)

# ============================================================
# (5) AFFORDABILITY+: origin-dest housing change, rent-to-income, opportunity set
# ============================================================
oh <- acs[, .(origin_bg=geoid_bg, year, o_homeval=med_home_value)]
oy <- merge(oy, oh, by=c("origin_bg","year"), all.x=TRUE)
oy[, `:=`(hcost_change = dest_homeval - o_homeval, rent_to_inc = (dest_rent*12)/pmax(o_inc0,1))]
mods <- list(
  base   = feols(as.formula(paste0("re_pm25 ~ grp + ",ctrl," | o_cbsa + year")), oy[out_flow>0], weights=~out_flow, cluster=~o_cbsa),
  destHV = feols(as.formula(paste0("re_pm25 ~ grp + ",ctrl," + log(dest_homeval) | o_cbsa + year")), oy[out_flow>0 & dest_homeval>0], weights=~out_flow, cluster=~o_cbsa),
  hchg   = feols(as.formula(paste0("re_pm25 ~ grp + ",ctrl," + hcost_change | o_cbsa + year")), oy[out_flow>0 & !is.na(hcost_change)], weights=~out_flow, cluster=~o_cbsa),
  r2i    = feols(as.formula(paste0("re_pm25 ~ grp + ",ctrl," + rent_to_inc | o_cbsa + year")), oy[out_flow>0 & is.finite(rent_to_inc)], weights=~out_flow, cluster=~o_cbsa),
  opp    = feols(as.formula(paste0("re_pm25 ~ grp + ",ctrl," + opp_share | o_cbsa + year")), oy[out_flow>0 & !is.na(opp_share)], weights=~out_flow, cluster=~o_cbsa))
b0 <- bcoef(mods$base)
affplus <- rbindlist(lapply(names(mods), function(n) data.table(model=n, black_gap=round(bcoef(mods[[n]]),4),
             black_se=round(tryCatch(coeftable(mods[[n]])["grpMajority Black","Std. Error"], error=function(e) NA_real_),4),
             pct_explained=round(100*(b0-bcoef(mods[[n]]))/b0,1))), fill=TRUE)
fwrite(affplus, file.path(out_dir,"affordability_plus.csv"))
cat("\n=== AFFORDABILITY+ : Black gap and % explained by each added control ===\n"); print(affplus)

# ============================================================
# Figures
# ============================================================
oppr[, race := factor(o_race, levels=c("Majority white (NH)","No majority","Majority Hispanic","Majority Black"))]
p1 <- ggplot(oppr[!is.na(race)], aes(race, opp, fill=race)) + geom_col() +
  scale_fill_manual(values=c(`Majority Black`="#1b7837",`Majority Hispanic`="#762a83",`Majority white (NH)`="#2166ac",`No majority`="#999999"), guide="none") +
  labs(title="Environmental opportunity is smaller for minority origins",
       subtitle="Flow-weighted share of clean (below-metro-median PM2.5) neighborhoods affordable to the origin",
       x=NULL, y="Affordable-clean opportunity set (share of metro BGs)") + coord_flip() + theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"))
ggsave(file.path(fig_dir,"fig_opportunity_captured.png"), p1, width=8, height=3.8, dpi=200)

if (nrow(segq)) { sq <- copy(segq)[, `:=`(lo=black_gap-1.96*se, hi=black_gap+1.96*se)]
  p2 <- ggplot(sq, aes(seg_q, black_gap)) + geom_hline(yintercept=0, linewidth=.3) +
    geom_pointrange(aes(ymin=lo,ymax=hi), color="#1b7837") +
    labs(title="The avoidance gap is larger in more segregated metros",
         subtitle="Conditional Black avoidance gap (ug/m3) by metropolitan Black-white dissimilarity quartile",
         x="Metro segregation (dissimilarity) quartile", y="Black avoidance gap (ug/m3)") + theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"))
  ggsave(file.path(fig_dir,"fig_segregation_gap.png"), p2, width=8, height=4.2, dpi=200) }

# ============================================================
# (S17) OPPORTUNITY-SET SENSITIVITY: income multiple, rent-based, clean threshold
# ============================================================
opp_variant <- function(mult=3, clean="median", cost="own"){
  a2[, {
    thr <- if (clean=="median") median(pm25, na.rm=TRUE) else quantile(pm25, 0.25, na.rm=TRUE)
    cl  <- pm25 < thr
    if (cost=="own") { ch <- sort(med_home_value[cl]); v <- if (length(ch)) findInterval(med_hh_income*mult, ch)/.N else rep(0,.N) }
    else             { cr <- sort(med_gross_rent[cl]); v <- if (length(cr)) findInterval(med_hh_income*0.30/12, cr)/.N else rep(0,.N) }
    list(origin_bg=geoid_bg, opp_v=v)
  }, by=.(cbsa,year)]
}
variants <- list(
  `2.5x income (own)`     = list(mult=2.5, clean="median", cost="own"),
  `3x income (own, main)` = list(mult=3.0, clean="median", cost="own"),
  `4x income (own)`       = list(mult=4.0, clean="median", cost="own"),
  `5x income (own)`       = list(mult=5.0, clean="median", cost="own"),
  `Rent <=30% income`     = list(mult=NA,  clean="median", cost="rent"),
  `Clean=cleanest quartile (own 3x)` = list(mult=3.0, clean="q25", cost="own"))
opp_sens <- rbindlist(lapply(names(variants), function(nm){
  v <- variants[[nm]]; ov <- opp_variant(v$mult, v$clean, v$cost)
  d <- merge(oy[, .(origin_bg, year, o_race, out_flow)], ov, by=c("origin_bg","year"), all.x=TRUE)
  a <- d[!is.na(o_race) & !is.na(opp_v), .(opp=wmean(opp_v,out_flow)), by=o_race]
  data.table(variant=nm,
             white =a[o_race=="Majority white (NH)"]$opp,
             black =a[o_race=="Majority Black"]$opp,
             hisp  =a[o_race=="Majority Hispanic"]$opp) }), fill=TRUE)
opp_sens[, `:=`(white=round(white,3), black=round(black,3), hisp=round(hisp,3))]
opp_sens[, black_white_ratio := round(black/white,2)]
fwrite(opp_sens, file.path(out_dir,"opportunity_sensitivity.csv"))
cat("\n=== (S17) OPPORTUNITY-SET SENSITIVITY (share by origin race across definitions) ===\n"); print(opp_sens)

# ============================================================
# (S18) RACE-THRESHOLD SENSITIVITY of the avoidance gap
# ============================================================
rsh <- acs[, .(origin_bg=geoid_bg, year, sh_black, sh_hisp, sh_white)]
oy2 <- merge(oy, rsh, by=c("origin_bg","year"), all.x=TRUE)[out_flow>0]
ctrl2 <- "ns(orig_pm25,4) + log(o_inc0) + o_pov + o_bach"
# continuous racial-composition specification (coef = gap per full 0->1 share)
mcont <- feols(as.formula(paste0("re_pm25 ~ sh_black + sh_hisp + ",ctrl2," | o_cbsa + year")), oy2, weights=~out_flow, cluster=~o_cbsa)
cc <- as.data.table(coeftable(mcont), keep.rownames="term")
thr_rows <- rbindlist(lapply(c(0.5,0.6,0.7), function(t){
  d <- copy(oy2); d[, g := relevel(factor(fifelse(sh_black>t,"Black", fifelse(sh_hisp>t,"Hispanic",
                          fifelse(sh_white>t,"White","No majority")))), ref="White")]
  f <- tryCatch(feols(as.formula(paste0("re_pm25 ~ g + ",ctrl2," | o_cbsa + year")), d, weights=~out_flow, cluster=~o_cbsa), error=function(e) NULL)
  if (is.null(f)) return(NULL)
  ct <- as.data.table(coeftable(f), keep.rownames="term")
  data.table(spec=paste0("majority > ",t*100,"%"),
             black_gap=ct[term=="gBlack"]$Estimate, black_se=ct[term=="gBlack"]$`Std. Error`,
             hisp_gap=ct[term=="gHispanic"]$Estimate, hisp_se=ct[term=="gHispanic"]$`Std. Error`) }), fill=TRUE)
thr_rows <- rbind(thr_rows, data.table(spec="continuous (per 0->1 share)",
             black_gap=cc[term=="sh_black"]$Estimate, black_se=cc[term=="sh_black"]$`Std. Error`,
             hisp_gap=cc[term=="sh_hisp"]$Estimate, hisp_se=cc[term=="sh_hisp"]$`Std. Error`))
thr_rows[, `:=`(black_gap=round(black_gap,4), black_se=round(black_se,4), hisp_gap=round(hisp_gap,4), hisp_se=round(hisp_se,4))]
fwrite(thr_rows, file.path(out_dir,"race_threshold_sensitivity.csv"))
cat("\n=== (S18) RACE-THRESHOLD SENSITIVITY of the avoidance gap ===\n"); print(thr_rows)

cat("\nDONE. Opportunity + heterogeneity + affordability-plus outputs in ", out_dir, "\n", sep="")
