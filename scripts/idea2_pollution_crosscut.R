# ============================================================
# CROSS-POLLUTANT: does the same mechanism hold for NO2 and
# wildfire PM2.5 as for total PM2.5?
# ------------------------------------------------------------
# For each pollutant in the panel (pm25, no2, wildfire_pm25) it
# re-estimates the three headline quantities and assembles a single
# cross-pollutant comparison, so the paper can show the mechanism is
# about environmental inequality broadly, not one pollutant:
#   (i)   CONDITIONAL avoidance gap by origin race (CBSA + year FE,
#         origin-pollutant spline) - the avoidance privilege;
#   (ii)  SORTING - destination pollutant on minority share of
#         in-migrant streams;
#   (iii) SITING - block-group change in the pollutant on baseline
#         minority share, net of baseline level.
#
# Self-contained: builds its inputs from the multi-pollutant panel,
# the gentrification file, the cached race pulls, and the CBSA
# crosswalk (from idea2_pm25_avoidance_v3.R). One MIGRATE pass per
# year serves all pollutants.
#
# Output: /Users/priyanka/MIGRATE_dir/pm25_avoidance_lag3_olag1/
#   cross_pollutant_estimates.csv, figures/fig10_cross_pollutant.png
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr); library(fixest); library(splines); library(ggplot2) })
base_dir <- "/Users/priyanka/MIGRATE_dir"; mig_dir <- base_dir
gent_dir <- file.path(base_dir,"acs_gentrification"); change_lag <- 3; origin_lag <- 1
tag <- "lag3_olag1"; migrate_years <- 2011:2019
out_dir <- file.path(base_dir, paste0("pm25_avoidance_",tag)); fig_dir <- file.path(out_dir,"figures")
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
panel_file <- file.path(base_dir,"pollution_bg","bg_pollution_panel.csv")
combined_file <- file.path(gent_dir, paste0("blockgroup_gentrification_2009_2019_lag",change_lag,"_with_HMDA.csv"))
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")
wmean  <- function(x,w){ok<-!is.na(x)&!is.na(w); if(!any(ok)) return(NA_real_); sum(x[ok]*w[ok])/sum(w[ok])}
PLAB <- c(pm25="Total PM2.5", no2="NO2", wildfire_pm25="Wildfire PM2.5", ozone="Ozone (warm-season MDA8)",
          so4="Sulfate (power)", no3="Nitrate (traffic/ag)", bc="Black carbon (diesel)", om="Organic matter (combustion)")

# ---- pollution panel (multi) + speciated components ----
stopifnot("multi-pollutant panel not found; run extract_pollution_multi_bg.R" = file.exists(panel_file))
poll <- fread(panel_file, colClasses=c(geoid_bg="character")); poll[, geoid_bg := pad_bg(geoid_bg)]
COMPONENTS_USE <- c("so4","no3","bc","om")   # combustion/industry/traffic components to analyze
comp_file <- file.path(base_dir,"pollution_bg","bg_pm25_components_panel.csv")
if (file.exists(comp_file)) {
  cp <- fread(comp_file, colClasses=c(geoid_bg="character")); cp[, geoid_bg := pad_bg(geoid_bg)]
  have_comp <- intersect(COMPONENTS_USE, names(cp))
  stopifnot("components requested but columns missing in components panel" = length(have_comp)==length(COMPONENTS_USE))
  poll <- merge(poll, cp[, c("geoid_bg","year",have_comp), with=FALSE], by=c("geoid_bg","year"), all.x=TRUE)
} else message("Component panel not found; skipping components (set RUN_COMPONENTS <- TRUE in run_all_pm25.R).")
# guard against duplicate (geoid_bg, year) rows (e.g. wildfire centroid-in-grid multi-matches)
ndup <- nrow(poll) - nrow(unique(poll, by=c("geoid_bg","year")))
if (ndup > 0) { message("Deduplicating pollution panel: ", ndup, " duplicate BG-year rows removed."); poll <- unique(poll, by=c("geoid_bg","year")) }
polls <- intersect(c("pm25","no2","wildfire_pm25","ozone", COMPONENTS_USE), names(poll))
polls <- polls[sapply(polls, function(p) mean(!is.na(poll[[p]]))>0.2)]   # keep pollutants with coverage
message("Pollutants analyzed: ", paste(polls, collapse=", "))
# national SD of each pollutant (across BG-years) -> standardize for cross-pollutant comparison
nsd <- sapply(polls, function(p) sd(poll[[p]], na.rm=TRUE))

# ---- composition + gentrification ----
hdr <- names(fread(combined_file, nrows=0))
cap_col <- if ("gentrifying_capital_top10_screened" %in% hdr) "gentrifying_capital_top10_screened" else "gentrifying_capital_top10"
acs <- fread(combined_file, colClasses=c(geoid_bg="character", state_county="character"),
  select=unique(c("year","geoid_bg","state_county","eligible_gentrifiable","gentrifying_transition_top10",
                  "med_hh_income","poverty_pct","bach_pct","total_pop")))
acs[, geoid_bg := pad_bg(geoid_bg)]
acs[, inc_tercile := cut(med_hh_income, quantile(med_hh_income,c(0,1/3,2/3,1),na.rm=TRUE), labels=c("Low income","Middle","High income"), include.lowest=TRUE), by=year]
acs[, gent_group := fifelse(gentrifying_transition_top10==1,"Gentrifying", fifelse(gentrifying_transition_top10==0 & eligible_gentrifiable==TRUE,"Eligible, not gentrifying",NA_character_))]
# race
rcs <- list.files(file.path(base_dir,"pollution_bg"), "^race_bg_.*csv$", full.names=TRUE)
rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}")); r<-fread(f,colClasses=c(geoid_bg="character")); r[,year:=y]; r }), fill=TRUE)
rc[, geoid_bg := pad_bg(geoid_bg)]
acs <- merge(acs, rc[, .(year,geoid_bg,sh_black,sh_hisp,sh_white)], by=c("year","geoid_bg"), all.x=TRUE)
acs[, `:=`(sh_min=sh_black+sh_hisp,
           race_majority=fifelse(sh_black>0.5,"Majority Black", fifelse(sh_hisp>0.5,"Majority Hispanic", fifelse(sh_white>0.5,"Majority white (NH)","No majority"))))]
acs <- merge(acs, poll, by=c("geoid_bg","year"), all.x=TRUE)
acs[, county := substr(geoid_bg,1,5)]
xw_file <- file.path(base_dir,"pollution_bg","county_cbsa_xwalk.csv")
if (file.exists(xw_file)) { cxw<-fread(xw_file, colClasses=c(county="character",cbsa="character")); acs<-merge(acs,unique(cxw,by="county"),by="county",all.x=TRUE); acs[is.na(cbsa),cbsa:=paste0("nonmetro_",county)] } else acs[, cbsa := county]
acs <- unique(acs, by=c("year","geoid_bg"))   # ensure one row per BG-year before keyed joins

# ---- origin / dest lookups (all pollutants attached) ----
# built from acs preserving row order, so acs[[p]] aligns column-wise
oa <- acs[, .(cy=year, origin_bg=geoid_bg, o_race=race_majority, o_inc=inc_tercile, o_gent=gent_group,
              o_inc0=med_hh_income, o_pov=poverty_pct, o_bach=bach_pct, o_min=sh_min, o_county=county, o_cbsa=cbsa)]
for (p in polls) oa[, (paste0("o_",p)) := acs[[p]]]
setkey(oa, cy, origin_bg)
da <- acs[, .(cy=year, dest_bg=geoid_bg, d_county=county, d_min=sh_min)]
for (p in polls) da[, (paste0("d_",p)) := acs[[p]]]
setkey(da, cy, dest_bg)

# ---- one MIGRATE pass: origin-agg (conditional) + dest-agg (sorting) ----
acc_o <- list(); acc_d <- list()
for (yr in migrate_years) {
  x <- fread(file.path(mig_dir, paste0(yr,".csv")), select=c("Origin","Destination","Flow"))
  x <- x[, .(origin_bg=pad_bg(Origin), dest_bg=pad_bg(Destination), flow=as.numeric(Flow))][flow>0 & origin_bg!=dest_bg]
  x[, cy:=yr]; setkey(x,cy,origin_bg); x <- oa[x, on=.(cy,origin_bg)]; setkey(x,cy,dest_bg); x <- da[x, on=.(cy,dest_bg)]
  for (p in polls) x[, (paste0("re_",p)) := get(paste0("d_",p)) - get(paste0("o_",p))]
  oexpr <- c(list(year=quote(yr), out_flow=quote(sum(flow,na.rm=TRUE))),
             setNames(lapply(polls, function(p) bquote(wmean(.(as.name(paste0("re_",p))),flow))), paste0("re_",polls)),
             setNames(lapply(polls, function(p) bquote(wmean(.(as.name(paste0("o_",p))),flow))), paste0("orig_",polls)))
  acc_o[[as.character(yr)]] <- x[, eval(as.call(c(quote(list),oexpr))),
        by=.(origin_bg, o_race, o_inc, o_gent, o_inc0, o_pov, o_bach, o_county, o_cbsa)]
  dexpr <- c(list(cy=quote(yr), in_flow=quote(sum(flow,na.rm=TRUE)), in_min=quote(wmean(o_min,flow))),
             setNames(lapply(polls, function(p) bquote(data.table::first(.(as.name(paste0("d_",p)))))), paste0("d_",polls)),
             list(d_county=quote(data.table::first(d_county))))
  acc_d[[as.character(yr)]] <- x[, eval(as.call(c(quote(list),dexpr))), by=dest_bg]
}
oy <- rbindlist(acc_o, fill=TRUE); inmig <- rbindlist(acc_d, fill=TRUE)

# ============================================================
# Estimate the three mechanisms per pollutant
# ============================================================
res <- list()
for (p in polls) {
  # (i) conditional avoidance gap by race, CBSA FE, origin-pollutant spline
  d <- oy[!is.na(o_race) & out_flow>0 & !is.na(o_inc0) & o_inc0>0 & !is.na(get(paste0("orig_",p)))]
  d[, grp := relevel(factor(o_race), ref="Majority white (NH)")]
  f <- as.formula(paste0("re_",p," ~ grp + ns(orig_",p,",4) + log(o_inc0) + o_pov + o_bach | o_cbsa + year"))
  fit <- tryCatch(feols(f, data=d, weights=~out_flow, cluster=~o_cbsa), error=function(e) NULL)
  if (!is.null(fit)) { ct<-as.data.table(coeftable(fit),keep.rownames="term")[grepl("^grp",term)]
    res[[paste("cond",p)]] <- ct[, .(mechanism="Conditional avoidance gap", pollutant=p, group=sub("^grpMajority ","",term),
                                     est=Estimate, se=`Std. Error`, p=`Pr(>|t|)`)] }
  # (ii) sorting: destination pollutant on in-migrant minority share
  im <- inmig[!is.na(get(paste0("d_",p))) & !is.na(in_min)]
  fs <- as.formula(paste0("d_",p," ~ in_min | d_county + cy"))
  fit2 <- tryCatch(feols(fs, data=im, weights=~in_flow, cluster=~d_county), error=function(e) NULL)
  if (!is.null(fit2)) { ct<-as.data.table(coeftable(fit2),keep.rownames="term")[term=="in_min"]
    res[[paste("sort",p)]] <- ct[, .(mechanism="Sorting (minority in-migration)", pollutant=p, group="Minority in-migrant share",
                                     est=Estimate, se=`Std. Error`, p=`Pr(>|t|)`)] }
  # (iii) siting: BG change in pollutant on baseline minority share, net of baseline
  b0 <- acs[year==min(migrate_years), .(geoid_bg, sh_black0=sh_black, sh_hisp0=sh_hisp, p0=get(p), pop0=total_pop, county)]
  b1 <- acs[year==max(migrate_years), .(geoid_bg, p1=get(p))]
  bb <- merge(b0,b1,by="geoid_bg")[!is.na(p0)&!is.na(p1)&!is.na(sh_black0)]; bb[, dp:=p1-p0]
  fit3 <- tryCatch(feols(dp ~ sh_black0 + sh_hisp0 + ns(p0,4) | county, data=bb, weights=~pop0, cluster=~county), error=function(e) NULL)
  if (!is.null(fit3)) { ct<-as.data.table(coeftable(fit3),keep.rownames="term")[term %in% c("sh_black0","sh_hisp0")]
    res[[paste("site",p)]] <- ct[, .(mechanism="Unequal cleanup (persistence)", pollutant=p,
                                     group=ifelse(term=="sh_black0","Black share","Hispanic share"),
                                     est=Estimate, se=`Std. Error`, p=`Pr(>|t|)`)] }
}
cross <- rbindlist(res, fill=TRUE)
cross[, pollutant_lab := factor(PLAB[pollutant], levels=unname(PLAB[polls]))]
# standardized estimate (national SDs) so pollutants in different units are comparable
cross[, `:=`(est_sd = est/nsd[pollutant], se_sd = se/nsd[pollutant])]
fwrite(cross, file.path(out_dir, "cross_pollutant_estimates.csv"))
cat("\n=== Cross-pollutant estimates (three mechanisms; native units + national-SD) ===\n")
print(cross[, .(mechanism, pollutant, group, est=round(est,3), est_sd=round(est_sd,3), p=signif(p,2))])

# ============================================================
# Figure: same mechanism across pollutants (conditional gap, national-SD units)
# ============================================================
cc <- cross[mechanism=="Conditional avoidance gap" & group %in% c("Black","Hispanic")]
cc[, `:=`(lo=est_sd-1.96*se_sd, hi=est_sd+1.96*se_sd)]
p10 <- ggplot(cc, aes(est_sd, pollutant_lab, color=group)) +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_point(position=position_dodge(.5), size=2.8) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=.2, position=position_dodge(.5)) +
  scale_color_manual(values=c(Black="#1b7837", Hispanic="#762a83")) +
  labs(title="The avoidance privilege holds across pollutants",
       subtitle="Conditional gap: destination pollutant reached by minority- vs white-origin movers, net of origin level (CBSA + year FE)",
       x="Adjusted gap vs white-origin movers (national SDs; + = reached dirtier air)", y=NULL, color=NULL) +
  theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"))
ggsave(file.path(fig_dir, "fig10_cross_pollutant.png"), p10, width=8.5, height=4.2, dpi=200)

# ---- 3-mechanism signature figure (all pollutants x avoidance/sorting/siting) ----
c3 <- copy(cross)
c3[, grp2 := fifelse(grepl("Black",group),"Black", fifelse(grepl("Hispanic",group),"Hispanic","Minority in-migration"))]
c3[, `:=`(lo=est_sd-1.96*se_sd, hi=est_sd+1.96*se_sd)]
ord <- intersect(c("no2","bc","om","pm25","no3","so4","ozone","wildfire_pm25"), polls)
c3[, plab := factor(PLAB[pollutant], levels=rev(unname(PLAB[ord])))]
c3[, mech := factor(mechanism, levels=c("Conditional avoidance gap","Sorting (minority in-migration)","Unequal cleanup (persistence)"),
                    labels=c("A. Avoidance privilege","B. Sorting","C. Unequal cleanup"))]
pmech <- ggplot(c3[!is.na(plab)], aes(est_sd, plab, color=grp2)) +
  geom_vline(xintercept=0, linetype=2, color="grey60") +
  geom_pointrange(aes(xmin=lo, xmax=hi), position=position_dodge(.5), size=.35) +
  facet_wrap(~mech, scales="free_x") +
  scale_color_manual(values=c(Black="#1b7837", Hispanic="#762a83", `Minority in-migration`="#b2182b"), name=NULL) +
  labs(title="Environmental-inequality mobility mechanisms concentrate in traffic/combustion pollution",
       subtitle="Effect in national SDs; + = minority-linked movement toward dirtier air / slower cleanup",
       x="Effect (national SDs)", y=NULL) +
  theme_minimal(base_size=10) + theme(plot.title=element_text(face="bold"), legend.position="bottom")
ggsave(file.path(fig_dir, "fig10_cross_pollutant_3panel.png"), pmech, width=13, height=4.8, dpi=200)
cat("\nDONE. Cross-pollutant outputs in ", out_dir, "\n", sep="")
