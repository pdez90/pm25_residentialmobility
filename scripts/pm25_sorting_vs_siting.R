# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# SORTING vs SITING: why is PM2.5 exposure racially disparate?
# ------------------------------------------------------------
# The classic environmental-justice question (Been & Gupta 1997,
# "Coming to the Nuisance or Going to the Barrios?"; Depro, Timmins &
# O'Neill 2015; Hernandez-Cortes & Lopez-Feldman 2022): is the
# demographic-pollution correlation maintained because (1) minorities
# selectively move INTO high-pollution neighborhoods / fail to leave
# them (SORTING), or because (2) pollution rises or persists where
# minorities already live (SITING / disparate pollution dynamics)?
#
# Prior work infers this from net demographic change and pollution
# change. Here we observe the actual migration FLOWS, so we can test
# both channels directly:
#   Section 1  SITING:  does PM2.5 decline less / rise where minorities
#              live? BG-level dPM2.5 (2011->2019) on baseline minority
#              share, net of baseline level, county FE.
#   Section 2  SORTING: do higher-minority migrant streams flow into
#              dirtier destinations, and do dirty places attract
#              disproportionately minority in-migration relative to
#              their own residents?  (The mover-level "fail to escape"
#              side is the conditional avoidance privilege, oy.)
#   Section 3  Which dominates: the shift-share (sorting = demographics
#              via migration; cleanup = pollution dynamics) from v3.
#
# Reuses poll/acs/oy from idea2_pm25_avoidance.R.
# Output: $MIGRATE_DATA_DIR/pm25_avoidance_lag3_olag1/
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr); library(fixest); library(splines); library(ggplot2) })
base_dir <- migrate_data_dir(); mig_dir <- base_dir
tag <- "lag3_olag1"; migrate_years <- 2011:2019
out_dir <- file.path(base_dir, paste0("pm25_avoidance_", tag)); fig_dir <- file.path(out_dir,"figures")
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")
wmean  <- function(x,w){ok<-!is.na(x)&!is.na(w); if(!any(ok)) return(NA_real_); sum(x[ok]*w[ok])/sum(w[ok])}

if (!exists("acs") || !exists("poll")) source(file.path(base_dir,"scripts","idea2_pm25_avoidance.R"))
# attach race shares to acs if absent
if (!all(c("sh_black","sh_hisp") %in% names(acs))) {
  rcs <- list.files(file.path(base_dir,"pollution_bg"), "^race_bg_.*csv$", full.names=TRUE)
  rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}"))
        r<-fread(f,colClasses=c(geoid_bg="character")); r[,year:=y]; r }), fill=TRUE)
  rc[, geoid_bg := pad_bg(geoid_bg)]
  acs <- merge(acs, rc[, .(year,geoid_bg,sh_black,sh_hisp,sh_white)], by=c("year","geoid_bg"), all.x=TRUE)
}
acs[, sh_min := sh_black + sh_hisp]

# ============================================================
# 1. SITING: pollution dynamics on baseline minority share
# ============================================================
bg0 <- acs[year==min(migrate_years), .(geoid_bg, sh_black0=sh_black, sh_hisp0=sh_hisp, sh_min0=sh_min,
                                        pm0=pm25, pop0=total_pop, county=substr(geoid_bg,1,5))]
bg1 <- acs[year==max(migrate_years), .(geoid_bg, pm1=pm25)]
bg  <- merge(bg0, bg1, by="geoid_bg")[!is.na(pm0)&!is.na(pm1)&!is.na(sh_min0)]
bg[, dpm := pm1 - pm0]
site_raw <- feols(dpm ~ sh_black0 + sh_hisp0 | county, data=bg, weights=~pop0, cluster=~county)
site_adj <- feols(dpm ~ sh_black0 + sh_hisp0 + ns(pm0,4) | county, data=bg, weights=~pop0, cluster=~county)
siting <- rbindlist(list(
  cbind(spec="raw (county FE)",       as.data.table(coeftable(site_raw), keep.rownames="term")),
  cbind(spec="net of baseline PM2.5", as.data.table(coeftable(site_adj), keep.rownames="term"))), fill=TRUE)
siting <- siting[grepl("sh_", term)]
fwrite(siting, file.path(out_dir,"siting_pollution_dynamics.csv"))
cat("\n=== SITING: dPM2.5 (2011->2019) on baseline minority share ",
    "(+ = pollution declined LESS / rose where minorities live) ===\n"); print(siting)

# ---- SITING ROBUSTNESS: vary baseline-PM2.5 control + fixed effects ----
xw_file <- file.path(base_dir, "pollution_bg", "county_cbsa_xwalk.csv")
if (file.exists(xw_file)) {
  cxw <- fread(xw_file, colClasses = c(county="character", cbsa="character"))
  bg <- merge(bg, cxw, by = "county", all.x = TRUE); bg[is.na(cbsa), cbsa := paste0("nonmetro_", county)]
}
site_specs <- list(
  list(name="(1) No baseline control, county FE",   rhs="sh_black0 + sh_hisp0",              fe="county"),
  list(name="(2) Linear baseline PM2.5, county FE",  rhs="sh_black0 + sh_hisp0 + pm0",        fe="county"),
  list(name="(3) Spline df=2, county FE",            rhs="sh_black0 + sh_hisp0 + ns(pm0,2)",  fe="county"),
  list(name="(4) Spline df=4, county FE (main)",     rhs="sh_black0 + sh_hisp0 + ns(pm0,4)",  fe="county"),
  list(name="(5) Spline df=6, county FE",            rhs="sh_black0 + sh_hisp0 + ns(pm0,6)",  fe="county"),
  list(name="(6) Spline df=4, CBSA FE",              rhs="sh_black0 + sh_hisp0 + ns(pm0,4)",  fe="cbsa"),
  list(name="(7) Spline df=4, national (no FE)",     rhs="sh_black0 + sh_hisp0 + ns(pm0,4)",  fe="none"))
run_site <- function(sp) {
  if (sp$fe == "cbsa" && !"cbsa" %in% names(bg)) return(NULL)
  clu <- switch(sp$fe, county=~county, cbsa=~cbsa, none=~county)
  f <- if (sp$fe == "none") as.formula(paste0("dpm ~ ", sp$rhs))
       else as.formula(paste0("dpm ~ ", sp$rhs, " | ", sp$fe))
  fit <- tryCatch(feols(f, data=bg, weights=~pop0, cluster=clu), error=function(e) NULL); if (is.null(fit)) return(NULL)
  ct <- as.data.table(coeftable(fit), keep.rownames="term")[term %in% c("sh_black0","sh_hisp0")]
  ct[, spec := sp$name][]
}
site_rob <- rbindlist(lapply(site_specs, run_site), fill=TRUE)
fwrite(site_rob, file.path(out_dir, "siting_robustness.csv"))
cat("\n=== SITING ROBUSTNESS: coef on baseline minority share across specs (+ = declined less where minorities lived) ===\n")
print(site_rob[, .(spec, term, Estimate=round(Estimate,3), p=signif(`Pr(>|t|)`,2))])
sr <- copy(site_rob)[, `:=`(lo=Estimate-1.96*`Std. Error`, hi=Estimate+1.96*`Std. Error`,
                            grp=ifelse(grepl("black",term),"Black share","Hispanic share"))]
sr[, spec := factor(spec, levels=rev(vapply(site_specs, `[[`, character(1), "name")))]
p9 <- ggplot(sr, aes(Estimate, spec, color=grp)) +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_point(position=position_dodge(.5), size=2.4) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=.2, position=position_dodge(.5)) +
  scale_color_manual(values=c(`Black share`="#1b7837", `Hispanic share`="#762a83")) +
  labs(title="Unequal cleanup robustness: pollution persistence in minority neighborhoods",
       subtitle="Coefficient on baseline minority share in dPM2.5 (2011-2019), across baseline-PM2.5 controls and fixed effects",
       x=expression("Coefficient ("*mu*g/m^3*"); + = PM2.5 declined less where minorities lived"), y=NULL, color=NULL) +
  theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"))
ggsave(file.path(fig_dir, "fig9_siting_robustness.png"), p9, width=9, height=4.8, dpi=200)

# ============================================================
# 2. SORTING: destination pull of minority in-migration
# ============================================================
oa_min <- acs[, .(cy=year, origin_bg=geoid_bg, o_min=sh_min)]; setkey(oa_min, cy, origin_bg)
dd     <- acs[, .(cy=year, dest_bg=geoid_bg, d_pm25=pm25, d_min=sh_min, d_county=substr(geoid_bg,1,5))]; setkey(dd, cy, dest_bg)
acc <- list()
for (yr in migrate_years){
  x <- fread(file.path(mig_dir, paste0(yr,".csv")), select=c("Origin","Destination","Flow"))
  x <- x[, .(origin_bg=pad_bg(Origin), dest_bg=pad_bg(Destination), flow=as.numeric(Flow))][flow>0 & origin_bg!=dest_bg]
  x[, cy:=yr]; setkey(x,cy,origin_bg); x <- oa_min[x, on=.(cy,origin_bg)]
  # aggregate to destination-year: flow-weighted origin minority share of IN-migrants
  acc[[as.character(yr)]] <- x[, .(in_flow=sum(flow,na.rm=TRUE), in_min=wmean(o_min,flow)), by=.(cy,dest_bg)]
}
inmig <- rbindlist(acc, fill=TRUE); setkey(inmig,cy,dest_bg); inmig <- dd[inmig, on=.(cy,dest_bg)]
inmig <- inmig[!is.na(d_pm25) & !is.na(in_min)]
# do dirtier destinations pull more-minority in-migrant streams?
sort_pull <- feols(d_pm25 ~ in_min | d_county + cy, data=inmig, weights=~in_flow, cluster=~d_county)
# do dirty places attract MORE-minority inflows than their own residents? (in_min - d_min)
inmig[, minority_pull := in_min - d_min]
sort_excess <- feols(minority_pull ~ ns(d_pm25,4) | d_county + cy, data=inmig, weights=~in_flow, cluster=~d_county)
cat("\n=== SORTING: destination PM2.5 ~ minority share of in-migrant streams (+ = minorities move to dirtier) ===\n")
print(as.data.table(coeftable(sort_pull), keep.rownames="term")[grepl("in_min",term)])
fwrite(as.data.table(coeftable(sort_pull), keep.rownames="term"), file.path(out_dir,"sorting_destination_pull.csv"))
# binned scatter: dest PM2.5 quintile vs in-migrant minority share and resident minority share
inmig[, pm_q := cut(d_pm25, quantile(d_pm25, 0:5/5, na.rm=TRUE), include.lowest=TRUE, labels=paste0("Q",1:5))]
binned <- inmig[!is.na(pm_q), .(in_migrant_minority = wmean(in_min,in_flow),
                                resident_minority = wmean(d_min,in_flow),
                                mean_pm25 = wmean(d_pm25,in_flow)), by=pm_q][order(pm_q)]
fwrite(binned, file.path(out_dir,"sorting_by_pm25_quintile.csv"))
cat("\n=== In-migrant vs resident minority share, by destination PM2.5 quintile ===\n"); print(binned)

# ============================================================
# 3. Figures
# ============================================================
# siting coefficients
sp <- copy(siting); sp[, `:=`(lo=Estimate-1.96*`Std. Error`, hi=Estimate+1.96*`Std. Error`,
                              grp=ifelse(grepl("black",term),"Black share","Hispanic share"))]
p1 <- ggplot(sp, aes(Estimate, grp, color=spec)) +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_point(position=position_dodge(.4), size=2.6) +
  geom_errorbarh(aes(xmin=lo,xmax=hi), height=.2, position=position_dodge(.4)) +
  scale_color_manual(values=c("#999999","#b2182b")) +
  labs(title="Unequal cleanup: did PM2.5 decline less where minorities lived?",
       subtitle="BG change in PM2.5 2011-2019 on baseline minority share; county FE, pop-weighted. + = pollution stayed/came.",
       x=expression(Coefficient~on~baseline~minority~share~"("~mu*g/m^3~")"), y=NULL, color=NULL) +
  theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"))
ggsave(file.path(fig_dir,"fig7_siting_pollution_dynamics.png"), p1, width=8, height=3.8, dpi=200)

bl <- melt(binned[, .(pm_q, `In-migrants`=in_migrant_minority, Residents=resident_minority)],
           id.vars="pm_q", variable.name="group", value.name="minority_share")
p2 <- ggplot(bl, aes(pm_q, minority_share, fill=group)) + geom_col(position="dodge") +
  scale_fill_manual(values=c(`In-migrants`="#d6604d", Residents="#4d4d4d")) +
  labs(title="SORTING: dirtier destinations draw more-minority in-migration",
       subtitle="Flow-weighted minority share of in-migrant streams vs resident minority share, by destination PM2.5 quintile",
       x="Destination PM2.5 quintile (Q5 = dirtiest)", y="Minority population share", fill=NULL) +
  theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"))
ggsave(file.path(fig_dir,"fig8_sorting_by_pm25_quintile.png"), p2, width=8, height=4.2, dpi=200)

cat("\nDONE. Siting + sorting outputs and figures in ", out_dir, "\n", sep="")
cat("Read alongside the v3 shift-share: sorting = population redistribution across neighborhoods (migration, plus births, deaths, and other demographic change);\n",
    "cleanup = pollution dynamics. Together they quantify which channel maintains the disparity.\n", sep="")
