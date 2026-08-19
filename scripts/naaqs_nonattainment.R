# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# NAAQS NONATTAINMENT as an EXOGENOUS CLEANUP SHOCK
# ------------------------------------------------------------
# The Clean Air Act imposes stricter regulation on counties whose
# pollution crosses a federal threshold ("nonattainment"), producing
# quasi-random differential cleanup (Chay & Greenstone 2005; Sullivan
# 2017). We use that shock three ways:
#   (A) FIRST STAGE / VALIDITY: an event study showing nonattainment
#       counties clean faster after designation (the instrument moves
#       pollution).
#   (B) CAPTURE OF CLEANUP (origin streams): among movers LEAVING
#       counties subject to regulation-driven cleanup, do minority-origin
#       streams capture LESS of the improvement than white-origin streams?
#       (reduced-form interaction, and an IV where origin nonattainment
#       instruments the origin's realized pollution decline.)
#   (C) CAPTURE OF CLEANUP (destination in-migration): do regulation-
#       cleaned counties draw relatively FEWER minority in-migrants, i.e.
#       is the exogenous improvement captured disproportionately by
#       white in-migration?
#
# WHAT TO DOWNLOAD:
#   EPA Green Book nonattainment area listings (epa.gov/green-book).
#   Prepare a county-level CSV -> MIGRATE_dir/naaqs/greenbook_nonattainment.csv
#   with columns:
#       county        5-digit county FIPS
#       pollutant     e.g. "PM2.5" or "Ozone"
#       design_year   year the county was designated nonattainment
#   (One row per county x pollutant designation.)
#
# Output: pm25_avoidance_lag3_olag1/
#   naaqs_firststage_eventstudy.csv + figures/fig_naaqs_firststage.png
#   naaqs_capture_origin.csv        + figures/fig_naaqs_capture.png
#   naaqs_capture_destination.csv
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr); library(fixest); library(splines); library(ggplot2) })
base_dir <- migrate_data_dir(); mig_dir <- base_dir
out_dir  <- file.path(base_dir,"pm25_avoidance_lag3_olag1"); fig_dir <- file.path(out_dir,"figures")
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
migrate_years <- 2011:2019; POLL <- "PM2.5"
pad_bg  <- function(x) str_pad(as.character(x),12,pad="0")
pad_cty <- function(x) str_pad(as.character(x),5,pad="0")

# ---- Green Book nonattainment (hard error if missing) ----
gbf <- file.path(base_dir,"naaqs","greenbook_nonattainment.csv")
stopifnot("Green Book file not found: naaqs/greenbook_nonattainment.csv (see header)" = file.exists(gbf))
gb <- fread(gbf); setnames(gb, tolower(names(gb)))
stopifnot("Green Book needs columns county, pollutant, design_year" = all(c("county","pollutant","design_year") %in% names(gb)))
gb[, county:=pad_cty(county)]
gb <- gb[grepl(gsub("\\.","", tolower(POLL)), gsub("\\.","", tolower(pollutant)))]      # PM2.5 rows
stopifnot("no nonattainment rows for chosen pollutant" = nrow(gb)>0)
naa <- unique(gb[, .(design_year=min(as.integer(design_year), na.rm=TRUE)), by=county])  # first designation
naa[, treat:=1L]

# ---- county-year PM2.5 (population-weighted from BG panel) ----
poll <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character"))
poll[, geoid_bg:=pad_bg(geoid_bg)]; poll <- unique(poll[,.(geoid_bg,year,pm25)], by=c("geoid_bg","year"))
pop <- fread(file.path(base_dir,"acs_gentrification","blockgroup_gentrification_2009_2019_lag3_with_HMDA.csv"),
             colClasses=c(geoid_bg="character"), select=c("year","geoid_bg","total_pop")); pop[, geoid_bg:=pad_bg(geoid_bg)]
poll <- merge(poll, pop, by=c("geoid_bg","year"), all.x=TRUE); poll[is.na(total_pop), total_pop:=1]
poll[, county:=substr(geoid_bg,1,5)]
cty <- poll[!is.na(pm25), .(pm25=sum(pm25*total_pop)/sum(total_pop)), by=.(county,year)]
cty <- merge(cty, naa, by="county", all.x=TRUE); cty[is.na(treat), treat:=0L]

# ============================================================
# (A) FIRST STAGE: event study of county PM2.5 around designation
# ============================================================
# DiD event study using NEVER-designated counties as controls. A treated-
# ONLY event study is collinear with the year FE when all treated share one
# designation year (the 2012 PM2.5 standard designates every county in 2015),
# so we identify the effect against untreated counties instead.
cty[, treat := as.integer(!is.na(design_year))]
dy  <- sort(unique(cty[treat==1 & !is.na(design_year)]$design_year))
REF <- (if (length(dy)==1) dy else round(median(dy))) - 1          # last pre-designation year
esd <- cty[year >= min(dy)-4 & year <= max(dy)+4]
m_es <- tryCatch(feols(pm25 ~ i(year, treat, ref=REF) | county + year, esd, cluster=~county), error=function(e) NULL)
if (!is.null(m_es)) {
  et <- as.data.table(coeftable(m_es), keep.rownames="term")[grepl("year::",term)]
  et[, yr := as.integer(sub(".*year::(-?[0-9]+).*","\\1",term))]
  fwrite(et, file.path(out_dir,"naaqs_firststage_eventstudy.csv"))
  cat("\n=== FIRST STAGE (DiD): treated vs never-designated county PM2.5, ref=", REF, " ===\n", sep="")
  print(et[order(yr), .(year=yr, est=round(Estimate,3), se=round(`Std. Error`,3), p=signif(`Pr(>|t|)`,2))])
  eb <- rbind(et[, .(yr, est=Estimate, se=`Std. Error`)], data.table(yr=REF, est=0, se=0))
  gA <- ggplot(eb[order(yr)], aes(yr, est)) + geom_hline(yintercept=0, linewidth=.3) + geom_vline(xintercept=REF+.5, linetype=2, color="grey50") +
    geom_pointrange(aes(ymin=est-1.96*se, ymax=est+1.96*se), color="#2166ac") +
    labs(title="Nonattainment designation is followed by faster PM2.5 decline",
         subtitle=paste0("DiD event study: treated (nonattainment) vs never-designated counties, reference year ", REF),
         x="Year", y=expression("PM2.5 vs reference year ("*mu*g/m^3*")")) + theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"))
  ggsave(file.path(fig_dir,"fig_naaqs_firststage.png"), gA, width=8, height=4.2, dpi=200)
} else cat("First-stage event study not estimable.\n")

# ============================================================
# (B) CAPTURE OF CLEANUP at the ORIGIN (stream-level)
# ============================================================
oyf <- file.path(out_dir,"origin_year_affordability.csv")
stopifnot("run affordability_robustness.R (stage 7) first" = file.exists(oyf))
oy <- fread(oyf, colClasses=c(origin_bg="character")); oy[, origin_bg:=pad_bg(origin_bg)]
oy[, o_county:=pad_cty(o_county)]
oy <- merge(oy, naa[,.(o_county=county, treat_o=treat, o_design=design_year)], by="o_county", all.x=TRUE)
oy[is.na(treat_o), treat_o:=0L]
oy[, grp := relevel(factor(o_race), ref="Majority white (NH)")]
ctrl <- "ns(orig_pm25,4) + log(o_inc0) + o_pov + o_bach"
# reduced form: in regulation-cleaned (nonattainment) origins, do minority streams capture less?
m_cap <- feols(as.formula(paste0("re_pm25 ~ treat_o*grp + ",ctrl," | o_cbsa + year")),
               oy[out_flow>0], weights=~out_flow, cluster=~o_cbsa)
captab <- as.data.table(coeftable(m_cap), keep.rownames="term")[grepl("treat_o|grp",term)]
fwrite(captab, file.path(out_dir,"naaqs_capture_origin.csv"))
cat("\n=== CAPTURE OF CLEANUP (origin): reduced-form treat x race on exposure reduction ===\n")
print(captab[, .(term, est=round(Estimate,3), se=round(`Std. Error`,3), p=signif(`Pr(>|t|)`,2))])
# IV: origin nonattainment instruments the origin's realized pollution decline; does captured cleanup differ by race
oy[, o_change := orig_pm25]   # origin baseline pollution; instrumented cleanup enters via treat_o (post-designation) x baseline
m_iv <- tryCatch(feols(as.formula(paste0("re_pm25 ~ grp + log(o_inc0) + o_pov + o_bach | o_cbsa + year | orig_pm25 ~ treat_o")),
                       oy[out_flow>0], weights=~out_flow, cluster=~o_cbsa), error=function(e) NULL)
if (!is.null(m_iv)) { cat("\nIV (nonattainment -> origin pollution) first-stage F and 2SLS coefficients:\n"); print(summary(m_iv)) }

# capture figure: exposure reduction by race x nonattainment origin
agg <- oy[out_flow>0 & o_race!="No majority", .(re=weighted.mean(re_pm25, out_flow, na.rm=TRUE)), by=.(treat_o, o_race)]
agg[, `:=`(Regime=fifelse(treat_o==1,"Nonattainment origin","Attainment origin"),
           race=factor(o_race, levels=c("Majority white (NH)","Majority Hispanic","Majority Black")))]
gB <- ggplot(agg[!is.na(race)], aes(race, re, fill=Regime)) + geom_col(position="dodge") + geom_hline(yintercept=0, linewidth=.3) +
  scale_fill_manual(values=c(`Nonattainment origin`="#1b7837",`Attainment origin`="#bdbdbd")) +
  labs(title="Exposure change by origin race and PM2.5 nonattainment status",
       subtitle="Flow-weighted exposure reduction of out-migrant streams, by origin race and nonattainment status",
       x=NULL, y=expression("Exposure reduction, "*Delta*" PM2.5 ("*mu*g/m^3*")"), fill=NULL) +
  coord_flip() + theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"), legend.position="top")
ggsave(file.path(fig_dir,"fig_naaqs_capture.png"), gB, width=8, height=4.4, dpi=200)

# ============================================================
# (C) CAPTURE OF CLEANUP at the DESTINATION (in-migration composition)
# ============================================================
# origin racial majority (for classifying in-migrants)
rcs <- list.files(file.path(base_dir,"pollution_bg"), "^race_bg_.*csv$", full.names=TRUE)
rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}")); r<-fread(f,colClasses=c(geoid_bg="character")); r[,year:=y]; r }), fill=TRUE)
rc[, geoid_bg:=pad_bg(geoid_bg)]
rc[, minor := as.integer(sh_black>0.5 | sh_hisp>0.5)]
rc <- unique(rc[, .(year, geoid_bg, minor)], by=c("year","geoid_bg"))
acc <- list()
for (yr in migrate_years) {
  x <- fread(file.path(mig_dir, paste0(yr,".csv")), select=c("Origin","Destination","Flow"))
  x <- x[, .(origin_bg=pad_bg(Origin), dest_bg=pad_bg(Destination), flow=as.numeric(Flow))][flow>0 & origin_bg!=dest_bg]
  x[, `:=`(year=yr, dcounty=substr(dest_bg,1,5))]
  x <- merge(x, rc[year==yr, .(origin_bg=geoid_bg, minor)], by="origin_bg", all.x=TRUE)
  acc[[as.character(yr)]] <- x[!is.na(minor), .(inflow=sum(flow), minor_inflow=sum(flow*minor)), by=.(dcounty,year)]
}
dcap <- rbindlist(acc)[, .(inflow=sum(inflow), minor_inflow=sum(minor_inflow)), by=.(dcounty,year)]
dcap[, minor_share := minor_inflow/inflow]
dcap <- merge(dcap, naa[,.(dcounty=county, treat_d=treat)], by="dcounty", all.x=TRUE); dcap[is.na(treat_d), treat_d:=0L]
dcap <- merge(dcap, cty[,.(dcounty=county, year, dpm=pm25)], by=c("dcounty","year"), all.x=TRUE)
m_dc <- feols(minor_share ~ treat_d + dpm | year, dcap[inflow>0], weights=~inflow, cluster=~dcounty)
dctab <- as.data.table(coeftable(m_dc), keep.rownames="term")
fwrite(dctab, file.path(out_dir,"naaqs_capture_destination.csv"))
cat("\n=== CAPTURE OF CLEANUP (destination): minority in-migrant share vs nonattainment status ===\n")
print(dctab[, .(term, est=round(Estimate,4), se=round(`Std. Error`,4), p=signif(`Pr(>|t|)`,2))])

cat("\nDONE. NAAQS nonattainment outputs in ", out_dir, "\n", sep="")
