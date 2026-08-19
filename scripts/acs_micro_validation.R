# ============================================================
# HOUSEHOLD-LEVEL VALIDATION with ACS MICRODATA (IPUMS USA)
# ------------------------------------------------------------
# Answers the ecological-inference objection with PUBLIC individual-level
# data. The ACS asks every respondent where they lived ONE YEAR AGO
# (MIGPUMA1 + MIGPLAC1) and where they live now (PUMA + STATEFIP). We
# attach PM2.5 to both endpoints and, for movers, compute the change in
# exposure. If minority movers reach systematically dirtier air than
# white movers who START equally polluted -- the same estimand as the
# ecological paper, now at the individual level -- the neighborhood-stream
# result is corroborated for actual households.
#
# Geography note: current residence is a PUMA; prior residence is a
# MIGPUMA (an aggregation of PUMAs). To compare like with like we compute
# pollution at the MIGPUMA level for BOTH endpoints (mapping the current
# PUMA to its MIGPUMA). Resolution is coarser than block groups but far
# finer than state, and requires no restricted data.
#
# WHAT TO DOWNLOAD:
#  (a) IPUMS USA extract (usa.ipums.org), ACS 1-year samples 2012-2019,
#      PERSON records, variables:
#        YEAR SAMPLE SERIAL PERNUM PERWT STATEFIP PUMA
#        MIGRATE1 MIGPLAC1 MIGPUMA1 RACE RACED HISPAN AGE INCTOT EDUCD
#      Export CSV -> MIGRATE_dir/ipums/usa_acs_migration.csv
#  (b) PUMA -> MIGPUMA correspondence (2010 vintage), columns
#      statefip,puma,migpuma  ->  MIGRATE_dir/ipums/puma_migpuma_xwalk.csv
#      (Census "2010 PUMA/MIGPUMA" correspondence, or IPUMS' MIGPUMA
#       composition file.)
#
# Output: pm25_avoidance_lag3_olag1/
#   acs_household_gap.csv       individual-level avoidance gap by race
#   figures/fig_acs_household_gap.png
# ============================================================
suppressPackageStartupMessages({ library(sf); library(data.table); library(stringr); library(fixest); library(splines); library(ggplot2); library(tigris) })
options(tigris_use_cache=TRUE); sf::sf_use_s2(FALSE)
base_dir <- "/Users/priyanka/MIGRATE_dir"
out_dir  <- file.path(base_dir,"pm25_avoidance_lag3_olag1"); fig_dir <- file.path(out_dir,"figures")
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")

# ---- required inputs (hard error if missing) ----
ipf <- file.path(base_dir,"ipums","usa_acs_migration.csv")
xwf <- file.path(base_dir,"ipums","puma_migpuma_xwalk.csv")
stopifnot("IPUMS ACS extract not found: ipums/usa_acs_migration.csv (see header)" = file.exists(ipf))
stopifnot("PUMA->MIGPUMA crosswalk not found: ipums/puma_migpuma_xwalk.csv (see header)" = file.exists(xwf))
xw <- fread(xwf, colClasses="character")
setnames(xw, tolower(names(xw)))
stopifnot("xwalk needs columns statefip,puma,migpuma" = all(c("statefip","puma","migpuma") %in% names(xw)))
xw[, `:=`(statefip=str_pad(statefip,2,pad="0"), puma=str_pad(puma,5,pad="0"), migpuma=str_pad(migpuma,5,pad="0"))]
xw <- unique(xw, by=c("statefip","puma"))

# ---- BG PM2.5 -> MIGPUMA (population-weighted) ----
poll <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character"))
poll[, geoid_bg:=pad_bg(geoid_bg)]; poll <- unique(poll[,.(geoid_bg,year,pm25)], by=c("geoid_bg","year"))
pop <- fread(file.path(base_dir,"acs_gentrification","blockgroup_gentrification_2009_2019_lag3_with_HMDA.csv"),
             colClasses=c(geoid_bg="character"), select=c("year","geoid_bg","total_pop")); pop[, geoid_bg:=pad_bg(geoid_bg)]
# BG -> PUMA by centroid-in-polygon
bg_cache <- file.path(base_dir,"pollution_bg","cache_multi","bg_national.rds")
if (!file.exists(bg_cache)) bg_cache <- file.path(base_dir,"pollution_bg","cache","bg_national.rds")
stopifnot("BG geometry cache not found" = file.exists(bg_cache))
bg <- readRDS(bg_cache); bg <- st_transform(bg, 4326); bg$geoid_bg <- pad_bg(bg$geoid_bg)
cent <- st_centroid(bg[,"geoid_bg"])
pum <- pumas(cb=TRUE, year=2019, class="sf"); pum <- st_transform(pum, 4326)
pn <- names(pum); pcol <- pn[grepl("PUMACE|PUMA",pn)][1]; scol <- pn[grepl("STATEFP",pn)][1]
pum$statefip <- str_pad(pum[[scol]],2,pad="0"); pum$puma <- str_pad(pum[[pcol]],5,pad="0")
j <- st_join(cent, pum[,c("statefip","puma")], join=st_within)
bx <- as.data.table(st_drop_geometry(j))[!is.na(puma), .(geoid_bg, statefip, puma)]
bx <- merge(bx, xw[,.(statefip,puma,migpuma)], by=c("statefip","puma"), all.x=TRUE)
# MIGPUMA-year pollution
pm <- merge(poll, pop, by=c("geoid_bg","year"), all.x=TRUE); pm[is.na(total_pop), total_pop:=1]
pm <- merge(pm, bx[,.(geoid_bg,statefip,migpuma)], by="geoid_bg", allow.cartesian=TRUE)
mp <- pm[!is.na(pm25) & !is.na(migpuma), .(pm25 = sum(pm25*total_pop)/sum(total_pop)), by=.(statefip,migpuma,year)]
setkey(mp, statefip, migpuma, year)

# ---- ACS movers ----
ac <- fread(ipf)
setnames(ac, toupper(names(ac)))
req <- c("YEAR","PERWT","STATEFIP","PUMA","MIGRATE1","MIGPLAC1","MIGPUMA1","RACE","HISPAN")
stopifnot("IPUMS extract missing required columns" = all(req %in% names(ac)))
ac[, `:=`(statefip=str_pad(STATEFIP,2,pad="0"), puma=str_pad(PUMA,5,pad="0"),
          mstate=str_pad(MIGPLAC1,2,pad="0"), migpuma1=str_pad(MIGPUMA1,5,pad="0"))]
# movers within the US who changed house in the last year (MIGRATE1: 2=moved within state?,3=moved between states; 1=same house)
ac <- ac[MIGRATE1 %in% c(2,3) & MIGPLAC1>=1 & MIGPLAC1<=56 & MIGPUMA1>0]
# race groups (IPUMS RACE: 1 White,2 Black; HISPAN>0 = Hispanic)
# Race categories are mutually exclusive: Hispanic (any race) first, then
# NON-HISPANIC SINGLE-RACE Black (RACE==2) and White (RACE==1); all else "Other".
# NOTE: this excludes multiracial (Black alone-or-in-combination) respondents from
# the Black category, so the contrast is NH single-race Black vs NH single-race White.
ac[, race3 := fifelse(HISPAN>0 & HISPAN<9,"Hispanic", fifelse(RACE==2,"Black", fifelse(RACE==1,"White","Other")))]
# destination MIGPUMA = current PUMA -> migpuma
ac <- merge(ac, xw[,.(statefip,puma,dmigpuma=migpuma)], by=c("statefip","puma"), all.x=TRUE)
# attach pollution: destination at survey YEAR, origin at YEAR-1
ac[, yr:=YEAR]
ac <- merge(ac, mp[,.(statefip,migpuma,year,dpm=pm25)], by.x=c("statefip","dmigpuma","yr"),
            by.y=c("statefip","migpuma","year"), all.x=TRUE)
ac[, oyr:=YEAR-1]
ac <- merge(ac, mp[,.(mstate=statefip,migpuma1=migpuma,oyr=year,opm=pm25)],
            by=c("mstate","migpuma1","oyr"), all.x=TRUE)
mv <- ac[!is.na(dpm) & !is.na(opm)]
mv[, delta_pm := dpm - opm]                       # negative = moved to cleaner air
mv[, grp := relevel(factor(race3), ref="White")]
if (!"INCTOT" %in% names(mv)) mv[, INCTOT:=NA_real_]
cat("\nACS movers by race (MIGPUMA changers with pollution both ends):\n"); print(mv[, .N, by=race3])

# ---- individual-level avoidance gap (PERSON records, PERWT-weighted; mirrors the ecological spec) ----
mv[, income_valid := is.finite(INCTOT) & INCTOT>0 & INCTOT<9999999]
# diagnostics (reviewer): prior-state coding match + interstate-mover retention
cat("MIGPLAC1->mstate match to statefip xwalk:",
    round(100*mean(mv$mstate %chin% unique(xw$statefip)),1), "% of movers | interstate movers =",
    mv[mstate!=statefip, .N], "| valid-income share =", round(100*mean(mv$income_valid),1), "%\n")
m_raw <- feols(delta_pm ~ grp, mv, weights=~PERWT, cluster=~dmigpuma)
# primary conditional: origin-PM2.5 spline + year FE (full sample; income NOT required)
m_adj <- feols(delta_pm ~ grp + ns(opm,4) + factor(yr), mv, weights=~PERWT, cluster=~dmigpuma)
# additionally conditional on income, on valid-income records only (sample change made explicit)
m_inc <- if (sum(mv$income_valid) > 500)
  feols(delta_pm ~ grp + ns(opm,4) + log(INCTOT) + factor(yr), mv[income_valid==TRUE], weights=~PERWT, cluster=~dmigpuma) else NULL
grab <- function(fit,lab){ ct<-as.data.table(coeftable(fit),keep.rownames="term")[grepl("grp",term)]
  ct[, term:=sub("grp","",term)][, .(spec=lab, race=term, est=Estimate, se=`Std. Error`, p=`Pr(>|t|)`)] }
res <- rbindlist(c(list(grab(m_raw,"raw"), grab(m_adj,"conditional")),
                   if (!is.null(m_inc)) list(grab(m_inc,"conditional+income"))), fill=TRUE)
res[, `:=`(est=round(est,3), se=round(se,3), p=signif(p,2))]
fwrite(res, file.path(out_dir,"acs_household_gap.csv"))
cat("\n=== ACS INDIVIDUAL-LEVEL avoidance gap (positive = reached dirtier air than white movers) ===\n"); print(res)

# ---- figure ----
pf <- res[spec=="conditional"]; pf[, race:=factor(race, levels=c("Black","Hispanic","Other"))]
p <- ggplot(pf[!is.na(race)], aes(race, est, fill=race)) +
  geom_col() + geom_errorbar(aes(ymin=est-1.96*se, ymax=est+1.96*se), width=.2) + geom_hline(yintercept=0, linewidth=.3) +
  scale_fill_manual(values=c(Black="#1b7837",Hispanic="#762a83",Other="#999999"), guide="none") +
  labs(title="Individual-level validation (ACS microdata): minority movers reached dirtier air",
       subtitle="Change in MIGPUMA PM2.5 on moving, vs white movers, net of an origin-PM2.5 spline and year (income-adjusted version in the table)",
       x="Race (ref: White)", y=expression(Delta*" PM2.5 vs white movers ("*mu*g/m^3*")")) +
  theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"))
ggsave(file.path(fig_dir,"fig_acs_household_gap.png"), p, width=7.5, height=4.2, dpi=200)

# ============================================================
# (S19) ACS SUBSAMPLE SENSITIVITY: same-state vs interstate, income, tenure
# ============================================================
mv[, interstate := as.integer(mstate != statefip)]
med_inc <- suppressWarnings(median(mv$INCTOT[is.finite(mv$INCTOT) & mv$INCTOT>0], na.rm=TRUE))
mv[, hi_inc := as.integer(is.finite(INCTOT) & INCTOT > med_inc)]
subs <- list(`same-state`=mv[interstate==0], `interstate`=mv[interstate==1],
             `below-median income`=mv[hi_inc==0], `above-median income`=mv[hi_inc==1])
if ("OWNERSHP" %in% names(mv)) { subs[["owners"]] <- mv[OWNERSHP==1]; subs[["renters"]] <- mv[OWNERSHP==2] }
covs2 <- "ns(opm,4) + factor(yr)"   # subsample robustness of the gap; income not required (avoids the log(pmax(.,1)) hack)
subres <- rbindlist(lapply(names(subs), function(nm){
  d <- subs[[nm]]; if (nrow(d) < 500) return(NULL)
  f <- tryCatch(feols(as.formula(paste("delta_pm ~ grp +", covs2)), d, weights=~PERWT, cluster=~dmigpuma), error=function(e) NULL)
  if (is.null(f)) return(NULL)
  ct <- as.data.table(coeftable(f), keep.rownames="term")[grepl("grp",term)]
  ct[grepl("Black|Hispanic", term), .(subsample=nm, race=sub("grp","",term),
       est=round(Estimate,3), se=round(`Std. Error`,3), p=signif(`Pr(>|t|)`,2), n=nrow(d))] }), fill=TRUE)
fwrite(subres, file.path(out_dir,"acs_household_subsamples.csv"))
cat("\n=== (S19) ACS SUBSAMPLE SENSITIVITY (conditional gap) ===\n"); print(subres)
if (!"OWNERSHP" %in% names(mv))
  cat("(tenure split skipped: add OWNERSHP to the IPUMS extract via acs_ipums_pull.R and re-pull.)\n")

# ============================================================
# (S21) MULTIPOLLUTANT household-level gaps (all eight pollutants)
# ============================================================
poll2 <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character"))
poll2[, geoid_bg:=pad_bg(geoid_bg)]; poll2 <- unique(poll2, by=c("geoid_bg","year"))
cf2 <- file.path(base_dir,"pollution_bg","bg_pm25_components_panel.csv")
if (file.exists(cf2)) { cp2<-fread(cf2,colClasses=c(geoid_bg="character")); cp2[,geoid_bg:=pad_bg(geoid_bg)]
  cc2<-intersect(c("so4","no3","bc","om"),names(cp2))
  poll2<-merge(poll2, unique(cp2[,c("geoid_bg","year",cc2),with=FALSE],by=c("geoid_bg","year")), by=c("geoid_bg","year"), all.x=TRUE) }
P2 <- intersect(c("pm25","no2","ozone","wildfire_pm25","so4","no3","bc","om"), names(poll2))
pm2 <- merge(poll2, pop, by=c("geoid_bg","year"), all.x=TRUE); pm2[is.na(total_pop), total_pop:=1]
pm2 <- merge(pm2, bx[,.(geoid_bg,statefip,migpuma)], by="geoid_bg", allow.cartesian=TRUE)
mp2 <- pm2[, c(lapply(.SD, function(v) sum(v*total_pop, na.rm=TRUE)/sum(total_pop[!is.na(v)])),
               list(.dummy=1)), by=.(statefip,migpuma,year), .SDcols=P2][, .dummy:=NULL]
res2 <- list()
for (p in P2) {
  dsub <- mp2[, c("statefip","migpuma","year",p), with=FALSE]
  mvp <- copy(mv[, .(mstate,migpuma1,oyr,statefip,dmigpuma,yr,grp,PERWT,INCTOT,dmig=dmigpuma)])
  mvp <- merge(mvp, setnames(copy(dsub), c("statefip","migpuma","year",p), c("statefip","dmigpuma","yr","dp")),
               by=c("statefip","dmigpuma","yr"), all.x=TRUE)
  mvp <- merge(mvp, setnames(copy(dsub), c("statefip","migpuma","year",p), c("mstate","migpuma1","oyr","op")),
               by=c("mstate","migpuma1","oyr"), all.x=TRUE)
  mvp <- mvp[!is.na(dp) & !is.na(op)][, delta:=dp-op]
  fit <- tryCatch(feols(delta ~ grp + ns(op,4) + log(pmax(INCTOT,1)) + factor(yr), mvp, weights=~PERWT, cluster=~dmig), error=function(e) NULL)
  if (is.null(fit)) next
  ct <- as.data.table(coeftable(fit), keep.rownames="term")[grepl("grpBlack|grpHispanic",term)]
  res2[[p]] <- ct[, .(pollutant=p, race=sub("grp","",term), est=Estimate, se=`Std. Error`, p=`Pr(>|t|)`)]
}
hh_mp <- rbindlist(res2)
fwrite(hh_mp, file.path(out_dir,"acs_household_gap_multipollutant.csv"))
cat("\n=== (S21) ACS HOUSEHOLD-LEVEL gaps across pollutants (conditional; vs white movers) ===\n")
print(hh_mp[, .(pollutant, race, est=round(est,4), se=round(se,4), p=signif(p,2))])

cat("\nDONE. ACS individual-level validation in ", out_dir, "\n", sep="")
