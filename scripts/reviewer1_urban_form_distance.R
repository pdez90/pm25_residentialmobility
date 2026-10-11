# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# REVIEWER 1 (revision, Oct 2026): urban form and move distance
# ------------------------------------------------------------
# R1.1  Do the clean-air gains from moving occur in the periphery of large
#       urban centres? -> classify every block group by position within its
#       metropolitan area (urban core / inner suburb / periphery / micropolitan
#       / nonmetro), tabulate origin-type x destination-type flows and the
#       mean PM2.5 reallocation in each cell by origin race/ethnicity, and
#       re-estimate the conditional avoidance gap WITHIN origin-destination
#       type (O-D-type fixed effects).
# R1.2  How far do movers go, and does the gap depend on distance? ->
#       flow-weighted distance distribution by origin group, and the
#       conditional gap estimated within distance bands
#       (0-10, 10-50, 50-200, >200 km), plus a pooled interaction model.
#
# Urban-form classes (2019 cb geographies):
#   Urban core    = block group inside the CBSA's first-named principal city
#   Inner suburb  = other metro BG in the principal city's county, or with
#                   centroid within 15 km of the principal city's centroid
#   Periphery     = any other block group in a metropolitan (M1) CBSA
#   Micropolitan  = block group in a micropolitan (M2) CBSA
#   Nonmetro      = not in any CBSA
# Principal city = the first city named in the CBSA title (e.g. "Denver" in
# "Denver-Aurora-Lakewood, CO"), matched to a tigris place in the CBSA's
# states. CBSAs whose principal city cannot be matched are reported and
# their block groups labelled "Metro (unmatched)".
#
# Inputs : bg_pollution_panel.csv, the ACS/HMDA panel, race_bg_*.csv,
#          county_cbsa_xwalk.csv, bg_national.rds, MIGRATE <year>.csv
# Outputs: pm25_avoidance_lag3_olag1/r1_*.csv and figures/r1_*.png
# Runtime: one pass over the nine MIGRATE files (~20-25 min).
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr); library(fixest); library(splines)
  library(ggplot2); library(sf); library(tigris) })
sf::sf_use_s2(FALSE); options(tigris_use_cache=TRUE)
base_dir <- migrate_data_dir(); mig_dir <- base_dir
gent_dir <- file.path(base_dir,"acs_gentrification"); change_lag <- 3
out_dir <- file.path(base_dir,"pm25_avoidance_lag3_olag1"); fig_dir <- file.path(out_dir,"figures")
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
migrate_years <- 2011:2019
INNER_KM <- 15
BANDS <- c(0,10,50,200,Inf); BAND_LAB <- c("0-10 km","10-50 km","50-200 km",">200 km")
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")
wmean  <- function(x,w){ok<-!is.na(x)&!is.na(w); if(!any(ok)) return(NA_real_); sum(x[ok]*w[ok])/sum(w[ok])}
wq     <- function(x,w,p){ok<-!is.na(x)&!is.na(w)&w>0; x<-x[ok]; w<-w[ok]; o<-order(x); x<-x[o]; w<-w[o]; cw<-cumsum(w)/sum(w); x[which(cw>=p)[1]]}
hav <- function(lo1,la1,lo2,la2){ R<-6371; dlo<-(lo2-lo1)*pi/180; dla<-(la2-la1)*pi/180
  a<-sin(dla/2)^2+cos(la1*pi/180)*cos(la2*pi/180)*sin(dlo/2)^2; 2*R*asin(pmin(1,sqrt(a))) }
ct <- function(fit, term="grpMajority Black") tryCatch({ x<-coeftable(fit)[term,]; data.table(est=x[["Estimate"]], se=x[["Std. Error"]], p=x[["Pr(>|t|)"]]) },
                                                       error=function(e) data.table(est=NA_real_,se=NA_real_,p=NA_real_))

# ============================================================
# 1. Block-group panel (same construction as affordability_robustness.R)
# ============================================================
poll <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character"), select=c("geoid_bg","year","pm25"))
poll[, geoid_bg:=pad_bg(geoid_bg)]; poll <- unique(poll, by=c("geoid_bg","year"))
cf  <- file.path(gent_dir, paste0("blockgroup_gentrification_2009_2019_lag",change_lag,"_with_HMDA.csv"))
acs <- fread(cf, colClasses=c(geoid_bg="character", state_county="character"),
             select=c("year","geoid_bg","state_county","med_hh_income","poverty_pct","bach_pct","total_pop"))
acs[, geoid_bg:=pad_bg(geoid_bg)]
rcs <- list.files(file.path(base_dir,"pollution_bg"), "^race_bg_.*csv$", full.names=TRUE)
rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}")); r<-fread(f,colClasses=c(geoid_bg="character")); r[,year:=y]; r }), fill=TRUE)
rc[, geoid_bg:=pad_bg(geoid_bg)]
acs <- merge(acs, rc[,.(year,geoid_bg,sh_black,sh_hisp,sh_white)], by=c("year","geoid_bg"), all.x=TRUE)
acs[, race_majority := fifelse(sh_black>0.5,"Majority Black", fifelse(sh_hisp>0.5,"Majority Hispanic", fifelse(sh_white>0.5,"Majority white (NH)","No majority")))]
acs <- merge(acs, poll, by=c("geoid_bg","year"), all.x=TRUE)
acs[, county := substr(geoid_bg,1,5)]
cxw <- fread(file.path(base_dir,"pollution_bg","county_cbsa_xwalk.csv"), colClasses=c(county="character",cbsa="character"))
acs <- merge(acs, unique(cxw,by="county"), by="county", all.x=TRUE); acs[is.na(cbsa), cbsa:=paste0("nonmetro_",county)]
acs <- unique(acs, by=c("year","geoid_bg"))

# ---- centroids ----
bg_cache <- file.path(base_dir,"pollution_bg","cache_multi","bg_national.rds")
if (!file.exists(bg_cache)) bg_cache <- file.path(base_dir,"pollution_bg","cache","bg_national.rds")
stopifnot("bg_national.rds not found" = file.exists(bg_cache))
bg <- readRDS(bg_cache); bg <- st_transform(bg, 4326)
cc <- st_coordinates(st_centroid(bg)); cent <- data.table(geoid_bg=pad_bg(bg$geoid_bg), lon=cc[,1], lat=cc[,2])
cent[, county := substr(geoid_bg,1,5)]
cent <- merge(cent, unique(cxw,by="county"), by="county", all.x=TRUE)

# ============================================================
# 2. Urban-form classification
# ============================================================
cat("\n=== Urban-form classification (2019 cb geographies) ===\n")
cb <- st_transform(core_based_statistical_areas(cb=TRUE, year=2019), 4326)
cbd <- data.table(cbsa=cb$CBSAFP, cbsa_name=cb$NAME, lsad=cb$LSAD)
cbd[, pc_name := str_trim(str_split_fixed(cbsa_name, "[-,/]", 2)[,1])]
cbd[, states := str_trim(str_split_fixed(cbsa_name, ",", 2)[,2])]
# places in every state (cb=TRUE), matched to CBSA by containing county -> cbsa
sts <- unique(as.data.table(tigris::fips_codes)[as.integer(state_code)<=56][, .(state_code, state)])
pl <- do.call(rbind, lapply(sts$state_code, function(s) tryCatch(places(state=s, cb=TRUE, year=2019)[,c("GEOID","NAME","STATEFP")], error=function(e) NULL)))
pl <- st_transform(pl, 4326)
plc <- st_coordinates(st_centroid(pl))
pld <- data.table(place_geoid=pl$GEOID, pname=pl$NAME, statefp=pl$STATEFP, plon=plc[,1], plat=plc[,2])
pld <- merge(pld, sts[, .(statefp=state_code, st=state)], by="statefp")
# a place belongs to the CBSA of the county its centroid falls in
pc_sf <- st_as_sf(pld, coords=c("plon","plat"), crs=4326, remove=FALSE)
co <- st_transform(counties(cb=TRUE, year=2019)[,"GEOID"], 4326)
wi <- st_within(pc_sf, co); pld$county <- co$GEOID[vapply(wi, function(i) if (length(i)) i[1] else NA_integer_, 1L)]
pld <- merge(pld, unique(cxw,by="county"), by="county", all.x=TRUE)
# match principal city: exact NAME within the CBSA, state among the CBSA's states
cbd[, st1 := str_trim(str_split_fixed(states,"-",2)[,1])]
pm <- merge(cbd[lsad=="M1"], pld[, .(cbsa, pname, st, place_geoid, plon, plat, pc_county=county)], by="cbsa", allow.cartesian=TRUE)
pm <- pm[tolower(pname)==tolower(pc_name) & st==st1]
pm <- pm[, .SD[1], by=cbsa]            # one principal city per CBSA
cat(sprintf("  metropolitan CBSAs: %d; principal city matched by name: %d\n", cbd[lsad=="M1",.N], nrow(pm)))
unmatched <- cbd[lsad=="M1" & !cbsa %in% pm$cbsa]
if (nrow(unmatched)) cat("  unmatched (labelled 'Metro (unmatched)'):", paste(head(unmatched$cbsa_name,15), collapse="; "), if(nrow(unmatched)>15) "..." else "", "\n")
fwrite(pm[, .(cbsa, cbsa_name, principal_city=pname, place_geoid, pc_county)], file.path(out_dir,"r1_principal_cities.csv"))

# BG inside the principal city polygon?
pc_poly <- pl[pl$GEOID %in% pm$place_geoid, ]
bgc_sf  <- st_as_sf(cent, coords=c("lon","lat"), crs=4326, remove=FALSE)
cent[, in_pc := lengths(st_within(bgc_sf, pc_poly)) > 0]
cent <- merge(cent, cbd[, .(cbsa, lsad)], by="cbsa", all.x=TRUE)
cent <- merge(cent, pm[, .(cbsa, plon, plat, pc_county)], by="cbsa", all.x=TRUE)
cent[, d_pc := hav(lon,lat,plon,plat)]
cent[, uform := fcase(
  is.na(lsad),                                              "Nonmetro",
  lsad %in% "M2",                                           "Micropolitan",
  lsad %in% "M1" & is.na(plon),                             "Metro (unmatched)",
  lsad %in% "M1" & in_pc,                                   "Urban core",
  lsad %in% "M1" & ((!is.na(pc_county) & county==pc_county) | (!is.na(d_pc) & d_pc<=INNER_KM)), "Inner suburb",
  default =                                                 "Periphery")]
cent[, uform := factor(uform, levels=c("Urban core","Inner suburb","Periphery","Micropolitan","Nonmetro","Metro (unmatched)"))]
cat("\n  block groups by urban form:\n"); print(cent[, .N, by=uform][order(uform)])
fwrite(cent[, .(geoid_bg, cbsa, uform, dist_to_principal_city_km=round(d_pc,1))], file.path(out_dir,"r1_bg_urban_form.csv"))

# ============================================================
# 3. MIGRATE pass: origin-year x destination-type and x distance band
# ============================================================
oa <- acs[, .(cy=year, origin_bg=geoid_bg, o_race=race_majority, o_inc0=med_hh_income, o_pov=poverty_pct,
              o_bach=bach_pct, o_county=county, o_cbsa=cbsa, o_pm25=pm25)]
oa <- merge(oa, cent[,.(origin_bg=geoid_bg,o_lon=lon,o_lat=lat,o_uform=uform)], by="origin_bg", all.x=TRUE); setkey(oa, cy, origin_bg)
da <- acs[, .(cy=year, dest_bg=geoid_bg, d_pm25=pm25, d_cbsa=cbsa)]
da <- merge(da, cent[,.(dest_bg=geoid_bg,d_lon=lon,d_lat=lat,d_uform=uform)], by="dest_bg", all.x=TRUE); setkey(da, cy, dest_bg)

acc_u <- list(); acc_d <- list(); acc_q <- list()
for (yr in migrate_years) {
  cat("  MIGRATE", yr, "...\n")
  x <- fread(file.path(mig_dir, paste0(yr,".csv")), select=c("Origin","Destination","Flow"))
  x <- x[, .(origin_bg=pad_bg(Origin), dest_bg=pad_bg(Destination), flow=as.numeric(Flow))][flow>0 & origin_bg!=dest_bg]
  x[, cy:=yr]; setkey(x,cy,origin_bg); x<-oa[x,on=.(cy,origin_bg)]; setkey(x,cy,dest_bg); x<-da[x,on=.(cy,dest_bg)]
  x <- x[!is.na(o_race) & !is.na(o_pm25) & !is.na(d_pm25)]
  x[, re := d_pm25 - o_pm25]
  x[, dist := hav(o_lon,o_lat,d_lon,d_lat)]
  x[, band := cut(dist, BANDS, labels=BAND_LAB, right=TRUE, include.lowest=TRUE)]
  byv <- c("origin_bg","o_race","o_inc0","o_pov","o_bach","o_county","o_cbsa","o_uform")
  acc_u[[as.character(yr)]] <- x[!is.na(d_uform), .(year=yr, flow=sum(flow), re_pm25=wmean(re,flow), orig_pm25=wmean(o_pm25,flow), mean_dist=wmean(dist,flow)), by=c(byv,"d_uform")]
  acc_d[[as.character(yr)]] <- x[!is.na(band),   .(year=yr, flow=sum(flow), re_pm25=wmean(re,flow), orig_pm25=wmean(o_pm25,flow)), by=c(byv,"band")]
  # distance quantiles by origin group (flow-weighted), computed on the flow level
  acc_q[[as.character(yr)]] <- x[!is.na(dist), .(year=yr, flow=sum(flow), p25=wq(dist,flow,.25), p50=wq(dist,flow,.50), p75=wq(dist,flow,.75), p90=wq(dist,flow,.90),
                                                 share_le10=wmean(dist<=10,flow), share_le50=wmean(dist<=50,flow), share_within_cbsa=wmean(o_cbsa==d_cbsa,flow)), by=o_race]
  rm(x); gc()
}
ou <- rbindlist(acc_u); od <- rbindlist(acc_d); oq <- rbindlist(acc_q)
ou <- ou[flow>0 & !is.na(o_inc0) & o_inc0>0]; od <- od[flow>0 & !is.na(o_inc0) & o_inc0>0]
ou[, grp := relevel(factor(o_race), ref="Majority white (NH)")]; od[, grp := relevel(factor(o_race), ref="Majority white (NH)")]
fwrite(ou, file.path(out_dir,"r1_origin_year_by_dest_type.csv")); fwrite(od, file.path(out_dir,"r1_origin_year_by_distance_band.csv"))
ctrl <- "+ ns(orig_pm25,4) + log(o_inc0) + o_pov + o_bach"

# ============================================================
# 4. R1.1  Origin-type x destination-type flow matrix and reallocation
# ============================================================
cat("\n=== R1.1 Origin type x destination type: share of out-flow and mean PM2.5 reallocation ===\n")
mat <- ou[!is.na(o_uform), .(flow=sum(flow), re_pm25=wmean(re_pm25,flow)), by=.(o_race, o_uform, d_uform)]
mat[, share_of_origin_flow := flow/sum(flow), by=.(o_race, o_uform)]
mat <- mat[order(o_race, o_uform, d_uform)]
fwrite(mat, file.path(out_dir,"r1_od_type_matrix.csv"))
# headline: core origins, where did they go, by group
core <- mat[o_uform=="Urban core" & o_race %in% c("Majority white (NH)","Majority Black","Majority Hispanic")]
print(dcast(core, o_race ~ d_uform, value.var="share_of_origin_flow", fun.aggregate=function(v) round(100*v,1)))
cat("  (percent of out-flow from URBAN-CORE origins reaching each destination type)\n")
print(dcast(core, o_race ~ d_uform, value.var="re_pm25", fun.aggregate=function(v) round(v,3)))
cat("  (mean PM2.5 reallocation, ug/m3, of those moves; negative = cleaner)\n")

# ---- conditional gap with O-D-type fixed effects ----
cat("\n=== R1.1 Conditional avoidance gap: main vs. within origin-destination type ===\n")
ou[, od_type := paste(o_uform, d_uform, sep=" -> ")]
ouM <- ou[!o_uform %in% "Metro (unmatched)" & !d_uform %in% "Metro (unmatched)"]
m_main <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," | o_cbsa + year")),          ouM, weights=~flow, cluster=~o_cbsa)
m_od   <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," | o_cbsa + year + od_type")), ouM, weights=~flow, cluster=~o_cbsa)
m_dt   <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," | o_cbsa + year + d_uform")), ouM, weights=~flow, cluster=~o_cbsa)
odfe <- rbindlist(lapply(c("grpMajority Black","grpMajority Hispanic","grpNo majority"), function(t) rbindlist(list(
  cbind(spec="CBSA + year FE (main, on O-D-type panel)", group=sub("grp","",t), ct(m_main,t)),
  cbind(spec="+ destination-type FE",                     group=sub("grp","",t), ct(m_dt,t)),
  cbind(spec="+ origin x destination type FE",            group=sub("grp","",t), ct(m_od,t))))))
odfe[, `:=`(lo=est-1.96*se, hi=est+1.96*se)]
print(odfe[, .(spec, group, est=round(est,4), lo=round(lo,4), hi=round(hi,4), p=signif(p,2))])
fwrite(odfe, file.path(out_dir,"r1_gap_od_type_fe.csv"))
# gap by destination type (separate regressions)
cat("\n--- Black gap by destination type (separate models, CBSA + year FE) ---\n")
gdt <- rbindlist(lapply(levels(droplevels(ouM$d_uform)), function(dt){ d<-ouM[d_uform==dt]; if (nrow(d)<500) return(NULL)
  f <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," | o_cbsa + year")), d, weights=~flow, cluster=~o_cbsa)
  cbind(dest_type=dt, n=nrow(d), flow=sum(d$flow), ct(f)) }))
gdt[, `:=`(lo=est-1.96*se, hi=est+1.96*se)]; print(gdt[, .(dest_type, n, est=round(est,4), lo=round(lo,4), hi=round(hi,4), p=signif(p,2))])
fwrite(gdt, file.path(out_dir,"r1_black_gap_by_dest_type.csv"))

# ============================================================
# 5. R1.2  Move distance
# ============================================================
cat("\n=== R1.2 Move distance by origin group (flow-weighted, pooled 2011-2019) ===\n")
dq <- oq[, .(p25_km=wmean(p25,flow), median_km=wmean(p50,flow), p75_km=wmean(p75,flow), p90_km=wmean(p90,flow),
             share_le10km=wmean(share_le10,flow), share_le50km=wmean(share_le50,flow), share_within_cbsa=wmean(share_within_cbsa,flow)), by=o_race]
# (flow-weighted means of the yearly quantiles; exact pooled quantiles would need all flow-level rows in memory)
print(dq[, lapply(.SD, function(v) if (is.numeric(v)) round(v,3) else v)])
fwrite(dq, file.path(out_dir,"r1_distance_distribution_by_group.csv"))

cat("\n=== R1.2 Conditional Black gap within distance bands (separate models, CBSA + year FE) ===\n")
gb <- rbindlist(lapply(BAND_LAB, function(b){ d<-od[band==b]; if (nrow(d)<500) return(NULL)
  f <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," | o_cbsa + year")), d, weights=~flow, cluster=~o_cbsa)
  rbindlist(lapply(c("grpMajority Black","grpMajority Hispanic"), function(t) cbind(band=b, group=sub("grp","",t), n=nrow(d), flow=sum(d$flow), ct(f,t)))) }))
gb[, `:=`(lo=est-1.96*se, hi=est+1.96*se)]
print(gb[, .(band, group, n, est=round(est,4), lo=round(lo,4), hi=round(hi,4), p=signif(p,2))])
fwrite(gb, file.path(out_dir,"r1_gap_by_distance_band.csv"))
# pooled interaction model: does the gap differ by band?
m_int <- feols(as.formula(paste0("re_pm25 ~ grp*band",ctrl," | o_cbsa + year")), od, weights=~flow, cluster=~o_cbsa)
w <- tryCatch(wald(m_int, "grpMajority Black:band"), error=function(e) NULL)
if (!is.null(w)) cat(sprintf("  Wald test, Black gap equal across bands: F = %.2f, p = %.3g\n", w$stat, w$p))

# ============================================================
# 6. Figures
# ============================================================
pdat <- rbind(gb[group=="Majority Black", .(panel="By move distance", x=band, est, lo, hi)],
              gdt[, .(panel="By destination type", x=dest_type, est, lo, hi)])
pdat[, x := factor(x, levels=c(BAND_LAB, levels(cent$uform)))]
p1 <- ggplot(pdat, aes(x, est)) + geom_hline(yintercept=0, linetype=2, colour="grey50") +
  geom_pointrange(aes(ymin=lo, ymax=hi)) + facet_wrap(~panel, scales="free_x") +
  labs(x=NULL, y=expression("Conditional Black avoidance gap ("*mu*"g/m"^3*")"),
       title="Conditional PM2.5 avoidance gap, majority-Black vs majority-white origins",
       subtitle="Separate models within each subset; CBSA + year FE, origin-PM2.5 spline, SES controls; 95% CI") +
  theme_minimal(base_size=11) + theme(axis.text.x=element_text(angle=20, hjust=1))
ggsave(file.path(fig_dir,"r1_gap_by_distance_and_dest_type.png"), p1, width=9, height=4.5, dpi=200)
m2 <- mat[o_race %in% c("Majority white (NH)","Majority Black","Majority Hispanic") & !o_uform %in% "Metro (unmatched)" & !d_uform %in% "Metro (unmatched)"]
p2 <- ggplot(m2, aes(d_uform, o_uform, fill=100*share_of_origin_flow)) + geom_tile() +
  geom_text(aes(label=sprintf("%.0f%%\n%+.2f", 100*share_of_origin_flow, re_pm25)), size=2.6) +
  facet_wrap(~o_race) + scale_fill_gradient(low="white", high="steelblue", name="% of out-flow") +
  labs(x="Destination type", y="Origin type", title="Where moves go, and what they gain: share of out-flow and mean PM2.5 reallocation (ug/m3)") +
  theme_minimal(base_size=10) + theme(axis.text.x=element_text(angle=25, hjust=1))
ggsave(file.path(fig_dir,"r1_od_type_matrix.png"), p2, width=12, height=4.8, dpi=200)
cat("\nDONE. Outputs: r1_*.csv in", out_dir, "and figures/r1_*.png\n")
