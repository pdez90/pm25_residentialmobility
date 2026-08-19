# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# MULTI-POLLUTANT ROBUSTNESS: affordability, selection, opportunity,
# segregation — for ALL EIGHT pollutants (stage 8b)
# ------------------------------------------------------------
# Extends the PM2.5-only robustness suite so ozone, NO2, wildfire smoke,
# and the PM2.5 components do not "disappear" after Figure 5:
#   (A) Gelbach affordability: gap_base vs gap + destination home value/rent
#   (B) Selection: all / within-CBSA / short-distance (<=50 km) / +distance
#   (C) Opportunity set: clean = below metro-median OF THAT pollutant
#   (D) Segregation: Black gap x standardized dissimilarity interaction
# Outputs (pm25_avoidance_lag3_olag1/):
#   multipollutant_affordability.csv  multipollutant_selection.csv
#   multipollutant_opportunity.csv    multipollutant_segregation.csv
#   figures/figS_mp_affordability.png figures/figS_mp_selection.png
#   figures/figS_mp_opportunity.png   figures/figS_mp_segregation.png
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr); library(fixest); library(splines); library(ggplot2); library(sf) })
base_dir <- migrate_data_dir(); mig_dir <- base_dir
out_dir <- file.path(base_dir,"pm25_avoidance_lag3_olag1"); fig_dir <- file.path(out_dir,"figures")
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
migrate_years <- 2011:2019; AFFORD_MULT <- 3; SHORT_KM <- 50
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")
wmean <- function(x,w){ok<-!is.na(x)&!is.na(w)&w>0; if(!any(ok)) NA_real_ else sum(x[ok]*w[ok])/sum(w[ok])}
PLAB <- c(pm25="PM2.5", no2="NO2", ozone="Ozone", bc="Black carbon", om="Organic matter",
          so4="Sulfate", no3="Nitrate", wildfire_pm25="Wildfire PM2.5")

# ---- pollutant panel (all eight) ----
poll <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character"))
poll[, geoid_bg:=pad_bg(geoid_bg)]; poll <- unique(poll, by=c("geoid_bg","year"))
cf <- file.path(base_dir,"pollution_bg","bg_pm25_components_panel.csv")
stopifnot("components panel missing" = file.exists(cf))
cp <- fread(cf, colClasses=c(geoid_bg="character")); cp[, geoid_bg:=pad_bg(geoid_bg)]
cc <- intersect(c("so4","no3","bc","om"), names(cp))
poll <- merge(poll, unique(cp[, c("geoid_bg","year",cc), with=FALSE], by=c("geoid_bg","year")), by=c("geoid_bg","year"), all.x=TRUE)
POLLS <- intersect(names(PLAB), names(poll))
cat("Pollutants:", paste(POLLS, collapse=", "), "\n")

# ---- covariates: income/pov/bach/housing + race + cbsa + centroids ----
gf <- file.path(base_dir,"acs_gentrification","blockgroup_gentrification_2009_2019_lag3_with_HMDA.csv")
acs <- fread(gf, colClasses=c(geoid_bg="character"),
             select=c("year","geoid_bg","med_hh_income","poverty_pct","bach_pct","med_home_value","med_gross_rent","total_pop"))
acs[, geoid_bg:=pad_bg(geoid_bg)]
rcs <- list.files(file.path(base_dir,"pollution_bg"), "^race_bg_.*csv$", full.names=TRUE)
rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}")); r<-fread(f,colClasses=c(geoid_bg="character")); r[,year:=y]; r }), fill=TRUE)
rc[, geoid_bg:=pad_bg(geoid_bg)]
rc[, o_race:=fifelse(sh_black>0.5,"Majority Black", fifelse(sh_hisp>0.5,"Majority Hispanic",
              fifelse(sh_white>0.5,"Majority white (NH)","No majority")))]
acs <- merge(acs, rc[,.(year,geoid_bg,o_race,sh_black,sh_white)], by=c("year","geoid_bg"), all.x=TRUE)
acs[, county:=substr(geoid_bg,1,5)]
xw <- fread(file.path(base_dir,"pollution_bg","county_cbsa_xwalk.csv"), colClasses=c(county="character",cbsa="character"))
acs <- merge(acs, unique(xw,by="county"), by="county", all.x=TRUE); acs[is.na(cbsa), cbsa:=paste0("nonmetro_",county)]
acs <- unique(acs, by=c("geoid_bg","year"))
# centroids (cache)
cenf <- file.path(base_dir,"pollution_bg","bg_centroids.csv")
if (!file.exists(cenf)) {
  bgc <- file.path(base_dir,"pollution_bg","cache_multi","bg_national.rds")
  if (!file.exists(bgc)) bgc <- file.path(base_dir,"pollution_bg","cache","bg_national.rds")
  bg <- readRDS(bgc); sf::sf_use_s2(FALSE)
  ct <- suppressWarnings(st_coordinates(st_centroid(bg)))
  fwrite(data.table(geoid_bg=pad_bg(bg$geoid_bg), lon=ct[,1], lat=ct[,2]), cenf)
}
cen <- fread(cenf, colClasses=c(geoid_bg="character")); cen[, geoid_bg:=pad_bg(geoid_bg)]
hav <- function(lo1,la1,lo2,la2){ r<-6371; p<-pi/180
  a<-sin((la2-la1)*p/2)^2+cos(la1*p)*cos(la2*p)*sin((lo2-lo1)*p/2)^2; 2*r*asin(pmin(1,sqrt(a))) }

# ---- panel with everything, keyed for flow joins ----
pan <- merge(poll, acs, by=c("geoid_bg","year"), all.x=FALSE)
pan <- merge(pan, cen, by="geoid_bg", all.x=TRUE)
oa <- pan[, c("year","geoid_bg","o_race","med_hh_income","poverty_pct","bach_pct","cbsa","lon","lat", POLLS), with=FALSE]
setnames(oa, c("year","geoid_bg","med_hh_income","poverty_pct","bach_pct","cbsa","lon","lat", POLLS),
             c("cy","origin_bg","o_inc0","o_pov","o_bach","o_cbsa","o_lon","o_lat", paste0("o_",POLLS)))
setkey(oa, cy, origin_bg)
da <- pan[, c("year","geoid_bg","cbsa","med_home_value","med_gross_rent","lon","lat", POLLS), with=FALSE]
setnames(da, c("year","geoid_bg","cbsa","med_home_value","med_gross_rent","lon","lat", POLLS),
             c("cy","dest_bg","d_cbsa","d_homeval","d_rent","d_lon","d_lat", paste0("d_",POLLS)))
setkey(da, cy, dest_bg)

# ---- one MIGRATE pass -> origin-year reallocation for every pollutant ----
acc <- list()
for (yr in migrate_years) {
  x <- fread(file.path(mig_dir, paste0(yr,".csv")), select=c("Origin","Destination","Flow"))
  x <- x[, .(origin_bg=pad_bg(Origin), dest_bg=pad_bg(Destination), flow=as.numeric(Flow))][flow>0 & origin_bg!=dest_bg]
  x[, cy:=yr]; x <- oa[x, on=.(cy,origin_bg)]; x <- da[x, on=.(cy,dest_bg)]
  x <- x[!is.na(o_cbsa)]
  x[, dist := hav(o_lon,o_lat,d_lon,d_lat)]
  x[, wc := as.integer(!is.na(d_cbsa) & d_cbsa==o_cbsa)]
  x[, sd50 := as.integer(!is.na(dist) & dist<=SHORT_KM)]
  ag <- x[, c(list(out_flow=sum(flow),
                   dest_homeval=wmean(d_homeval,flow), dest_rent=wmean(d_rent,flow),
                   mean_dist=wmean(dist,flow),
                   out_flow_wc=sum(flow*wc), out_flow_sd=sum(flow*sd50)),
             setNames(lapply(POLLS, function(p) wmean(get(paste0("d_",p))-get(paste0("o_",p)), flow)), paste0("re_",POLLS)),
             setNames(lapply(POLLS, function(p) wmean(fifelse(wc==1, get(paste0("d_",p))-get(paste0("o_",p)), NA_real_), fifelse(wc==1,flow,NA_real_))), paste0("rewc_",POLLS)),
             setNames(lapply(POLLS, function(p) wmean(fifelse(sd50==1, get(paste0("d_",p))-get(paste0("o_",p)), NA_real_), fifelse(sd50==1,flow,NA_real_))), paste0("resd_",POLLS)),
             setNames(lapply(POLLS, function(p) get(paste0("o_",p))[1]), paste0("o_",POLLS))),
          by=.(origin_bg, cy, o_race, o_inc0, o_pov, o_bach, o_cbsa)]
  acc[[as.character(yr)]] <- ag
  cat("year",yr,"done\n")
}
oy <- rbindlist(acc); setnames(oy,"cy","year")
oy <- oy[!is.na(o_race) & o_inc0>0 & out_flow>0]
oy[, grp := relevel(factor(o_race), ref="Majority white (NH)")]
fwrite(oy, file.path(out_dir,"origin_year_multipollutant.csv"))
cat("origin-year rows:", nrow(oy), "\n")

gcoef <- function(fit,t) tryCatch(as.list(coeftable(fit)[t,c("Estimate","Std. Error","Pr(>|t|)")]), error=function(e) list(NA,NA,NA))

# ============================================================
# (A) AFFORDABILITY (Gelbach) per pollutant
# ============================================================
aff <- rbindlist(lapply(POLLS, function(p){
  d <- oy[!is.na(get(paste0("re_",p))) & !is.na(get(paste0("o_",p)))]
  f0 <- as.formula(paste0("re_",p," ~ grp + ns(o_",p,",4) + log(o_inc0) + o_pov + o_bach | o_cbsa + year"))
  f1 <- as.formula(paste0("re_",p," ~ grp + ns(o_",p,",4) + log(o_inc0) + o_pov + o_bach + log(dest_homeval) + log(dest_rent) | o_cbsa + year"))
  m0 <- feols(f0, d, weights=~out_flow, cluster=~o_cbsa)
  m1 <- feols(f1, d[dest_homeval>0 & dest_rent>0], weights=~out_flow, cluster=~o_cbsa)
  rbindlist(lapply(c("Majority Black","Majority Hispanic"), function(g){
    b0<-gcoef(m0,paste0("grp",g)); b1<-gcoef(m1,paste0("grp",g))
    data.table(pollutant=p, group=sub("Majority ","",g), gap_base=b0[[1]], se_base=b0[[2]],
               gap_housing=b1[[1]], pct_explained=100*(b0[[1]]-b1[[1]])/b0[[1]]) }))
}))
fwrite(aff, file.path(out_dir,"multipollutant_affordability.csv"))
cat("\n=== (A) AFFORDABILITY across pollutants (% of gap explained by destination housing values) ===\n")
print(aff[, .(pollutant, group, gap_base=round(gap_base,4), pct_explained=round(pct_explained,1))])

# ============================================================
# (B) SELECTION per pollutant
# ============================================================
sel <- rbindlist(lapply(POLLS, function(p) rbindlist(lapply(c("all","within_cbsa","short_dist","plus_dist"), function(s){
  yv <- switch(s, all=paste0("re_",p), within_cbsa=paste0("rewc_",p), short_dist=paste0("resd_",p), plus_dist=paste0("re_",p))
  wv <- switch(s, within_cbsa="out_flow_wc", short_dist="out_flow_sd", "out_flow")
  extra <- if (s=="plus_dist") " + mean_dist" else ""
  d <- oy[!is.na(get(yv)) & get(wv)>0 & !is.na(get(paste0("o_",p)))]
  f <- as.formula(paste0(yv," ~ grp + ns(o_",p,",4) + log(o_inc0) + o_pov + o_bach",extra," | o_cbsa + year"))
  m <- tryCatch(feols(f, d, weights=as.formula(paste0("~",wv)), cluster=~o_cbsa), error=function(e) NULL)
  if (is.null(m)) return(NULL)
  b <- gcoef(m,"grpMajority Black"); h <- gcoef(m,"grpMajority Hispanic")
  data.table(pollutant=p, spec=s, black=b[[1]], black_se=b[[2]], black_p=b[[3]], hisp=h[[1]], hisp_se=h[[2]], hisp_p=h[[3]])
}))))
fwrite(sel, file.path(out_dir,"multipollutant_selection.csv"))
cat("\n=== (B) SELECTION across pollutants (Black gap by subsample) ===\n")
print(dcast(sel, pollutant~spec, value.var="black")[, lapply(.SD, function(x) if(is.numeric(x)) round(x,4) else x)])

# ============================================================
# (C) OPPORTUNITY SET per pollutant (clean = below metro median of p)
# ============================================================
opp <- rbindlist(lapply(POLLS, function(p){
  a2 <- pan[!is.na(get(p)) & !is.na(med_home_value) & med_hh_income>0]
  a2[, opp_v := { thr<-median(get(p),na.rm=TRUE); ch<-sort(med_home_value[get(p)<thr])
                  if(!length(ch)) rep(0,.N) else findInterval(med_hh_income*AFFORD_MULT, ch)/.N }, by=.(cbsa,year)]
  m <- merge(oy[, .(origin_bg, year, o_race, out_flow)],
             a2[, .(origin_bg=geoid_bg, year, opp_v)], by=c("origin_bg","year"))
  m[, .(opp=wmean(opp_v,out_flow)), by=o_race][, pollutant:=p][]
}))
fwrite(opp, file.path(out_dir,"multipollutant_opportunity.csv"))
cat("\n=== (C) OPPORTUNITY SET across pollutants (share clean+affordable, by origin race/ethnicity) ===\n")
print(dcast(opp, pollutant~o_race, value.var="opp")[, lapply(.SD, function(x) if(is.numeric(x)) round(x,3) else x)])

# ============================================================
# (D) SEGREGATION interaction per pollutant
# ============================================================
acs[, `:=`(bpop=sh_black*total_pop, wpop=sh_white*total_pop)]
seg <- acs[cbsa %in% unique(oy$o_cbsa), .(D=0.5*sum(abs(bpop/sum(bpop,na.rm=TRUE)-wpop/sum(wpop,na.rm=TRUE)),na.rm=TRUE)), by=.(cbsa,year)][
  , .(seg=mean(D,na.rm=TRUE)), by=cbsa]
oy <- merge(oy, seg, by.x="o_cbsa", by.y="cbsa", all.x=TRUE); oy[, seg_z:=as.numeric(scale(seg))]
sgi <- rbindlist(lapply(POLLS, function(p){
  d <- oy[!is.na(get(paste0("re_",p))) & !is.na(seg_z) & !is.na(get(paste0("o_",p)))]
  f <- as.formula(paste0("re_",p," ~ grp*seg_z + ns(o_",p,",4) + log(o_inc0) + o_pov + o_bach | year"))
  m <- tryCatch(feols(f, d, weights=~out_flow, cluster=~o_cbsa), error=function(e) NULL)
  if (is.null(m)) return(NULL)
  b <- gcoef(m,"grpMajority Black:seg_z")
  data.table(pollutant=p, black_x_seg=b[[1]], se=b[[2]], p=b[[3]])
}))
fwrite(sgi, file.path(out_dir,"multipollutant_segregation.csv"))
cat("\n=== (D) SEGREGATION x Black gap interaction across pollutants ===\n")
print(sgi[, .(pollutant, black_x_seg=round(black_x_seg,4), se=round(se,4), p=signif(p,2))])

# ============================================================
# figures (one per analysis; pollutants ordered as in Fig 5)
# ============================================================
ordp <- intersect(c("no2","bc","om","pm25","no3","so4","ozone","wildfire_pm25"), POLLS)
plab <- function(d) d[, pl:=factor(PLAB[pollutant], levels=rev(unname(PLAB[ordp])))]
g1 <- ggplot(plab(copy(aff)), aes(pct_explained, pl, fill=group)) + geom_col(position="dodge") +
  geom_vline(xintercept=0, linewidth=.3) +
  scale_fill_manual(values=c(Black="#1b7837",Hispanic="#762a83"), name=NULL) +
  labs(title="Destination housing values explain little of the gap for any pollutant",
       subtitle="Gelbach share of the conditional avoidance gap explained by destination home value and rent",
       x="% of gap explained", y=NULL) + theme_minimal(base_size=10) +
  theme(plot.title=element_text(face="bold"), legend.position="bottom")
ggsave(file.path(fig_dir,"figS_mp_affordability.png"), g1, width=8, height=5, dpi=200)
s2 <- plab(copy(sel)); s2[, spec:=factor(spec, levels=c("all","within_cbsa","short_dist","plus_dist"))]
g2 <- ggplot(s2, aes(black, pl, color=spec)) + geom_vline(xintercept=0, linewidth=.3) +
  geom_pointrange(aes(xmin=black-1.96*black_se, xmax=black+1.96*black_se), position=position_dodge(.6), size=.3) +
  scale_color_brewer(palette="Dark2", name=NULL) +
  labs(title="The Black avoidance gap survives selection restrictions across pollutants",
       subtitle="Conditional gap (native units) by subsample; within-CBSA, short-distance (<=50 km), and +distance",
       x="Black avoidance gap", y=NULL) + theme_minimal(base_size=10) +
  theme(plot.title=element_text(face="bold"), legend.position="bottom")
ggsave(file.path(fig_dir,"figS_mp_selection.png"), g2, width=8.5, height=5.5, dpi=200)
o2 <- plab(copy(opp))[o_race!="No majority"]
g3 <- ggplot(o2, aes(opp, pl, fill=o_race)) + geom_col(position="dodge") +
  scale_fill_manual(values=c(`Majority Black`="#1b7837",`Majority Hispanic`="#762a83",`Majority white (NH)`="#2166ac"), name=NULL) +
  labs(title="The affordable-clean opportunity set is smaller for minority origins for every pollutant",
       subtitle="Clean = below metro-median of the pollutant; affordable = home value <= 3x origin income",
       x="Opportunity set (share of metro neighborhoods)", y=NULL) + theme_minimal(base_size=10) +
  theme(plot.title=element_text(face="bold"), legend.position="bottom")
ggsave(file.path(fig_dir,"figS_mp_opportunity.png"), g3, width=8.5, height=5.5, dpi=200)
g4 <- ggplot(plab(copy(sgi)), aes(black_x_seg, pl)) + geom_vline(xintercept=0, linewidth=.3) +
  geom_pointrange(aes(xmin=black_x_seg-1.96*se, xmax=black_x_seg+1.96*se), color="#1b7837") +
  labs(title="The Black gap widens with metropolitan segregation for combustion pollutants",
       subtitle="Interaction of the Black avoidance gap with standardized Black-white dissimilarity (native units)",
       x="Black gap x segregation (per SD of dissimilarity)", y=NULL) + theme_minimal(base_size=10) +
  theme(plot.title=element_text(face="bold"))
ggsave(file.path(fig_dir,"figS_mp_segregation.png"), g4, width=8, height=5, dpi=200)
cat("\nDONE. Multi-pollutant robustness outputs + 4 figures in ", out_dir, "\n", sep="")
