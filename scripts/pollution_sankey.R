# ============================================================
# SANKEY / ALLUVIAL: migration flows across the pollution surface
# ------------------------------------------------------------
# For the SI. For each pollutant, shows how movers redistribute across
# pollution quintiles (origin quintile -> destination quintile,
# flow-weighted), and how origin racial majority maps to the
# destination pollution tercile reached. Makes the sorting result
# visible as flows.
#
# Output: pm25_avoidance_lag3_olag1/figures/
#   sankey_<pollutant>_quintile.png        (origin -> destination quintile)
#   sankey_<pollutant>_race_to_dest.png    (origin race -> dest tercile)
#   sankey_flows_<pollutant>.csv
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr); library(ggplot2) })
if (!requireNamespace("ggalluvial", quietly=TRUE)) install.packages("ggalluvial")
library(ggalluvial)
base_dir <- "/Users/priyanka/MIGRATE_dir"; mig_dir <- base_dir
gent_dir <- file.path(base_dir,"acs_gentrification"); change_lag <- 3
out_dir <- file.path(base_dir,"pm25_avoidance_lag3_olag1"); fig_dir <- file.path(out_dir,"figures")
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
migrate_years <- 2011:2019
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")
PLAB <- c(pm25="PM2.5", no2="NO2", wildfire_pm25="Wildfire PM2.5", ozone="Ozone",
          so4="Sulfate", no3="Nitrate", bc="Black carbon", om="Organic matter")

# ---- pollutant panel (+ components) ----
panel <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character")); panel[, geoid_bg:=pad_bg(geoid_bg)]
comp_file <- file.path(base_dir,"pollution_bg","bg_pm25_components_panel.csv")
if (file.exists(comp_file)) { cp<-fread(comp_file,colClasses=c(geoid_bg="character")); cp[,geoid_bg:=pad_bg(geoid_bg)]
  cc<-intersect(c("so4","no3","bc","om"),names(cp)); panel<-merge(panel,cp[,c("geoid_bg","year",cc),with=FALSE],by=c("geoid_bg","year"),all.x=TRUE) }
polls <- intersect(names(PLAB), names(panel))

# ---- origin race (majority) from cached pulls ----
rcs <- list.files(file.path(base_dir,"pollution_bg"), "^race_bg_.*csv$", full.names=TRUE)
rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}")); r<-fread(f,colClasses=c(geoid_bg="character")); r[,year:=y]; r }), fill=TRUE)
rc[, geoid_bg:=pad_bg(geoid_bg)]
rc[, race:=fifelse(sh_black>0.5,"Majority Black", fifelse(sh_hisp>0.5,"Majority Hispanic", fifelse(sh_white>0.5,"Majority white","No majority")))]
panel <- merge(panel, rc[, .(year,geoid_bg,race)], by=c("year","geoid_bg"), all.x=TRUE)
panel <- unique(panel, by=c("geoid_bg","year"))   # guard against duplicate BG-year rows before keyed flow joins

# ---- assign national quantiles per pollutant per year, by RANK ----
# rank-based bins are robust to ties (e.g. wildfire smoke is mostly ~0,
# so value-based quantile breakpoints collapse and cut() fails)
rankbin <- function(x, k) {
  r <- frank(x, ties.method="average", na.last="keep") / sum(!is.na(x))
  cut(r, breaks=seq(0,1,length.out=k+1), include.lowest=TRUE, labels=FALSE)
}
for (p in polls) {
  panel[, (paste0(p,"_q")) := factor(rankbin(get(p),5), levels=1:5, labels=paste0("Q",1:5)), by=year]
  panel[, (paste0(p,"_t")) := factor(rankbin(get(p),3), levels=1:3, labels=c("Cleanest","Middle","Dirtiest")), by=year]
}
qcols <- unlist(lapply(polls, function(p) c(paste0(p,"_q"), paste0(p,"_t"))))
oa <- panel[, c("year","geoid_bg","race", qcols), with=FALSE]
setnames(oa, c("year","geoid_bg","race", qcols), c("cy","origin_bg","o_race", paste0("o_",qcols))); setkey(oa,cy,origin_bg)
da <- panel[, c("year","geoid_bg", qcols), with=FALSE]
setnames(da, c("year","geoid_bg", qcols), c("cy","dest_bg", paste0("d_",qcols))); setkey(da,cy,dest_bg)

# ---- one MIGRATE pass: quintile-transition and race->dest-tercile flows ----
acc_q <- list(); acc_r <- list()
for (yr in migrate_years) {
  x <- fread(file.path(mig_dir, paste0(yr,".csv")), select=c("Origin","Destination","Flow"))
  x <- x[, .(origin_bg=pad_bg(Origin), dest_bg=pad_bg(Destination), flow=as.numeric(Flow))][flow>0 & origin_bg!=dest_bg]
  x[, cy:=yr]; setkey(x,cy,origin_bg); x<-oa[x,on=.(cy,origin_bg)]; setkey(x,cy,dest_bg); x<-da[x,on=.(cy,dest_bg)]
  for (p in polls) {
    oq<-paste0("o_",p,"_q"); dq<-paste0("d_",p,"_q")
    acc_q[[paste(p,yr)]] <- x[!is.na(get(oq)) & !is.na(get(dq)), .(pollutant=p, flow=sum(flow,na.rm=TRUE)), by=c(oq,dq)][
      , .(pollutant, origin_q=get(oq), dest_q=get(dq), flow)]
    dt<-paste0("d_",p,"_t")
    acc_r[[paste(p,yr)]] <- x[!is.na(o_race) & !is.na(get(dt)), .(pollutant=p, flow=sum(flow,na.rm=TRUE)), by=c("o_race",dt)][
      , .(pollutant, o_race, dest_t=get(dt), flow)]
  }
}
qflow <- rbindlist(acc_q)[, .(flow=sum(flow)), by=.(pollutant,origin_q,dest_q)]
rflow <- rbindlist(acc_r)[, .(flow=sum(flow)), by=.(pollutant,o_race,dest_t)]

# ---- draw one alluvial per pollutant ----
for (p in polls) {
  dq <- qflow[pollutant==p]; fwrite(dq, file.path(out_dir, paste0("sankey_flows_",p,".csv")))
  g1 <- ggplot(dq, aes(axis1=origin_q, axis2=dest_q, y=flow)) +
    geom_alluvium(aes(fill=origin_q), alpha=.7) + geom_stratum(width=.3, fill="grey90", color="grey40") +
    geom_text(stat="stratum", aes(label=after_stat(stratum)), size=3) +
    scale_x_discrete(limits=c("Origin","Destination"), expand=c(.05,.05)) +
    scale_fill_brewer(palette="RdYlGn", direction=-1, guide="none") +
    labs(title=paste0(PLAB[p], ": migration across pollution quintiles"),
         subtitle="Flow-weighted moves, origin quintile to destination quintile (Q5 = dirtiest)", y="Moves", x=NULL) +
    theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"))
  ggsave(file.path(fig_dir, paste0("sankey_",p,"_quintile.png")), g1, width=7.5, height=5, dpi=200)

  dr <- rflow[pollutant==p & o_race!="No majority"]
  g2 <- ggplot(dr, aes(axis1=o_race, axis2=dest_t, y=flow)) +
    geom_alluvium(aes(fill=o_race), alpha=.7) + geom_stratum(width=.3, fill="grey90", color="grey40") +
    geom_text(stat="stratum", aes(label=after_stat(stratum)), size=3) +
    scale_x_discrete(limits=c("Origin race","Destination"), expand=c(.05,.05)) +
    scale_fill_manual(values=c(`Majority Black`="#1b7837",`Majority Hispanic`="#762a83",`Majority white`="#2166ac"), guide="none") +
    labs(title=paste0(PLAB[p], ": origin race to destination pollution level"),
         subtitle="Flow-weighted moves from majority-race origins to destination pollution tercile", y="Moves", x=NULL) +
    theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"))
  ggsave(file.path(fig_dir, paste0("sankey_",p,"_race_to_dest.png")), g2, width=7.5, height=5, dpi=200)
}
message("\nDONE. Sankey diagrams (", length(polls), " pollutants x 2) in ", fig_dir)
