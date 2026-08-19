# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# MULTI-POLLUTANT DISPARITY TRENDS + SHIFT-SHARE DECOMPOSITION
# ------------------------------------------------------------
# Extends the two DESCRIPTIVE results (absolute/relative disparity trends,
# Figure 2; cleanup-vs-sorting shift-share, Figure 3) beyond PM2.5 to NO2
# and warm-season ozone, so the reader can see the aggregate dynamics
# generalize (or not) beyond the anchor pollutant. The mechanism panel
# (Figure 5) is already multi-pollutant; this fills in the descriptive side.
#
# Exposure is population-weighted BY RACIAL GROUP: for group g,
#   mean_g(t) = sum_bg [ poll_bg(t) * pop_g_bg(t) ] / sum_bg pop_g_bg(t),
#   pop_g_bg = share_g * total_pop   (the standard exposure-disparity mean).
# Shift-share (2011->2019), exact additive split of each group's change:
#   cleanup_g = sum w_g,2011 * (poll_2019 - poll_2011)   [pollution change, fixed distribution]
#   sorting_g = sum (w_g,2019 - w_g,2011) * poll_2019     [redistribution, fixed pollution]
#
# Output: pm25_avoidance_lag3_olag1/
#   multipollutant_disparity_trends.csv   + figures/figS_multipollutant_trends.png
#   multipollutant_decomposition.csv      + figures/figS_multipollutant_decomp.png
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr); library(ggplot2) })
base_dir <- migrate_data_dir()
out_dir  <- file.path(base_dir,"pm25_avoidance_lag3_olag1"); fig_dir <- file.path(out_dir,"figures")
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
YEARS <- 2011:2019
POLLS <- c("pm25","no2","ozone","bc","om","so4","no3","wildfire_pm25")   # policy pollutants + combustion components + wildfire (completeness)
PLAB  <- c(pm25="PM2.5", no2="NO2", ozone="Ozone", bc="Black carbon", om="Organic matter",
           so4="Sulfate", no3="Nitrate", wildfire_pm25="Wildfire PM2.5")
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")

# ---- panel + race + population ----
poll <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character"))
poll[, geoid_bg:=pad_bg(geoid_bg)]; poll <- unique(poll, by=c("geoid_bg","year"))
comp_file <- file.path(base_dir,"pollution_bg","bg_pm25_components_panel.csv")
if (file.exists(comp_file)) {
  cp <- fread(comp_file, colClasses=c(geoid_bg="character")); cp[, geoid_bg:=pad_bg(geoid_bg)]
  cc <- intersect(c("so4","no3","bc","om"), names(cp))
  cp <- unique(cp[, c("geoid_bg","year",cc), with=FALSE], by=c("geoid_bg","year"))
  poll <- merge(poll, cp, by=c("geoid_bg","year"), all.x=TRUE)
}
POLLS <- intersect(POLLS, names(poll))
rcs <- list.files(file.path(base_dir,"pollution_bg"), "^race_bg_.*csv$", full.names=TRUE)
rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}")); r<-fread(f,colClasses=c(geoid_bg="character")); r[,year:=y]; r }), fill=TRUE)
rc[, geoid_bg:=pad_bg(geoid_bg)]
pop <- fread(file.path(base_dir,"acs_gentrification","blockgroup_gentrification_2009_2019_lag3_with_HMDA.csv"),
             colClasses=c(geoid_bg="character"), select=c("year","geoid_bg","total_pop")); pop[, geoid_bg:=pad_bg(geoid_bg)]
D <- Reduce(function(a,b) merge(a,b,by=c("geoid_bg","year"),all=FALSE),
            list(poll[, c("geoid_bg","year",POLLS), with=FALSE],
                 rc[, .(geoid_bg,year,sh_black,sh_hisp,sh_white)], pop))
D[is.na(total_pop), total_pop:=0]
for (g in c("black","hisp","white")) D[[paste0("pop_",g)]] <- D[[paste0("sh_",g)]] * D$total_pop

wmean <- function(x,w){ ok<-!is.na(x)&!is.na(w)&w>0; if(!any(ok)) NA_real_ else sum(x[ok]*w[ok])/sum(w[ok]) }

# ---- disparity trends by pollutant ----
tr <- rbindlist(lapply(POLLS, function(p) rbindlist(lapply(YEARS, function(y){
  d <- D[year==y]
  b<-wmean(d[[p]],d$pop_black); h<-wmean(d[[p]],d$pop_hisp); w<-wmean(d[[p]],d$pop_white)
  data.table(pollutant=p, year=y, Black=b, Hispanic=h, White=w,
             gap_BW=b-w, gap_HW=h-w, ratio_BW=b/w, ratio_HW=h/w)
}))))
fwrite(tr, file.path(out_dir,"multipollutant_disparity_trends.csv"))
cat("\n=== MULTI-POLLUTANT disparity trends (population-weighted by group) ===\n")
print(tr[year %in% c(2011,2019), .(pollutant,year,Black=round(Black,2),White=round(White,2),ratio_BW=round(ratio_BW,3),ratio_HW=round(ratio_HW,3))])

# ---- shift-share decomposition 2011 -> 2019 ----
dec <- rbindlist(lapply(POLLS, function(p){
  a <- D[year==2011, .(geoid_bg, p0=get(p), pb0=pop_black, ph0=pop_hisp, pw0=pop_white)]
  b <- D[year==2019, .(geoid_bg, p1=get(p), pb1=pop_black, ph1=pop_hisp, pw1=pop_white)]
  m <- merge(a,b,by="geoid_bg")[!is.na(p0)&!is.na(p1)]
  one <- function(pg0,pg1){
    w0<-m[[pg0]]/sum(m[[pg0]],na.rm=TRUE); w1<-m[[pg1]]/sum(m[[pg1]],na.rm=TRUE)
    list(cleanup=sum(w0*(m$p1-m$p0),na.rm=TRUE), sorting=sum((w1-w0)*m$p1,na.rm=TRUE),
         d=sum(w1*m$p1,na.rm=TRUE)-sum(w0*m$p0,na.rm=TRUE))
  }
  B<-one("pb0","pb1"); H<-one("ph0","ph1"); W<-one("pw0","pw1")
  data.table(pollutant=p,
    gap=c("Black - White","Hispanic - White"),
    total_change=c(B$d-W$d, H$d-W$d),
    cleanup=c(B$cleanup-W$cleanup, H$cleanup-W$cleanup),
    sorting=c(B$sorting-W$sorting, H$sorting-W$sorting))
}))
fwrite(dec, file.path(out_dir,"multipollutant_decomposition.csv"))
cat("\n=== MULTI-POLLUTANT gap decomposition (change in the exposure gap: cleanup vs sorting) ===\n")
print(dec[, .(pollutant, gap, total=round(total_change,3), cleanup=round(cleanup,3), sorting=round(sorting,3))])

# ---- endpoint robustness: 2-year-averaged endpoints (2011-12 vs 2018-19) ----
avg_ep <- function(p){
  a <- D[year %in% 2011:2012, .(p0=mean(get(p),na.rm=TRUE), pb0=mean(pop_black,na.rm=TRUE), ph0=mean(pop_hisp,na.rm=TRUE), pw0=mean(pop_white,na.rm=TRUE)), by=geoid_bg]
  b <- D[year %in% 2018:2019, .(p1=mean(get(p),na.rm=TRUE), pb1=mean(pop_black,na.rm=TRUE), ph1=mean(pop_hisp,na.rm=TRUE), pw1=mean(pop_white,na.rm=TRUE)), by=geoid_bg]
  m <- merge(a,b,by="geoid_bg")[!is.na(p0)&!is.na(p1)]
  one <- function(pg0,pg1){ w0<-m[[pg0]]/sum(m[[pg0]],na.rm=TRUE); w1<-m[[pg1]]/sum(m[[pg1]],na.rm=TRUE)
    list(cleanup=sum(w0*(m$p1-m$p0),na.rm=TRUE), sorting=sum((w1-w0)*m$p1,na.rm=TRUE)) }
  B<-one("pb0","pb1"); W<-one("pw0","pw1")
  data.table(pollutant=p, cleanup_avg=B$cleanup-W$cleanup, sorting_avg=B$sorting-W$sorting)
}
dec_rob <- merge(dec[gap=="Black - White", .(pollutant, cleanup_endpoint=cleanup, sorting_endpoint=sorting)],
                 rbindlist(lapply(POLLS, avg_ep)), by="pollutant")
fwrite(dec_rob, file.path(out_dir,"decomposition_endpoint_robustness.csv"))
cat("\n=== ENDPOINT ROBUSTNESS: single-year (2011/2019) vs 2-year-averaged (2011-12/2018-19), Black-White gap ===\n")
print(dec_rob[, .(pollutant, cleanup_endpoint=round(cleanup_endpoint,3), cleanup_avg=round(cleanup_avg,3),
                  sorting_endpoint=round(sorting_endpoint,3), sorting_avg=round(sorting_avg,3))])

# ---- figures ----
trl <- melt(tr[, .(pollutant,year,ratio_BW,ratio_HW)], id.vars=c("pollutant","year"),
            variable.name="grp", value.name="ratio")
trl[, grp:=fifelse(grp=="ratio_BW","Black / White","Hispanic / White")]
trl[, plab:=factor(PLAB[pollutant], levels=unname(PLAB[POLLS]))]
g1 <- ggplot(trl, aes(year, ratio, color=grp)) + geom_hline(yintercept=1, linewidth=.3, color="grey60") +
  geom_line(linewidth=.9) + geom_point(size=1.2) + facet_wrap(~plab, scales="free_y") +
  scale_color_manual(values=c(`Black / White`="#1b7837",`Hispanic / White`="#762a83"), name=NULL) +
  labs(title="Relative exposure disparities across pollutants, 2011-2019",
       subtitle="Group-to-white ratio of population-weighted exposure (>1 = higher exposure than white). Wildfire smoke, episodic and non-capitalized, is unstable by design.",
       x=NULL, y="Exposure ratio (group / white)") +
  theme_minimal(base_size=10) + theme(plot.title=element_text(face="bold"), legend.position="bottom")
ggsave(file.path(fig_dir,"figS_multipollutant_trends.png"), g1, width=11, height=6.5, dpi=200)

decl <- melt(dec, id.vars=c("pollutant","gap"), measure.vars=c("cleanup","sorting"),
             variable.name="channel", value.name="val")
decl[, plab:=factor(PLAB[pollutant], levels=unname(PLAB[POLLS]))]
g2 <- ggplot(decl[gap=="Black - White"], aes(plab, val, fill=channel)) +
  geom_col(position="dodge") + geom_hline(yintercept=0, linewidth=.3) +
  scale_fill_manual(values=c(cleanup="#2166ac", sorting="#b2182b"),
                    labels=c("In-place cleanup","Population redistribution (sorting)"), name=NULL) +
  labs(title="What changed the Black-white exposure gap: cleanup vs migration, by pollutant",
       subtitle="Shift-share decomposition of the 2011-2019 change in the Black-white exposure gap",
       x=NULL, y="Contribution to change in the gap") +
  theme_minimal(base_size=10) + theme(plot.title=element_text(face="bold"), legend.position="bottom") +
  coord_flip()
ggsave(file.path(fig_dir,"figS_multipollutant_decomp.png"), g2, width=9, height=5, dpi=200)

cat("\nDONE. Multi-pollutant disparity + decomposition in ", out_dir, "\n", sep="")
