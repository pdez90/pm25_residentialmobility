# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# BUILD METRO-LEVEL DATA FOR THE SHINY "Metro explorer" TAB
# ------------------------------------------------------------
# For the N_METRO largest metropolitan areas (by within-metro moves), one pass
# over the MIGRATE flows restricted to moves that start and end in the same
# metro, pooled over 2011-2019. Writes three small tables into shiny_app/data/:
#
#   metro_list.csv      cbsa, name, lon, lat, moves, n_bg
#   metro_bg.csv        one row per block group (origin side): location, 2019
#                       PM2.5 and NO2, racial/ethnic majority, urban-form
#                       class, total out-moves, flow-weighted mean PM2.5
#                       reallocation of its out-movers (negative = cleaner)
#   metro_dest.csv      cbsa x destination BG x origin group: moves arriving
#                       and the mean PM2.5 those movers left behind
#   metro_sankey.csv    cbsa x origin group x destination PM2.5 quintile
#                       (quintiles within metro-year, Q5 = dirtiest): moves,
#                       mean reallocation
#
# Sizes: ~30 metros x ~3-6k block groups -> metro_bg ~100-150k rows (~10 MB),
# metro_dest ~4x that but only for destinations with >= 5 moves (~15 MB).
# Runtime: ~20 min (the MIGRATE pass). Run AFTER run_main.R.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(stringr); library(sf) })
sf::sf_use_s2(FALSE)
base_dir <- migrate_data_dir(); mig_dir <- base_dir
out_dir  <- file.path(base_dir,"pm25_avoidance_lag3_olag1")
app_dir  <- file.path(repo_root, "shiny_app", "data"); dir.create(app_dir, showWarnings=FALSE, recursive=TRUE)
migrate_years <- 2011:2019; N_METRO <- 30; MIN_DEST_MOVES <- 5
pad_bg <- function(x) str_pad(as.character(x),12,pad="0")
wmean  <- function(x,w){ok<-!is.na(x)&!is.na(w); if(!any(ok)) return(NA_real_); sum(x[ok]*w[ok])/sum(w[ok])}

# ---- block-group attributes ----
poll <- fread(file.path(base_dir,"pollution_bg","bg_pollution_panel.csv"), colClasses=c(geoid_bg="character"), select=c("geoid_bg","year","pm25","no2"))
poll[, geoid_bg:=pad_bg(geoid_bg)]; poll <- unique(poll, by=c("geoid_bg","year"))
rcs <- list.files(file.path(base_dir,"pollution_bg"), "^race_bg_.*csv$", full.names=TRUE)
rc <- rbindlist(lapply(rcs, function(f){ y<-as.integer(str_extract(basename(f),"[0-9]{4}")); r<-fread(f,colClasses=c(geoid_bg="character")); r[,year:=y]; r }), fill=TRUE)
rc[, geoid_bg:=pad_bg(geoid_bg)]
rc[, race_majority := fifelse(sh_black>0.5,"Majority Black", fifelse(sh_hisp>0.5,"Majority Hispanic", fifelse(sh_white>0.5,"Majority white (NH)","No majority")))]
cxw <- fread(file.path(base_dir,"pollution_bg","county_cbsa_xwalk.csv"), colClasses=c(county="character",cbsa="character"))
bg_cache <- file.path(base_dir,"pollution_bg","cache_multi","bg_national.rds")
if (!file.exists(bg_cache)) bg_cache <- file.path(base_dir,"pollution_bg","cache","bg_national.rds")
bg <- st_transform(readRDS(bg_cache), 4326); cc <- st_coordinates(st_centroid(bg))
cent <- data.table(geoid_bg=pad_bg(bg$geoid_bg), lon=round(cc[,1],4), lat=round(cc[,2],4)); rm(bg)
cent[, county := substr(geoid_bg,1,5)]; cent <- merge(cent, unique(cxw,by="county"), by="county", all.x=TRUE)
uf <- file.path(out_dir,"r1_bg_urban_form.csv")
if (file.exists(uf)) cent <- merge(cent, fread(uf, colClasses=c(geoid_bg="character"))[, .(geoid_bg=pad_bg(geoid_bg), uform)], by="geoid_bg", all.x=TRUE) else cent[, uform:=NA_character_]
# per-year attributes for the flow pass
att <- merge(poll, rc[, .(geoid_bg, year, race_majority)], by=c("geoid_bg","year"))
att <- merge(att, cent[, .(geoid_bg, cbsa)], by="geoid_bg")
att <- att[!is.na(cbsa) & !grepl("^nonmetro_", cbsa)]
setkey(att, year, geoid_bg)

# ---- metro names ----
cb <- tryCatch(tigris::core_based_statistical_areas(cb=TRUE, year=2019), error=function(e) NULL)
cbn <- if (!is.null(cb)) data.table(cbsa=cb$CBSAFP, name=cb$NAME, lsad=cb$LSAD) else data.table(cbsa=character(), name=character(), lsad=character())

# ---- MIGRATE pass: within-metro flows only ----
acc_o <- list(); acc_d <- list(); acc_s <- list(); acc_m <- list()
for (yr in migrate_years) {
  cat("  MIGRATE", yr, "...\n")
  x <- fread(file.path(mig_dir, paste0(yr,".csv")), select=c("Origin","Destination","Flow"))
  x <- x[, .(origin_bg=pad_bg(Origin), dest_bg=pad_bg(Destination), flow=as.numeric(Flow))][flow>0 & origin_bg!=dest_bg]
  x[, year := yr]
  x <- att[, .(year, origin_bg=geoid_bg, o_cbsa=cbsa, o_pm25=pm25, o_race=race_majority)][x, on=.(year, origin_bg), nomatch=NULL]
  x <- att[, .(year, dest_bg=geoid_bg, d_cbsa=cbsa, d_pm25=pm25)][x, on=.(year, dest_bg), nomatch=NULL]
  x <- x[o_cbsa==d_cbsa & !is.na(o_pm25) & !is.na(d_pm25)]
  x[, re := d_pm25 - o_pm25]
  # destination quintile within metro-year (over destination BGs reached, flow-weighted)
  x[, d_q := { q <- quantile(d_pm25, probs=seq(0,1,0.2), na.rm=TRUE); as.integer(cut(d_pm25, unique(q), include.lowest=TRUE, labels=FALSE)) }, by=.(o_cbsa)]
  x[, d_q := pmin(pmax(d_q, 1L), 5L)]
  acc_m[[as.character(yr)]] <- x[, .(moves=sum(flow)), by=.(cbsa=o_cbsa)]
  acc_o[[as.character(yr)]] <- x[, .(out_moves=sum(flow), re_sum=sum(flow*re)), by=.(cbsa=o_cbsa, geoid_bg=origin_bg)]
  acc_d[[as.character(yr)]] <- x[, .(in_moves=sum(flow), opm_sum=sum(flow*o_pm25)), by=.(cbsa=d_cbsa, dest_bg, o_race)]
  acc_s[[as.character(yr)]] <- x[, .(year=yr, moves=sum(flow), re_sum=sum(flow*re), dpm_sum=sum(flow*d_pm25), opm_sum=sum(flow*o_pm25)), by=.(cbsa=o_cbsa, o_race, dest_q=paste0("Q", d_q))]
  rm(x); gc()
}
mm <- rbindlist(acc_m)[, .(moves=sum(moves)), by=cbsa][order(-moves)]
top <- mm[cbsa %in% cbn[lsad=="M1"]$cbsa | nrow(cbn)==0][1:N_METRO]
cat("Top metros by within-metro moves:\n"); print(merge(top, cbn[, .(cbsa, name)], by="cbsa", all.x=TRUE)[order(-moves)])

# ---- assemble ----
o  <- rbindlist(acc_o)[cbsa %in% top$cbsa, .(out_moves=sum(out_moves), re_pm25=sum(re_sum)/sum(out_moves)), by=.(cbsa, geoid_bg)]
a19 <- att[year==2019, .(geoid_bg, pm25_2019=round(pm25,2), no2_2019=round(no2,2), race_majority)]
mbg <- merge(cent[cbsa %in% top$cbsa, .(geoid_bg, cbsa, lon, lat, uform)], a19, by="geoid_bg", all.x=TRUE)
mbg <- merge(mbg, o, by=c("cbsa","geoid_bg"), all.x=TRUE)
mbg[is.na(out_moves), out_moves := 0]; mbg[, `:=`(out_moves=round(out_moves,1), re_pm25=round(re_pm25,3))]
fwrite(mbg, file.path(app_dir,"metro_bg.csv"))

d  <- rbindlist(acc_d)[cbsa %in% top$cbsa, .(in_moves=sum(in_moves), origin_pm25=sum(opm_sum)/sum(in_moves)), by=.(cbsa, dest_bg, o_race)]
d  <- d[in_moves >= MIN_DEST_MOVES]
d  <- merge(d, mbg[, .(dest_bg=geoid_bg, lon, lat, dest_pm25=pm25_2019)], by="dest_bg")
d[, `:=`(in_moves=round(in_moves,1), origin_pm25=round(origin_pm25,2))]
fwrite(d, file.path(app_dir,"metro_dest.csv"))

sy <- rbindlist(acc_s)[cbsa %in% top$cbsa]
s  <- sy[, .(moves=round(sum(moves),1), re_pm25=round(sum(re_sum)/sum(moves),3)), by=.(cbsa, o_race, dest_q)]
fwrite(s[order(cbsa, o_race, dest_q)], file.path(app_dir,"metro_sankey.csv"))
# by year: quintile flows, and a per-year summary per origin group (for the trend chart)
fwrite(sy[, .(moves=round(sum(moves),1), re_pm25=round(sum(re_sum)/sum(moves),3)), by=.(cbsa, year, o_race, dest_q)][order(cbsa, year, o_race, dest_q)], file.path(app_dir,"metro_sankey_year.csv"))
my <- sy[, .(moves=round(sum(moves),1), re_pm25=round(sum(re_sum)/sum(moves),3), dest_pm25=round(sum(dpm_sum)/sum(moves),2), orig_pm25=round(sum(opm_sum)/sum(moves),2),
             share_dirtiest2=round(sum(moves[dest_q %in% c("Q4","Q5")])/sum(moves),3), share_cleanest2=round(sum(moves[dest_q %in% c("Q1","Q2")])/sum(moves),3)), by=.(cbsa, year, o_race)]
fwrite(my[order(cbsa, year, o_race)], file.path(app_dir,"metro_year.csv"))

ml <- merge(top, cbn[, .(cbsa, name)], by="cbsa", all.x=TRUE)
ml <- merge(ml, mbg[, .(lon=mean(lon), lat=mean(lat), n_bg=.N), by=cbsa], by="cbsa")
ml[is.na(name), name := cbsa]; ml[, moves := round(moves)]
fwrite(ml[order(-moves)], file.path(app_dir,"metro_list.csv"))
# principal-city marker for the map (from the Reviewer-1 classification, if run)
pcf <- file.path(out_dir,"r1_principal_cities.csv")
if (file.exists(pcf)) {
  pcs <- fread(pcf, colClasses=c(cbsa="character"))[cbsa %in% top$cbsa, .(cbsa, principal_city)]
  pl  <- tryCatch({ x <- do.call(rbind, lapply(unique(substr(mbg$geoid_bg,1,2)), function(st) tigris::places(state=st, cb=TRUE, year=2019)[,c("GEOID","NAME")]))
                    x <- sf::st_transform(x, 4326); xy <- sf::st_coordinates(sf::st_centroid(x)); data.table(place_geoid=x$GEOID, lon=round(xy[,1],4), lat=round(xy[,2],4)) }, error=function(e) NULL)
  if (!is.null(pl)) { pcs <- merge(fread(pcf, colClasses=c(cbsa="character",place_geoid="character"))[cbsa %in% top$cbsa], pl, by="place_geoid")
    fwrite(pcs[, .(cbsa, principal_city, lon, lat)], file.path(app_dir,"metro_principal_city.csv")) }
}
cat(sprintf("\nDONE. metro_bg %d rows, metro_dest %d rows, metro_sankey %d rows -> %s\n", nrow(mbg), nrow(d), nrow(s), app_dir))
for (f in c("metro_bg.csv","metro_dest.csv","metro_sankey.csv","metro_sankey_year.csv","metro_year.csv","metro_list.csv")) cat(sprintf("  %-18s %6.1f MB\n", f, file.size(file.path(app_dir,f))/1e6))
