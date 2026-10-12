# ============================================================
# Air got cleaner. Moving didn't give everyone the same chance to escape pollution.
# Interactive results explorer for:
#   deSouza, P.N. "Residential Mobility and the Reproduction of Environmental-
#   Opportunity Inequality: Evidence from National Migration Flows and
#   Multi-Pollutant Exposure, 2011-2019" (Sustainability, in revision)
#
# Every number is read from the summary tables in data/, written by run_main.R,
# reviewer1_urban_form_distance.R and build_app_metro_data.R. Nothing is fetched
# at runtime. Deploy: see README.md.
# ============================================================
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(bsicons); library(data.table); library(ggplot2)
  library(plotly); library(DT); library(leaflet); library(viridisLite)
})

# ---------- data ----------
rd <- function(f) { p <- file.path("data", f); if (file.exists(p)) fread(p) else NULL }
D <- list(
  trends = rd("multipollutant_disparity_trends.csv"), decomp = rd("multipollutant_decomposition.csv"),
  expdec = rd("exposure_decomposition.csv"), cond = rd("conditional_models_all_fe.csv"),
  thresh = rd("race_threshold_sensitivity.csv"), select = rd("selection_robustness.csv"),
  meerr = rd("measurement_error_robustness.csv"), cross = rd("cross_pollutant_estimates.csv"),
  gelbach = rd("gelbach_decomposition.csv"), opp = rd("opportunity_by_race.csv"),
  oppsens = rd("opportunity_sensitivity.csv"), segq = rd("gap_by_segregation_quartile.csv"),
  segmp = rd("multipollutant_segregation.csv"), sub = rd("avoidance_by_subgroup.csv"),
  r1_band = rd("r1_gap_by_distance_band.csv"), r1_dest = rd("r1_black_gap_by_dest_type.csv"),
  r1_odfe = rd("r1_gap_od_type_fe.csv"), r1_mat = rd("r1_od_type_matrix.csv"), r1_dist = rd("r1_distance_distribution_by_group.csv"),
  mlist = rd("metro_list.csv"), mbg = rd("metro_bg.csv"), mdest = rd("metro_dest.csv"), msank = rd("metro_sankey.csv"), mpc = rd("metro_principal_city.csv"), mgap = rd("metro_gap.csv"), msy = rd("metro_sankey_year.csv"), myr = rd("metro_year.csv")
)
for (k in c("mlist","mbg","mdest","msank","mpc","mgap","msy","myr")) if (!is.null(D[[k]])) D[[k]][, cbsa := as.character(cbsa)]
if (!is.null(D$mbg)) D$mbg[, geoid_bg := as.character(geoid_bg)]
HAS_METRO <- !is.null(D$mbg) && !is.null(D$mlist) && !is.null(D$msank)
HAS_R1    <- !is.null(D$r1_band) && !is.null(D$r1_dest)
HAS_YEAR  <- !is.null(D$msy) && !is.null(D$myr)
# ---- national origin-destination tracing (optional; needs hosted Parquet + duckdb) ----
OD_BASE <- Sys.getenv("OD_PARQUET_BASE", "")
HAS_OD  <- nzchar(OD_BASE) && file.exists("data/bg_attr.parquet") && requireNamespace("duckdb", quietly=TRUE) && requireNamespace("arrow", quietly=TRUE)
if (HAS_OD) {
  BGA <- as.data.table(arrow::read_parquet("data/bg_attr.parquet")); BGA[, cbsa := as.character(cbsa)]
  od_con <- tryCatch({ con <- DBI::dbConnect(duckdb::duckdb()); DBI::dbExecute(con, "INSTALL httpfs; LOAD httpfs;"); con }, error=function(e) NULL)
  if (is.null(od_con)) HAS_OD <- FALSE
}
od_query <- function(origin, year) {
  q <- sprintf("SELECT Destination AS dest_bg, Flow AS flow FROM read_parquet('%s/od_%d.parquet') WHERE Origin = '%s'", OD_BASE, year, origin)
  as.data.table(DBI::dbGetQuery(od_con, q))
}
MAJ_LAB <- c(B="Majority Black", H="Majority Hispanic", W="Majority white (NH)", N="No majority")

POLL <- c("Total PM2.5"="pm25", "NO2"="no2", "Ozone (warm-season MDA8)"="ozone", "Wildfire smoke PM2.5"="wildfire_pm25",
          "Black carbon"="bc", "Organic matter"="om", "Sulfate"="so4", "Nitrate"="no3")
UNIT <- c(pm25="µg/m³", no2="ppb", ozone="ppb", wildfire_pm25="µg/m³", bc="µg/m³", om="µg/m³", so4="µg/m³", no3="µg/m³")
PGRP <- c(pm25="Traffic & combustion", no2="Traffic & combustion", bc="Traffic & combustion", om="Traffic & combustion", no3="Traffic & combustion",
          so4="Industrial / secondary", ozone="Secondary (ozone)", wildfire_pm25="Wildfire smoke")
GRP_COL <- c(Black="#B2403A", Hispanic="#D98E2B", White="#3B6FA0", "No majority"="#7A7A7A")
theme_app <- function() theme_minimal(base_size=13) + theme(panel.grid.minor=element_blank(), legend.position="bottom", plot.title.position="plot")
ci <- function(dt, ecol="est", scol="se") { dt[, lo := get(ecol) - 1.96*get(scol)]; dt[, hi := get(ecol) + 1.96*get(scol)]; dt }
fmt <- function(x, d=3) formatC(x, digits=d, format="f")
hleg <- function(p) p |> layout(legend=list(orientation="h", x=0, y=1.08, xanchor="left", yanchor="bottom"), margin=list(t=60)) |> config(displaylogo=FALSE, toImageButtonOptions=list(format="png", scale=2))

# ---------- headline numbers, read from the tables so they track the pipeline ----------
h <- list()
if (!is.null(D$cond)) { b <- D$cond[term=="grpMajority Black" & fe=="cbsa"]; h$gap <- b$Estimate; h$gap_lo <- b$Estimate-1.96*b$`Std. Error`; h$gap_hi <- b$Estimate+1.96*b$`Std. Error` }
if (!is.null(D$trends)) { t <- D$trends[pollutant=="pm25"]; h$b11 <- t[year==2011]$Black; h$b19 <- t[year==2019]$Black; h$w11 <- t[year==2011]$White; h$w19 <- t[year==2019]$White; h$r11 <- t[year==2011]$ratio_BW; h$r19 <- t[year==2019]$ratio_BW }
if (!is.null(D$decomp)) { d <- D$decomp[pollutant=="pm25" & gap=="Black - White"]; h$cleanup <- d$cleanup; h$sorting <- d$sorting }
if (!is.null(D$gelbach)) h$gelbach <- D$gelbach[group=="Black"]$pct_total
if (!is.null(D$opp)) { h$opp_b <- D$opp[o_race=="Majority Black"]$opp; h$opp_w <- D$opp[o_race=="Majority white (NH)"]$opp }
if (!is.null(D$sub)) { s <- D$sub[dimension=="o_race_majority"]; h$avoid_w <- s[subgroup=="Majority white (NH)"]$avoidance; h$orig_w <- s[subgroup=="Majority white (NH)"]$orig_pm25; h$orig_b <- s[subgroup=="Majority Black"]$orig_pm25 }
if (!is.null(D$segq)) { h$seg_lo <- D$segq[1]$black_gap; h$seg_hi <- D$segq[.N]$black_gap }

# ---------- reusable UI pieces ----------
how <- function(title, ...) accordion(open=FALSE, accordion_panel(title, icon=bs_icon("gear"), ...))
unitnote <- function(txt) p(class="text-muted small", bs_icon("info-circle"), " ", txt)
cite_block <- tagList(
  p(strong("Cite:"), " deSouza, P.N. (2026). Residential Mobility and the Reproduction of Environmental-Opportunity Inequality: Evidence from National Migration Flows and Multi-Pollutant Exposure, 2011–2019. ", em("Sustainability"), " (in revision)."),
  p("Code and data recipes: ", a("github.com/pdez90/pm25_residentialmobility", href="https://github.com/pdez90/pm25_residentialmobility", target="_blank"), ". Every chart has a camera icon (top right on hover) that downloads it as a PNG."))

# ---------- plot builders ----------
p_schematic <- function() {
  o <- if (length(h$orig_w)) round(mean(c(h$orig_w, h$orig_b)), 1) else 9
  aw <- if (length(h$avoid_w)) h$avoid_w else 0.05
  gap <- if (length(h$gap)) h$gap else 0.118
  d <- data.table(grp=c("Streams leaving\nmajority-white neighbourhoods","Streams leaving\nmajority-Black neighbourhoods"), x=c(1,2), origin=o, dest=c(o-aw, o-aw+gap))
  d[, grp := factor(grp, levels=grp)]
  ggplot(d) +
    geom_segment(aes(x=x, xend=x, y=origin, yend=dest, colour=grp), arrow=arrow(length=unit(0.3,"cm"), type="closed"), linewidth=1.6, show.legend=FALSE) +
    geom_point(aes(x=x, y=origin), size=5, colour="grey40") +
    geom_hline(yintercept=o, linetype=2, colour="grey60") +
    annotate("text", x=1.5, y=o+0.045, label=sprintf("Same starting pollution (%.1f µg/m³)", o), size=4, colour="grey30") +
    annotate("segment", x=2.35, xend=2.35, y=o-aw, yend=o-aw+gap, colour=GRP_COL[["Black"]], linewidth=0.8) +
    annotate("text", x=2.6, y=o-aw+gap/2, label=sprintf("the gap:\n%+.3f µg/m³", gap), size=4.2, colour=GRP_COL[["Black"]], fontface="bold") +
    scale_colour_manual(values=c(GRP_COL[["White"]], GRP_COL[["Black"]])) +
    scale_x_continuous(breaks=1:2, labels=levels(d$grp), limits=c(0.5, 2.9)) +
    scale_y_continuous(limits=c(o-aw-0.08, o+0.09)) +
    labs(x=NULL, y="PM2.5 at the average destination (µg/m³)", title="Same start, different finish",
         subtitle="Arrow length = average exposure change of movers, adjusted. Schematic, not observed levels.") +
    theme_app() + theme(axis.text.x=element_text(size=12))
}
p_trends <- function(poll, view) {
  t <- D$trends[pollutant==poll]; u <- UNIT[[poll]]
  if (view=="level") {
    m <- melt(t[, .(year, Black, Hispanic, White)], id.vars="year", variable.name="group", value.name="v")
    g <- ggplot(m, aes(year, v, colour=group)) + geom_line(linewidth=1.1) + geom_point(size=2) + scale_colour_manual(values=GRP_COL) + labs(x=NULL, y=sprintf("Population-weighted mean (%s)", u), colour=NULL)
  } else if (view=="ratio") {
    m <- melt(t[, .(year, `Black / White`=ratio_BW, `Hispanic / White`=ratio_HW)], id.vars="year", variable.name="group", value.name="v")
    g <- ggplot(m, aes(year, v, colour=group)) + geom_hline(yintercept=1, linetype=2, colour="grey50") + geom_line(linewidth=1.1) + geom_point(size=2) +
      scale_colour_manual(values=c(`Black / White`=GRP_COL[["Black"]], `Hispanic / White`=GRP_COL[["Hispanic"]])) + labs(x=NULL, y="Ratio of group means (1 = parity)", colour=NULL)
  } else {
    m <- melt(t[, .(year, `Black − White`=gap_BW, `Hispanic − White`=gap_HW)], id.vars="year", variable.name="group", value.name="v")
    g <- ggplot(m, aes(year, v, colour=group)) + geom_hline(yintercept=0, linetype=2, colour="grey50") + geom_line(linewidth=1.1) + geom_point(size=2) +
      scale_colour_manual(values=c(`Black − White`=GRP_COL[["Black"]], `Hispanic − White`=GRP_COL[["Hispanic"]])) + labs(x=NULL, y=sprintf("Gap (%s)", u), colour=NULL)
  }
  g + theme_app()
}
p_decomp <- function(poll) {
  d <- D$decomp[pollutant==poll]
  m <- melt(d[, .(gap, cleanup, sorting)], id.vars="gap", variable.name="term", value.name="v")
  m[, term := factor(term, levels=c("cleanup","sorting"), labels=c("In-place cleanup","Residential sorting"))]
  ggplot(m, aes(gap, v, fill=term)) + geom_col(position="stack") + geom_hline(yintercept=0, colour="grey40") +
    geom_point(data=d, aes(gap, total_change), inherit.aes=FALSE, shape=21, size=3.5, fill="white") +
    scale_fill_manual(values=c(`In-place cleanup`="#3B6FA0", `Residential sorting`="#D98E2B")) +
    labs(x=NULL, y=sprintf("Change in the exposure gap, 2011→2019 (%s); white dot = total", UNIT[[poll]]), fill=NULL) + theme_app()
}
p_cond <- function(fe_sel) {
  d <- copy(D$cond)[fe==fe_sel & dimension=="o_race_majority"]
  d[, group := sub("^grp", "", term)]; d[, `:=`(lo=Estimate-1.96*`Std. Error`, hi=Estimate+1.96*`Std. Error`)]
  d[, group := factor(group, levels=rev(c("Majority Black","Majority Hispanic","No majority")))]
  ggplot(d, aes(Estimate, group)) + geom_vline(xintercept=0, linetype=2, colour="grey50") +
    geom_errorbar(aes(xmin=lo, xmax=hi), width=0.2) + geom_point(size=3.5, colour="#B2403A") +
    labs(x="Extra PM2.5 at destination vs. movers from majority-white origins (µg/m³, 95% CI)", y=NULL) + theme_app()
}
p_cross <- function(grp_sel, sd=TRUE) {
  d <- copy(D$cross)[group==grp_sel]
  if (sd) { d[, `:=`(e=est_sd, s=se_sd)]; xl <- "Estimate, national SD units (comparable across pollutants)" } else { d[, `:=`(e=est, s=se)]; xl <- "Estimate, native units" }
  d[, `:=`(lo=e-1.96*s, hi=e+1.96*s, pgrp=PGRP[pollutant], sig=fifelse(p < 0.05, "p < 0.05", "not distinguishable from 0"))]
  ord <- c("no2","bc","om","pm25","no3","so4","ozone","wildfire_pm25")
  d[, pollutant_lab := factor(pollutant_lab, levels=rev(unique(D$cross[order(match(pollutant, ord))]$pollutant_lab)))]
  d[, pgrp := factor(pgrp, levels=c("Traffic & combustion","Industrial / secondary","Secondary (ozone)","Wildfire smoke"))]
  ggplot(d, aes(e, pollutant_lab, colour=sig)) + geom_vline(xintercept=0, linetype=2, colour="grey50") +
    geom_errorbar(aes(xmin=lo, xmax=hi), width=0.25) + geom_point(size=2.8) +
    facet_grid(pgrp ~ mechanism, scales="free", space="free_y") +
    scale_colour_manual(values=c(`p < 0.05`="#B2403A", `not distinguishable from 0`="grey55")) + labs(x=xl, y=NULL, colour=NULL) + theme_app() +
    theme(strip.text.y=element_text(angle=0, hjust=0))
}
p_gelbach <- function() {
  m <- melt(D$gelbach[, .(group, `Home value`=pct_homeval, `Gross rent`=pct_rent)], id.vars="group", variable.name="cov", value.name="pct")
  ggplot(m, aes(group, pct, fill=cov)) + geom_col() + scale_fill_manual(values=c(`Home value`="#3B6FA0", `Gross rent`="#D98E2B")) +
    scale_y_continuous(limits=c(0, 100)) + labs(x="Origin group", y="% of the gap explained (out of 100)", fill=NULL) + theme_app()
}
p_seg <- function() {
  d <- D$segq[, .(seg_q, est=black_gap, se)]; ci(d)
  ggplot(d, aes(seg_q, est)) + geom_hline(yintercept=0, linetype=2, colour="grey50") + geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2) + geom_point(size=3.5, colour="#B2403A") +
    labs(x="Metro segregation quartile (Black–white dissimilarity)", y="Black gap (µg/m³)") + theme_app()
}
p_opp <- function() {
  d <- copy(D$opp); d[, o_race := factor(o_race, levels=c("Majority white (NH)","No majority","Majority Hispanic","Majority Black"))]
  ggplot(d, aes(o_race, 100*opp, fill=o_race)) + geom_col(show.legend=FALSE) +
    scale_fill_manual(values=c(`Majority white (NH)`=GRP_COL[["White"]], `No majority`=GRP_COL[["No majority"]], `Majority Hispanic`=GRP_COL[["Hispanic"]], `Majority Black`=GRP_COL[["Black"]])) +
    labs(x="Origin neighbourhood", y="% of the metro's block groups that are\nclean AND affordable to that origin") + theme_app() + theme(legend.position="none")
}
p_trend <- function(cb) {
  y <- D$myr[cbsa==cb & o_race %in% c("Majority Black","Majority white (NH)","Majority Hispanic")]
  y[, o_race := factor(o_race, levels=c("Majority Black","Majority Hispanic","Majority white (NH)"), labels=c("Black-majority origins","Hispanic-majority origins","White-majority origins"))]
  m <- melt(y[, .(year, o_race, `% of moves to dirtiest two quintiles`=100*share_dirtiest2, `Mean Δ PM2.5 of movers (µg/m³)`=re_pm25)], id.vars=c("year","o_race"))
  ggplot(m, aes(year, value, colour=o_race)) + geom_line(linewidth=1) + geom_point(size=2) + facet_wrap(~variable, scales="free_y") +
    scale_colour_manual(values=c(`Black-majority origins`=GRP_COL[["Black"]], `Hispanic-majority origins`=GRP_COL[["Hispanic"]], `White-majority origins`=GRP_COL[["White"]])) +
    scale_x_continuous(breaks=2011:2019) + labs(x=NULL, y=NULL, colour=NULL) + theme_app()
}
p_r1 <- function() {
  gb <- copy(D$r1_band)[group=="Majority Black"]; gb[, band := factor(band, levels=unique(band))]
  gd <- copy(D$r1_dest); gd[, dest_type := factor(dest_type, levels=unique(dest_type))]
  pd <- rbind(gb[, .(panel="By distance moved", x=as.character(band), est, lo, hi)], gd[, .(panel="By destination type", x=as.character(dest_type), est, lo, hi)])
  pd[, x := factor(x, levels=c(levels(gb$band), levels(gd$dest_type)))]
  ggplot(pd, aes(x, est)) + geom_hline(yintercept=0, linetype=2, colour="grey50") + geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2) + geom_point(size=3, colour="#B2403A") +
    facet_wrap(~panel, scales="free_x") + labs(x=NULL, y="Black gap (µg/m³, 95% CI)") + theme_app() + theme(axis.text.x=element_text(angle=20, hjust=1))
}

sank_data <- function(cb, yr="all") { if (yr=="all" || is.null(D$msy)) D$msank[cbsa==cb] else D$msy[cbsa==cb & year==as.integer(yr)] }
p_quint <- function(cb, yr="all") {
  s <- sank_data(cb, yr)[o_race %in% c("Majority Black","Majority white (NH)","Majority Hispanic")]
  s[, share := moves/sum(moves), by=o_race]
  s[, o_race := factor(o_race, levels=c("Majority Black","Majority Hispanic","Majority white (NH)"), labels=c("Black-majority origins","Hispanic-majority origins","White-majority origins"))]
  s[, dest_q := factor(dest_q, levels=paste0("Q",1:5), labels=c("Q1 cleanest","Q2","Q3","Q4","Q5 dirtiest"))]
  ggplot(s, aes(dest_q, 100*share, fill=o_race, text=sprintf("%s → %s<br>%.0f%% of this group's moves<br>mean Δ PM2.5 %+.2f µg/m³", o_race, dest_q, 100*share, re_pm25))) +
    geom_col(position=position_dodge(width=0.75), width=0.7) +
    scale_fill_manual(values=c(`Black-majority origins`=GRP_COL[["Black"]], `Hispanic-majority origins`=GRP_COL[["Hispanic"]], `White-majority origins`=GRP_COL[["White"]])) +
    labs(x="Destination block group's PM2.5 quintile within this metro (2019)", y="% of the group's within-metro moves", fill=NULL) + theme_app()
}
metro_gap_text <- function(cb) {
  if (is.null(D$mgap)) return(NULL); g <- D$mgap[cbsa==cb & group=="Majority Black"]; if (!nrow(g)) return("Adjusted gap not estimable for this metro (too few majority-Black origins).")
  sig <- if (g$p < 0.05) "distinguishable from zero" else "not distinguishable from zero at the 5% level"
  sprintf("Adjusted gap in this metro: movers from majority-Black origins reached destinations %+.3f µg/m³ %s than comparable movers from majority-white origins (95%% CI %.3f to %.3f; %s; %s majority-Black origin-years). Same model as the national result, estimated within this metro only.",
          g$est, if (g$est >= 0) "more polluted" else "less polluted", g$lo, g$hi, sig, format(g$n_group, big.mark=","))
}

# ---------- metro helpers ----------
metro_insight <- function(cb, yr="all") {
  s <- sank_data(cb, yr); nm <- D$mlist[cbsa==cb]$name
  if (!nrow(s)) return(NULL)
  sh <- s[, .(dirty=sum(moves[dest_q %in% c("Q4","Q5")])/sum(moves), clean=sum(moves[dest_q %in% c("Q1","Q2")])/sum(moves), re=sum(moves*re_pm25)/sum(moves), n=sum(moves)), by=o_race]
  w <- sh[o_race=="Majority white (NH)"]; k <- sh[o_race=="Majority Black"]; hs <- sh[o_race=="Majority Hispanic"]
  parts <- sprintf("%s: %s within-metro moves, %s.", nm, format(round(sum(sh$n)), big.mark=","), if (yr=="all") "2011–2019" else yr)
  if (nrow(w) && nrow(k)) parts <- c(parts,
    sprintf("Of moves starting in majority-Black neighbourhoods, %.0f%% ended in the metro's two most polluted quintiles, versus %.0f%% of moves starting in majority-white neighbourhoods; %.0f%% and %.0f%% respectively reached the two cleanest.", 100*k$dirty, 100*w$dirty, 100*k$clean, 100*w$clean),
    sprintf("On average, movers from majority-white origins changed their PM2.5 by %+.2f µg/m³ and movers from majority-Black origins by %+.2f µg/m³ (raw averages, before adjusting for starting pollution).", w$re, k$re))
  if (nrow(hs)) parts <- c(parts, sprintf("Majority-Hispanic origins: %.0f%% to the two most polluted quintiles; %+.2f µg/m³ average change.", 100*hs$dirty, hs$re))
  paste(parts, collapse=" ")
}
base_map <- function(cb) {
  m <- leaflet(options=leafletOptions(preferCanvas=TRUE)) |> addProviderTiles(providers$Esri.WorldGrayCanvas, group="Light") |>
    addProviderTiles(providers$OpenStreetMap.Mapnik, group="Streets") |> addLayersControl(baseGroups=c("Light","Streets"), position="topright", options=layersControlOptions(collapsed=TRUE))
  if (!is.null(D$mpc)) { pc <- D$mpc[cbsa==cb]; if (nrow(pc)) m <- m |> addMarkers(data=pc, lng=~lon, lat=~lat, label=~principal_city, labelOptions=labelOptions(noHide=TRUE, direction="top")) }
  m
}
# pollution as a shaded background (every block-group centroid, low opacity), destinations as sized hollow circles
poll_layer <- function(m, cb) {
  b <- D$mbg[cbsa==cb & !is.na(pm25_2019)]; if (!nrow(b)) return(m)
  pal <- colorNumeric(c("#ffffd9","#c7e9b4","#41b6c4","#225ea8","#081d58"), domain=b$pm25_2019)
  m |> addCircleMarkers(data=b, lng=~lon, lat=~lat, radius=5, stroke=FALSE, fillOpacity=0.45, fillColor=~pal(pm25_2019), group="PM2.5 background",
         popup=~sprintf("<b>Block group %s</b><br>%s<br>PM2.5 2019: %.1f µg/m³<br>NO2 2019: %.1f ppb", geoid_bg, race_majority, pm25_2019, no2_2019)) |>
    addLegend("bottomleft", pal=pal, values=b$pm25_2019, title="PM2.5, 2019 (µg/m³)<br>shaded background", opacity=0.8)
}
dest_layer <- function(m, cb, grp, col) {
  d <- D$mdest[cbsa==cb & o_race==grp]; if (!nrow(d)) return(m)
  m |> addCircleMarkers(data=d, lng=~lon, lat=~lat, radius=~pmin(2 + sqrt(in_moves)/2.5, 16), stroke=TRUE, color=col, weight=1.5, fillColor=col, fillOpacity=0.25, group="Destinations",
         popup=~sprintf("<b>Destination block group %s</b><br>Arrivals from %s origins, 2011–19: %s<br>Destination PM2.5 (2019): %.1f µg/m³<br>Mean PM2.5 of those movers' origins: %.1f µg/m³", dest_bg, grp, format(round(in_moves)), dest_pm25, origin_pm25))
}
hl_layer <- function(m, cb, codes) {
  if (!length(codes)) return(m)
  lab <- c(B="Majority Black", W="Majority white (NH)", H="Majority Hispanic"); col <- c(B=GRP_COL[["Black"]], W=GRP_COL[["White"]], H=GRP_COL[["Hispanic"]])
  for (k in codes) { b <- D$mbg[cbsa==cb & race_majority==lab[[k]]]; if (nrow(b)) m <- m |> addCircleMarkers(data=b, lng=~lon, lat=~lat, radius=7, stroke=TRUE, color=col[[k]], weight=1.2, fill=FALSE, opacity=0.6, group=paste("Highlight", lab[[k]])) }
  m
}
GRP_OF <- c(black="Majority Black", white="Majority white (NH)", hisp="Majority Hispanic"); COL_OF <- c(black=GRP_COL[["Black"]], white=GRP_COL[["White"]], hisp=GRP_COL[["Hispanic"]])
dest_map <- function(cb, g, show_dest=TRUE, hl=character()) {
  m <- base_map(cb) |> poll_layer(cb) |> hl_layer(cb, hl)
  if (show_dest) m <- m |> dest_layer(cb, GRP_OF[[g]], COL_OF[[g]]) |> addControl(sprintf("<span style='background:%s;opacity:.6;display:inline-block;width:12px;height:12px;border-radius:6px'></span> circle size = arrivals from %s origins", COL_OF[[g]], tolower(sub(" \\(NH\\)","",GRP_OF[[g]]))), position="topleft")
  m
}

# ============================================================ UI
ui <- page_navbar(
  id = "main",
  title = "Mobility & environmental opportunity",
  window_title = "Air got cleaner. Moving didn't give everyone the same chance to escape pollution.",
  theme = bs_theme(version=5, bootswatch="flatly", primary="#3B6FA0", base_font=font_google("Source Sans 3"), heading_font=font_google("Source Serif 4")),
  fillable = FALSE,

  nav_panel("Start here",
    layout_columns(col_widths=c(7,5),
      card(class="border-0",
        h2("Air got cleaner. But moving didn't give everyone the same chance to escape pollution."),
        p(class="lead", "Between 2011 and 2019, air pollution fell across American neighbourhoods. Yet among neighbourhoods that started equally polluted, households moving out of majority-Black neighbourhoods tended to reach more polluted destinations than households moving out of majority-white ones."),
        p("This site walks through the evidence from a study of ", strong("every block-group-to-block-group move in the United States over nine years"), ", linked to maps of eight air pollutants. It is organised as six questions. Each page gives the answer in plain terms first and keeps the statistics behind ", em("How was this estimated?"), " panels."),
        tags$ol(class="fs-5",
          tags$li(actionLink("go1", "Did air get cleaner?"), " — yes, for everyone; the absolute gap narrowed, the relative gap did not."),
          tags$li(actionLink("go2", "Did moving help?"), " — the decline came from cleanup in place, not from people moving."),
          tags$li(actionLink("go3", "Who reached cleaner air?"), " — from the same starting point, movers from minority neighbourhoods captured less."),
          tags$li(actionLink("go4", "Why might this happen?"), " — not destination prices; it tracks segregation and traffic pollution."),
          tags$li(actionLink("go5", "Explore your city"), " — maps and flows for the 30 largest metros."),
          tags$li(actionLink("go6", "Methods & data"), " — how it was estimated, glossary, downloads.")),
        p(class="text-muted", "These are observational associations between neighbourhood types and the destinations their movers reached. The study documents the disparity and characterises it; it does not identify the mechanism that produces it, and it cannot see the race of any individual mover.")),
      card(card_header("The central comparison, in one picture"), plotlyOutput("schematic", height="430px"),
        unitnote("Illustrates the adjusted comparison in concept. The arrows use the study's average exposure change of movers and the estimated gap; they are not observed exposure levels for any neighbourhood."))),
    layout_columns(col_widths=c(3,3,3,3),
      value_box(title="Black-majority neighbourhoods, PM2.5", value=if (length(h$b11)) sprintf("%.1f → %.1f", h$b11, h$b19) else "—", p("µg/m³, 2011 → 2019"), theme="primary", showcase=bs_icon("arrow-down-right")),
      value_box(title="Relative gap to white-majority", value=if (length(h$r11)) sprintf("%.0f%% → %.0f%%", 100*(h$r11-1), 100*(h$r19-1)) else "—", p("above the white-majority level; did not close"), theme="secondary", showcase=bs_icon("dash-lg")),
      value_box(title="Extra pollution reached by movers from Black-majority origins", value=if (length(h$gap)) sprintf("+%s", fmt(h$gap)) else "—", p("µg/m³ vs. comparable white-majority origins"), theme="danger", showcase=bs_icon("signpost-split")),
      value_box(title="Clean & affordable neighbourhoods within reach", value=if (length(h$opp_b)) sprintf("%.0f%% vs %.0f%%", 100*h$opp_b, 100*h$opp_w) else "—", p("typical Black-majority vs white-majority origin — a price-based heuristic for reach, not a real choice set"), theme="info", showcase=bs_icon("map")))
  ),

  nav_panel("1 · Did air get cleaner?",
    layout_sidebar(sidebar=sidebar(width=300, selectInput("tr_poll", "Pollutant", POLL, selected="pm25"),
        radioButtons("tr_view", "Show", c("Exposure by neighbourhood type"="level", "Relative gap (ratio to white-majority)"="ratio", "Absolute gap"="gap"))),
      card(card_header("Yes — for everyone. The absolute gap narrowed; the relative gap did not."),
        p(if (length(h$b11)) sprintf("Residents of Black-majority neighbourhoods saw PM2.5 fall from %.1f to %.1f µg/m³, residents of white-majority neighbourhoods from %.1f to %.1f. The absolute Black–white difference shrank by about a fifth, but Black-majority neighbourhoods remained roughly %.0f%% more polluted at the end of the decade, as at the start. For NO2, a marker of traffic, the relative gap is larger still.", h$b11, h$b19, h$w11, h$w19, 100*(h$r19-1)) else ""),
        plotlyOutput("tr_plot", height="440px"),
        unitnote("Unit: population-weighted mean annual exposure of residents of block groups whose population is majority-Black, majority-Hispanic or majority-white (non-Hispanic), 2011–2019; ≈217,000 block groups per year."),
        how("How was this estimated?", p("Each pollutant surface (≈1 km grids; 10 km for smoke) is averaged to Census block groups by area weighting. Group means weight each block group by its population in that group and year, using ACS race/ethnicity tables matched to each year. Ratios and gaps are computed from those means. Source table: ", code("multipollutant_disparity_trends.csv"), "."))))
  ),

  nav_panel("2 · Did moving help?",
    layout_sidebar(sidebar=sidebar(width=300, selectInput("dc_poll", "Pollutant", POLL, selected="pm25")),
      card(card_header("Not much. Almost all of the improvement came from cleaner air in place, not from people moving."),
        p(if (length(h$cleanup)) sprintf("Of the %.3f µg/m³ narrowing in the Black–white PM2.5 gap, %.3f came from pollution falling where people already lived and only %.3f from people redistributing across neighbourhoods.", abs(h$cleanup+h$sorting), abs(h$cleanup), abs(h$sorting)) else ""),
        plotlyOutput("dc_gaps", height="400px"),
        unitnote("Unit: change in the gap between group means, 2011 → 2019, split into two terms that sum exactly to the total. Negative = gap narrowed."),
        card(class="bg-light", card_body(strong("But this does not mean moving was unimportant."), " This chart asks how much redistribution shifted a group's ", em("average"), " exposure. Some moves went to cleaner air and some to dirtier air, and those largely cancel when averaged over a group. It does not ask whether the households who moved fared equally well — that is the next question.")),
        how("How was this estimated?", p("A shift-share decomposition. Each group's exposure change is written as Σ w(g,i,2011)[P(i,2019) − P(i,2011)] (in-place cleanup: the 2011 residential distribution held fixed) plus Σ [w(g,i,2019) − w(g,i,2011)] P(i,2019) (sorting: pollution held at 2019 values), where w is the share of the group living in block group i. The two terms form an accounting identity, not a model. Source table: ", code("multipollutant_decomposition.csv"), "."))))
  ),

  nav_panel("3 · Who reached cleaner air?",
    card(card_header("From the same starting point, movers from minority neighbourhoods captured less of it."),
      layout_columns(col_widths=c(6,6),
        tagList(
          p(if (length(h$gap)) sprintf("Compare streams of movers leaving neighbourhoods with the same pollution, income, poverty and education, in the same metropolitan area and year. Movers from majority-Black neighbourhoods ended up at destinations %s µg/m³ more polluted (95%% CI %s to %s) than movers from majority-white neighbourhoods. That is small against total exposure (about 1.5%%), but it is more than double the average exposure improvement that moving delivers at all, and about one-seventh of the entire Black–white gap in 2019.", fmt(h$gap), fmt(h$gap_lo), fmt(h$gap_hi)) else ""),
          p("The raw comparison looks the opposite way: minority neighbourhoods start dirtier, so their movers record larger raw reductions. The gap appears only once starting pollution is held equal."),
          radioButtons("cg_fe", "Compare movers within the same…", c("Metropolitan area (main result)"="cbsa", "County"="county", "Nation (year only)"="none"), inline=TRUE),
          plotlyOutput("cg_plot", height="300px"),
          unitnote("Unit: origin block-group-years (≈1.9 million), each weighted by its number of out-moves, so these are statements about the average move. Confidence intervals use standard errors clustered on the metropolitan area."),
          card(class="bg-light", card_body(strong("What this contrast is — and is not."), " It compares the destinations reached by ", em("streams of movers leaving different kinds of neighbourhood"), ". It is not an estimate of any household's intention or ability to avoid pollution, of discrimination, or of the options open to a specific family. The study cannot see the race of individual movers; it sees where the out-flow from majority-Black neighbourhoods went compared with the out-flow from otherwise similar majority-white ones.")),
          if (!is.null(D$sub)) tagList(p(strong("The raw numbers behind it:")), DTOutput("cg_raw"), unitnote("Flow-weighted means over all moves 2011–2019, before adjustment. Minority origins start dirtier and record larger raw reductions; the gap appears once starting pollution is held equal."))),
        tagList(
          p(strong("Does it hold up?"), " Three checks, all giving the same answer:"),
          card(card_header("Different cut-offs for 'majority'"), plotlyOutput("cg_thresh", height="190px")),
          card(card_header("Only local moves (selection into moving)"), plotlyOutput("cg_select", height="190px")),
          card(card_header("Down-weighting uncertain pollution estimates"), plotlyOutput("cg_meerr", height="190px")))),
      how("How was this estimated?",
        p("For each origin block group and year, the exposure reallocation R(o,t) is the flow-weighted average over destinations of [pollution at destination − pollution at origin]. It is regressed on indicators for the origin's racial/ethnic majority, a natural spline (4 df) in origin pollution that absorbs mean reversion, log median income, poverty rate and share with a bachelor's degree, with metropolitan-area and year fixed effects, weighted by out-flow, standard errors clustered on metro. The coefficient on 'majority Black' is the gap shown. Source tables: ", code("conditional_models_all_fe.csv"), ", ", code("selection_robustness.csv"), ", ", code("race_threshold_sensitivity.csv"), ", ", code("measurement_error_robustness.csv"), "."),
        p(strong("Worked example of exposure reallocation."), " A block group at 10.0 µg/m³ sends 60 households to a destination at 9.0 and 40 to one at 10.5. R = (60×(9.0−10.0) + 40×(10.5−10.0)) / 100 = −0.40 µg/m³: its movers reached cleaner air by 0.40 on average."))),
    if (HAS_R1) card(card_header("Where in the metro does the gap arise, and how far do people move?"),
      layout_columns(col_widths=c(7,5),
        tagList(plotlyOutput("r1_plot", height="380px"),
          unitnote("Unit: origin block-group-years split by the type of destination reached (left) or by distance moved (right); separate models in each subset, same controls as above.")),
        tagList(
          p("Moves out of the urban core deliver the largest exposure gains for every group — but movers from majority-Black core neighbourhoods make those moves less often (66% of their within-metro moves stay in the core, versus 52% for majority-white). Comparing only movers making the same ", em("kind"), " of move (origin type × destination type) cuts the gap from 0.119 to 0.069 µg/m³: roughly 40% of the disparity is lower decentralisation; the rest is reaching more polluted suburbs when they do decentralise."),
          p("Movers from majority-Black neighbourhoods also move shorter distances (median 8 km vs 15 km). The gap is present for local moves under 10 km, where a new job is least likely to be the reason, and largest for moves over 200 km, which mostly cross metropolitan boundaries."),
          DTOutput("r1_odfe"))),
      how("How was this estimated?", p("Block groups are classified as urban core (inside the metro's principal city), inner suburb (same county as the principal city or within 15 km of its centre), periphery (rest of the metro), micropolitan or non-metropolitan. The main model is re-estimated with destination-type and origin×destination-type fixed effects, and separately within distance bands using great-circle distance between block-group centroids. Source tables: ", code("r1_*.csv"), "; script ", code("scripts/reviewer1_urban_form_distance.R"), ".")))
    else NULL
  ),

  nav_panel("4 · Why might this happen?",
    card(card_header("Not destination prices. The pattern tracks traffic pollution and metropolitan segregation."),
      p("The study cannot see why any household moved. What it can do is rule out one simple story and show which patterns are consistent with others.")),
    layout_columns(col_widths=c(6,6),
      card(card_header("Differences in destination home values and rents explain very little of the gap"),
        plotlyOutput("af_gelbach", height="300px"),
        p(if (length(h$gelbach)) sprintf("Adding the median home value and rent of the destinations reached accounts for %.2f%% of the Black gap. At the same destination price, movers from Black-majority origins still reached dirtier air.", h$gelbach) else ""),
        card(class="bg-light", card_body(strong("What this does not say."), " It does not rule out housing-market mechanisms. Discrimination and steering, rental availability, access to credit, search and information constraints, and differences within neighbourhoods can all shape where households end up without showing in median prices. It rules out ", em("destination price levels, as measured by Census medians"), " as a sufficient explanation — nothing broader.")),
        how("How was this estimated?", p("A Gelbach decomposition on the common sample of origin-years where both destination home value and rent are observed (home value is missing for 4.3% of block groups and rent for 15.8%, disproportionately rental neighbourhoods, so part of the rental market is not covered). The change in the group coefficient when the two housing variables are added is apportioned exactly between them. Descriptive, not causal mediation: destination price is chosen after the move. Source table: ", code("gelbach_decomposition.csv"), "."))),
      card(card_header("The gap grows with metropolitan segregation"),
        plotlyOutput("af_seg", height="300px"),
        p(if (length(h$seg_lo)) sprintf("In the least segregated quarter of metros the Black gap is %.3f µg/m³; in the most segregated it is %.3f.", h$seg_lo, h$seg_hi) else ""),
        p(class="text-muted small", "Consistent with accounts that emphasise discrimination and steering — and equally with co-ethnic preferences, information networks, and the spatial structure of housing supply in segregated markets. The data cannot separate these."),
        plotlyOutput("af_opp", height="280px"),
        unitnote("'Clean and affordable' = below the metro's median PM2.5 and median home value no more than 3× the origin's median income. A price-based heuristic for reach, not any household's real choice set; results are similar for 2.5×–5× income, rent-based rules and a 'cleanest quartile' definition."))),
    card(card_header("Three dimensions of environmental mobility inequality, across eight pollutants"),
      layout_sidebar(sidebar=sidebar(width=280, radioButtons("mp_group", "Origin group", c("Black","Hispanic")), checkboxInput("mp_sd", "Standardised units (compare across pollutants)", TRUE),
          p(class="small", strong("Reading the chart:"), " each row is one pollutant, grouped by type; columns are three separate questions. Red points are distinguishable from zero at p < 0.05.")),
        plotlyOutput("mp_plot", height="580px")),
      p("The disparity concentrates in pollutants tied to fixed urban infrastructure — NO2 from traffic above all, then black carbon, organic matter and total PM2.5. It is absent for wildfire smoke, which is regional and episodic and not durably priced into neighbourhoods, and it reverses for ozone, a secondary pollutant that forms downwind of its precursors and is typically higher in whiter exurban areas. That is the pattern a housing-market account predicts; it is also compatible with other accounts in which the spatial structure of a pollutant matters, so the study treats it as a consistency check rather than proof."),
      how("How was this estimated?", p("Three models per pollutant. Avoidance gap: the model from page 3. Sorting: destination pollution regressed on the share of its in-migrants who came from majority-Black or majority-Hispanic origins, with destination-county and year effects. Unequal cleanup: the 2011–2019 change in each block group's pollution regressed on its 2011 Black and Hispanic shares, a spline in 2011 pollution, and county effects. Standardised units divide by each pollutant's national standard deviation. Source table: ", code("cross_pollutant_estimates.csv"), ".")))
  ),

  nav_panel("5 · Explore your city",
    card(class="bg-light", card_body(h4("What we learned"), p("Cleaner air did not mean equal environmental opportunity. Pollution declined across U.S. neighbourhoods between 2011 and 2019, but the destinations reached through residential mobility remained unequal: movers leaving majority-Black neighbourhoods tended to reach more polluted destinations than those leaving otherwise comparable majority-white neighbourhoods. The disparity was most pronounced for traffic and combustion pollutants and in more segregated metropolitan areas, and destination housing prices alone explained little of it. The explorer below shows how this plays out in the 30 largest metros."))),
    if (!HAS_METRO) card(p("Metro-level maps are forthcoming."))
    else layout_sidebar(
      sidebar = sidebar(width=320,
        selectizeInput("me_cbsa", "Metropolitan area (type to search)", setNames(D$mlist$cbsa, D$mlist$name)),
        radioButtons("me_grp", "Movers from", c("Black-majority origins"="black", "White-majority origins"="white", "Hispanic-majority origins"="hisp")),
        if (HAS_YEAR) selectInput("me_year", "Year (charts and Sankey)", c("All years 2011–2019"="all", as.character(2011:2019))) else NULL,
        checkboxGroupInput("me_hl", "Highlight neighbourhoods on the map", c("Majority-Black"="B", "Majority-white"="W", "Majority-Hispanic"="H")),
        accordion(open=FALSE, accordion_panel("Advanced options", icon=bs_icon("sliders"),
          checkboxInput("me_showdest", "Show destinations on the map", TRUE),
          checkboxInput("me_compare", "Side-by-side comparison (Black- vs white-majority origins)", FALSE),
          checkboxInput("me_origins", "Origin view: colour each block group by the exposure change of its out-movers", FALSE),
          sliderInput("me_min", "Origin view: hide block groups with fewer out-moves than", 0, 200, 20, 10))),
        p(class="small text-muted", "Within-metro moves only, pooled 2011–2019. Dots are block-group centroids; click one for its numbers. The layer switcher at the top right of a map turns on street labels.")),
      card(card_header(textOutput("me_name")),
        plotlyOutput("me_quint", height="330px"),
        unitnote("Descriptive: share of each group's within-metro moves whose destination falls in each fifth of this metro's block groups ranked by 2019 PM2.5. Not adjusted for starting pollution.")),
      layout_columns(col_widths=c(6,6),
        card(card_header("What stands out here"), uiOutput("me_insight")),
        card(class="border-danger", card_header("Adjusted comparison, this metro only"), uiOutput("me_gap"))),
      if (HAS_YEAR) card(card_header("Year by year"), plotlyOutput("me_trend", height="320px"),
        unitnote("Each year's within-metro moves. Left: share of each group's moves ending in the metro's two most polluted quintiles. Right: average PM2.5 change experienced by movers (negative = reached cleaner air). Raw, not adjusted.")) else NULL,
      conditionalPanel("!input.me_compare",
        card(full_screen=TRUE, card_header(textOutput("me_maptitle")), leafletOutput("me_map", height="560px"),
          unitnote("Background shading: PM2.5 in 2019 for every block group (light = cleaner, dark = dirtier). Circles: destination block groups, sized by the number of arrivals from the selected origin group. Reading the map: are the circles sitting on light or dark ground?"))),
      conditionalPanel("input.me_compare",
        layout_columns(col_widths=c(6,6),
          card(card_header("Destinations of movers from Black-majority origins"), leafletOutput("me_mapA", height="480px")),
          card(card_header("Destinations of movers from white-majority origins"), leafletOutput("me_mapB", height="480px"))),
        unitnote("Same pollution background and the same circle-size scale on both maps.")),
      card(card_header(textOutput("me_sanktitle")), plotlyOutput("me_sankey", height="380px"),
        unitnote("Same data as the bar chart, drawn as flows. Hover a ribbon for moves and the mean change in PM2.5 those movers experienced.")),
      card(card_header("This metro in numbers"), DTOutput("me_table")))
  ),

  nav_panel("Trace a neighbourhood",
    if (!HAS_OD) card(card_header("Trace a neighbourhood"), p("Coming soon: click any block group in the country, pick a year, and see where its movers went."))
    else layout_sidebar(
      sidebar = sidebar(width=320,
        selectizeInput("od_cbsa", "Metropolitan area", setNames(D$mlist$cbsa, D$mlist$name)),
        selectInput("od_year", "Year", as.character(2011:2019), selected="2019"),
        checkboxGroupInput("od_hl", "Highlight neighbourhoods", c("Majority-Black"="B", "Majority-white"="W", "Majority-Hispanic"="H")),
        checkboxInput("od_showopp", "After clicking: outline the block groups that were clean AND affordable from the selected origin", TRUE),
        p(class="small text-muted", "Click any block group on the map to see where the households that left it in the selected year went. Origins are coloured by PM2.5 in that year; destinations appear as circles sized by movers, with lines to the origin (the 400 largest destinations are drawn; the charts and numbers use all of them).")),
      card(full_screen=TRUE, card_header(textOutput("od_title")), leafletOutput("od_map", height="600px")),
      layout_columns(col_widths=c(5,7),
        card(card_header("This neighbourhood's movers"), uiOutput("od_summary")),
        card(card_header("Reachable vs. reached"), plotlyOutput("od_reach", height="340px"),
          unitnote("Each curve: the share of neighbourhoods (or of movers) at or below a given PM2.5. Grey: every block group in the metro. Blue: the ones affordable to this origin (median home value ≤ 3× its median income — the paper's rule of thumb). Red: where its movers actually went, weighted by movers. Dashed lines: the origin's own PM2.5 and the metro median that defines 'clean'. A red curve to the right of the blue one means movers landed in dirtier places than the affordable set offered."))),
      layout_columns(col_widths=c(5,7),
        card(card_header("Clean and affordable from here"), uiOutput("od_oppset"),
          unitnote("'Clean' = below the metro's median PM2.5 in that year; 'affordable' = median home value within 3× this origin's median income. A price-based heuristic for reach, not any household's real choice set.")),
        card(card_header("Where they went: destination PM2.5"), plotlyOutput("od_hist", height="300px"),
          unitnote("Distribution of movers by their destination's PM2.5 in the same year, with the origin's own PM2.5 marked. Left of the line = reached cleaner air."))))
  ),

  nav_panel("6 · Methods & data",
    layout_columns(col_widths=c(5,7),
      tagList(
        card(card_header("Design in one paragraph"),
          p("National block-group origin–destination migration flows (MIGRATE, 2011–2019; ≈217,000 origins per year) are linked to eight pollutant surfaces aggregated to block groups. For each origin and year the exposure reallocation is the flow-weighted difference between destination and origin pollution. The conditional avoidance gap compares that quantity across origin neighbourhood types holding starting pollution, socioeconomic characteristics, metropolitan area and year fixed. A shift-share decomposition separates in-place cleanup from sorting; a Gelbach decomposition asks how much destination prices explain; cross-pollutant, segregation, opportunity-set, urban-form and distance analyses characterise the pattern; and ACS microdata provide a coarser individual-level check."),
          p(strong("What the data can and cannot say."), " MIGRATE flows are not disaggregated by race, so every contrast is about streams leaving majority-Black, majority-Hispanic or majority-white ", em("neighbourhoods"), ", not about individual movers. Pollution surfaces at ≈1 km smooth near-road gradients, which if anything attenuates the contrasts. Results describe 2011–2019; wildfire smoke, remote work and housing markets have changed since."),
          cite_block),
        card(card_header("Glossary"),
          tags$dl(
            tags$dt("Exposure reallocation"), tags$dd("For the households leaving a neighbourhood in a year, the flow-weighted average of (pollution at destination − pollution at origin). Negative = moved to cleaner air."),
            tags$dt("Conditional avoidance gap"), tags$dd("The difference in exposure reallocation between streams leaving minority-majority and white-majority neighbourhoods that started equally polluted, with the same income, poverty and education, in the same metro and year."),
            tags$dt("Shift-share decomposition"), tags$dd("Splits a change in a group's average exposure into the part due to pollution changing where people already lived and the part due to people changing where they live. The two sum exactly to the total."),
            tags$dt("Gelbach decomposition"), tags$dd("Apportions the change in a coefficient, when extra variables are added to a model, exactly among those variables, independent of the order they are added."),
            tags$dt("Opportunity set"), tags$dd("Share of a metro's block groups that are both cleaner than the metro median and affordable to a given origin (home value ≤ 3× origin median income). A price-based proxy for reach."),
            tags$dt("Dissimilarity index"), tags$dd("Share of one group that would have to move for two groups to be evenly distributed across a metro's block groups; 0 = integrated, 1 = complete segregation."),
            tags$dt("Block group"), tags$dd("Census unit of roughly 600–3,000 people; ≈217,000 in the study. Migration PUMA: county-based unit of ≥100,000 people used for the microdata check.")))),
      tagList(
        card(card_header("Sources"),
          tags$ul(tags$li("Migration flows: MIGRATE (Agostini et al., Nat. Commun. 2025)"), tags$li("PM2.5 and components, with uncertainty: WashU SatPM V6.NA.01 (Shen et al., 2026)"),
            tags$li("NO2 (~1 km): Anenberg et al. 2022"), tags$li("Ozone (warm-season MDA8, 0.01°): Liu et al. 2026, Yale Dataverse"), tags$li("Wildfire smoke PM2.5 (10 km daily): Childs et al. 2022"),
            tags$li("Neighbourhood covariates: ACS 5-year via Census API and NHGIS; HMDA"), tags$li("Individual movers: IPUMS USA, ACS 2012–2019"))),
        card(card_header("Browse and download the output tables"),
          selectInput("tb_file", "Table", choices=sort(list.files("data", pattern="\\.csv$")), width="70%"),
          DTOutput("tb_dt"), downloadButton("tb_dl", "Download CSV"))))
  )
)

# ============================================================ SERVER
server <- function(input, output, session) {
  observeEvent(input$go1, nav_select("main", "1 · Did air get cleaner?")); observeEvent(input$go2, nav_select("main", "2 · Did moving help?"))
  observeEvent(input$go3, nav_select("main", "3 · Who reached cleaner air?")); observeEvent(input$go4, nav_select("main", "4 · Why might this happen?"))
  observeEvent(input$go5, nav_select("main", "5 · Explore your city")); observeEvent(input$go6, nav_select("main", "6 · Methods & data"))

  output$schematic <- renderPlotly({ ggplotly(p_schematic(), tooltip="none") |> hleg() |> layout(showlegend=FALSE) })
  output$tr_plot <- renderPlotly({ req(D$trends); ggplotly(p_trends(input$tr_poll, input$tr_view)) |> hleg() })
  output$dc_gaps <- renderPlotly({ req(D$decomp); ggplotly(p_decomp(input$dc_poll)) |> hleg() })
  output$cg_plot <- renderPlotly({ req(D$cond); ggplotly(p_cond(input$cg_fe)) |> hleg() })
  small <- function(d) { ci(d); ggplot(d, aes(est, spec)) + geom_vline(xintercept=0, linetype=2, colour="grey50") + geom_errorbar(aes(xmin=lo, xmax=hi), width=0.2) + geom_point(size=2.5, colour="#B2403A") + labs(x="Black gap (µg/m³)", y=NULL) + theme_app() }
  output$cg_thresh <- renderPlotly({ req(D$thresh); ggplotly(small(D$thresh[, .(spec, est=black_gap, se=black_se)])) |> hleg() })
  output$cg_select <- renderPlotly({ req(D$select); d <- D$select[group=="Black"]; d[, spec := factor(spec, levels=c("all","within_cbsa","short_dist","plus_dist"), labels=c("All moves","Within metro","≤ 50 km","+ distance control"))]; ggplotly(small(d)) |> hleg() })
  output$cg_raw <- renderDT({ req(D$sub); d <- D$sub[dimension=="o_race_majority", .(`Origin neighbourhood`=subgroup, `Mean origin PM2.5`=round(orig_pm25,2), `Mean change at destination`=round(mean_change,3), `Moves (millions)`=round(flow/1e6,1))]; datatable(d, options=list(dom="t"), rownames=FALSE) })
  output$cg_meerr <- renderPlotly({ req(D$meerr); ggplotly(small(D$meerr[, .(spec, est=Black, se=Black_se)])) |> hleg() })
  output$r1_plot <- renderPlotly({ req(HAS_R1); ggplotly(p_r1()) |> hleg() })
  output$r1_odfe <- renderDT({ req(D$r1_odfe); datatable(D$r1_odfe[group=="Majority Black", .(Specification=spec, `Black gap`=round(est,3), `95% CI`=sprintf("%.3f to %.3f", lo, hi))], options=list(dom="t"), rownames=FALSE) })
  output$mp_plot <- renderPlotly({ req(D$cross); ggplotly(p_cross(input$mp_group, input$mp_sd)) |> hleg() })
  output$af_gelbach <- renderPlotly({ req(D$gelbach); ggplotly(p_gelbach()) |> hleg() })
  output$af_seg <- renderPlotly({ req(D$segq); ggplotly(p_seg()) |> hleg() })
  output$af_opp <- renderPlotly({ req(D$opp); ggplotly(p_opp()) |> hleg() })

  # metro explorer
  output$me_name <- renderText({ req(HAS_METRO, input$me_cbsa); D$mlist[cbsa==input$me_cbsa]$name })
  me_yr <- reactive({ if (HAS_YEAR && !is.null(input$me_year)) input$me_year else "all" })
  output$me_insight <- renderUI({ req(HAS_METRO, input$me_cbsa); p(metro_insight(input$me_cbsa, me_yr())) })
  output$me_trend <- renderPlotly({ req(HAS_YEAR, input$me_cbsa); ggplotly(p_trend(input$me_cbsa)) |> hleg() })
  output$me_sanktitle <- renderText({ req(input$me_cbsa); sprintf("Flows from origin group to destination pollution quintile, %s", if (me_yr()=="all") "all years" else me_yr()) })
  output$me_gap <- renderUI({ req(HAS_METRO, input$me_cbsa); t <- metro_gap_text(input$me_cbsa); if (is.null(t)) p(class="text-muted small", "Metro-specific adjusted gaps not yet built (scripts/build_app_metro_gaps.R).") else p(strong("Adjusted comparison. "), t) })
  output$me_quint <- renderPlotly({ req(HAS_METRO, input$me_cbsa); ggplotly(p_quint(input$me_cbsa, me_yr()), tooltip="text") |> hleg() })
  output$me_maptitle <- renderText({ req(input$me_cbsa); if (isTRUE(input$me_origins)) "Origin view: average PM2.5 change experienced by households leaving each block group" else sprintf("Where movers from %s went, on this metro's pollution map", c(black="Black-majority origins", white="white-majority origins", hisp="Hispanic-majority origins")[[input$me_grp]]) })
  output$me_map <- renderLeaflet({
    req(HAS_METRO, input$me_cbsa)
    if (!isTRUE(input$me_origins)) return(dest_map(input$me_cbsa, input$me_grp, isTRUE(input$me_showdest), input$me_hl))
    b <- D$mbg[cbsa==input$me_cbsa & out_moves >= input$me_min]; req(nrow(b) > 0)
    lim <- quantile(abs(b$re_pm25), 0.95, na.rm=TRUE); b[, v := pmax(pmin(re_pm25, lim), -lim)]
    pal <- colorNumeric(c("#2166AC","#F7F7F7","#B2182B"), domain=c(-lim, lim))
    base_map(input$me_cbsa) |> hl_layer(input$me_cbsa, input$me_hl) |> addCircleMarkers(data=b, lng=~lon, lat=~lat, radius=~pmin(2 + sqrt(out_moves)/6, 9), stroke=FALSE, fillOpacity=0.75, fillColor=~pal(v),
        popup=~sprintf("<b>Block group %s</b><br>%s<br>Out-moves 2011–19: %s<br>Mean PM2.5 change of out-movers: %+.2f µg/m³<br>PM2.5 2019: %.1f µg/m³", geoid_bg, race_majority, format(round(out_moves)), re_pm25, pm25_2019)) |>
      addLegend("bottomright", pal=pal, values=c(-lim, lim), title="Δ PM2.5 of out-movers<br>(blue = reached cleaner air)", opacity=0.9)
  })
  output$me_mapA <- renderLeaflet({ req(HAS_METRO, input$me_cbsa, isTRUE(input$me_compare)); dest_map(input$me_cbsa, "black", TRUE, input$me_hl) })
  output$me_mapB <- renderLeaflet({ req(HAS_METRO, input$me_cbsa, isTRUE(input$me_compare)); dest_map(input$me_cbsa, "white", TRUE, input$me_hl) })
  output$me_sankey <- renderPlotly({
    req(HAS_METRO, input$me_cbsa); s <- copy(sank_data(input$me_cbsa, me_yr()))[o_race %in% c("Majority white (NH)","Majority Black","Majority Hispanic","No majority")]; req(nrow(s) > 0)
    og <- c("Majority white (NH)","Majority Black","Majority Hispanic","No majority"); dq <- paste0("Q", 1:5)
    nodes <- c(og, paste0(dq, c(" (cleanest)","","","", " (dirtiest)")))
    s[, `:=`(src = match(o_race, og) - 1L, tgt = length(og) + match(dest_q, dq) - 1L)]; s[, share := moves / sum(moves), by=o_race]
    lc <- c(`Majority white (NH)`="rgba(59,111,160,0.45)", `Majority Black`="rgba(178,64,58,0.45)", `Majority Hispanic`="rgba(217,142,43,0.45)", `No majority`="rgba(122,122,122,0.45)")
    plot_ly(type="sankey", orientation="h", arrangement="snap",
      node=list(label=nodes, pad=18, thickness=18, color=c(GRP_COL[["White"]], GRP_COL[["Black"]], GRP_COL[["Hispanic"]], GRP_COL[["No majority"]], rev(viridis(5))), line=list(width=0.5, color="grey")),
      link=list(source=s$src, target=s$tgt, value=s$moves, color=lc[s$o_race],
                customdata=sprintf("%s → %s<br>%s moves (%.0f%% of this group's within-metro moves)<br>mean Δ PM2.5: %+.2f µg/m³", s$o_race, s$dest_q, format(round(s$moves), big.mark=","), 100*s$share, s$re_pm25),
                hovertemplate="%{customdata}<extra></extra>")) |> layout(font=list(size=12), margin=list(l=10, r=10, t=10, b=10)) |> config(displaylogo=FALSE)
  })
  output$me_table <- renderDT({
    req(HAS_METRO, input$me_cbsa); b <- D$mbg[cbsa==input$me_cbsa & out_moves > 0]
    t <- b[, .(`Block groups`=.N, `Out-moves`=format(round(sum(out_moves)), big.mark=","), `Mean PM2.5 2019`=round(weighted.mean(pm25_2019, out_moves, na.rm=TRUE), 2),
              `Mean Δ PM2.5 of out-movers`=round(weighted.mean(re_pm25, out_moves, na.rm=TRUE), 3)), by=.(`Origin majority`=race_majority)][order(`Origin majority`)]
    datatable(t, options=list(dom="t"), rownames=FALSE)
  })

  # trace a neighbourhood (national OD, remote Parquet)
  od_sel <- reactiveVal(NULL)
  od_bgs <- reactive({ req(HAS_OD, input$od_cbsa); BGA[cbsa==input$od_cbsa] })
  output$od_title <- renderText({ req(HAS_OD); if (is.null(od_sel())) "Click a block group to trace its movers" else sprintf("Movers leaving block group %s in %s", od_sel(), input$od_year) })
  output$od_map <- renderLeaflet({
    req(HAS_OD); b <- copy(od_bgs()); yr <- input$od_year; pmc <- paste0("pm25_", yr); mjc <- paste0("maj_", yr)
    b[, pm := get(pmc)]; b[, mj := MAJ_LAB[get(mjc)]]; b <- b[!is.na(pm)]
    pal <- colorNumeric(c("#ffffd9","#c7e9b4","#41b6c4","#225ea8","#081d58"), domain=b$pm)
    m <- base_map(input$od_cbsa) |> addCircleMarkers(data=b, lng=~lon, lat=~lat, layerId=~geoid_bg, radius=5, stroke=FALSE, fillOpacity=0.55, fillColor=~pal(pm),
            label=~sprintf("%s · %s · PM2.5 %.1f", geoid_bg, mj, pm)) |>
      addLegend("bottomleft", pal=pal, values=b$pm, title=sprintf("PM2.5, %s (µg/m³)", yr), opacity=0.8)
    for (k in input$od_hl) { hb <- b[get(mjc)==k]; if (nrow(hb)) m <- m |> addCircleMarkers(data=hb, lng=~lon, lat=~lat, radius=7, stroke=TRUE, color=c(B=GRP_COL[["Black"]], W=GRP_COL[["White"]], H=GRP_COL[["Hispanic"]])[[k]], weight=1.2, fill=FALSE, opacity=0.6) }
    m
  })
  observeEvent(input$od_map_marker_click, { id <- input$od_map_marker_click$id; if (!is.null(id) && grepl("^[0-9]{12}$", id)) od_sel(id) })
  od_flows <- reactive({
    req(HAS_OD, od_sel()); yr <- as.integer(input$od_year)
    f <- tryCatch(od_query(od_sel(), yr), error=function(e) NULL); req(!is.null(f))
    f <- merge(f, BGA[, .(dest_bg=geoid_bg, lon, lat, cbsa, pm=get(paste0("pm25_", yr)), mj=MAJ_LAB[get(paste0("maj_", yr))])], by="dest_bg", all.x=TRUE)
    f[order(-flow)]
  })
  observe({
    f <- od_flows(); o <- BGA[geoid_bg==od_sel()]; yr <- input$od_year; req(nrow(o))
    proxy <- leafletProxy("od_map") |> clearGroup("trace")
    if (!nrow(f)) return()
    fw <- f[!is.na(lon)][order(-flow)][1:min(.N, 400)]        # draw the 400 largest destinations; the summary and charts use all of them
    for (i in seq_len(min(nrow(fw), 150))) proxy <- proxy |> addPolylines(lng=c(o$lon, fw$lon[i]), lat=c(o$lat, fw$lat[i]), weight=pmax(0.5, log1p(fw$flow[i])), color="#444", opacity=0.3, group="trace")
    proxy |> addCircleMarkers(data=fw, lng=~lon, lat=~lat, radius=~pmin(2 + 2*sqrt(flow), 14), stroke=TRUE, color="#B2403A", weight=1.5, fillColor="#B2403A", fillOpacity=0.3, group="trace",
               popup=~sprintf("<b>Destination %s</b><br>%s<br>Movers: %.1f<br>PM2.5 %s: %.1f µg/m³", dest_bg, mj, flow, yr, pm)) |>
      addCircleMarkers(lng=o$lon, lat=o$lat, radius=10, stroke=TRUE, color="black", weight=3, fill=FALSE, group="trace")
  })
  od_reach <- reactive({
    req(HAS_OD, od_sel()); yr <- input$od_year; o <- BGA[geoid_bg==od_sel()]; req(nrow(o))
    inc <- o[[paste0("inc_", yr)]]; metro <- BGA[cbsa==o$cbsa & !is.na(get(paste0("pm25_", yr)))]
    metro[, `:=`(pm=get(paste0("pm25_", yr)), hv=get(paste0("hv_", yr)))]
    med <- median(metro$pm, na.rm=TRUE)
    metro[, affordable := !is.na(hv) & !is.na(inc) & hv <= 3*inc]; metro[, clean := pm < med]
    list(metro=metro[, .(geoid_bg, lon, lat, pm, affordable, clean)], med=med, inc=inc, opm=o[[paste0("pm25_", yr)]])
  })
  observe({ r <- od_reach(); proxy <- leafletProxy("od_map") |> clearGroup("oppset"); if (isTRUE(input$od_showopp)) { ca <- r$metro[affordable & clean]; if (nrow(ca)) proxy |> addCircleMarkers(data=ca, lng=~lon, lat=~lat, radius=6, stroke=TRUE, color="#2E8B57", weight=1.5, fill=FALSE, opacity=0.8, group="oppset") } })
  output$od_reach <- renderPlotly({
    r <- od_reach(); f <- od_flows()[!is.na(pm)]; req(nrow(f))
    ecdf_w <- function(x, w) { o <- order(x); x <- x[o]; w <- w[o]; data.table(pm=x, share=cumsum(w)/sum(w)) }
    d <- rbind(cbind(set="All neighbourhoods in the metro", ecdf_w(r$metro$pm, rep(1, nrow(r$metro)))),
               cbind(set="Affordable to this origin", ecdf_w(r$metro[affordable==TRUE]$pm, rep(1, sum(r$metro$affordable)))),
               cbind(set="Where its movers went", ecdf_w(f$pm, f$flow)))
    d[, set := factor(set, levels=c("All neighbourhoods in the metro","Affordable to this origin","Where its movers went"))]
    g <- ggplot(d, aes(pm, 100*share, colour=set)) + geom_step(linewidth=1.1) +
      geom_vline(xintercept=r$opm, linetype=2, colour="grey30") + geom_vline(xintercept=r$med, linetype=3, colour="#2E8B57") +
      annotate("text", x=r$opm, y=5, label="origin", hjust=-0.1, colour="grey30", size=3.5) + annotate("text", x=r$med, y=95, label="metro median ('clean')", hjust=-0.05, colour="#2E8B57", size=3.5) +
      scale_colour_manual(values=c(`All neighbourhoods in the metro`="grey55", `Affordable to this origin`="#3B6FA0", `Where its movers went`="#B2403A")) +
      labs(x=sprintf("PM2.5, %s (µg/m³)", input$od_year), y="% at or below this level", colour=NULL) + theme_app()
    ggplotly(g) |> hleg()
  })
  output$od_oppset <- renderUI({
    r <- od_reach(); f <- od_flows(); req(nrow(f))
    m <- r$metro; n_all <- nrow(m); n_aff <- sum(m$affordable); n_ca <- sum(m$affordable & m$clean)
    fk <- merge(f, m[, .(dest_bg=geoid_bg, affordable, clean)], by="dest_bg", all.x=TRUE)
    s_ca <- sum(fk$flow[fk$affordable & fk$clean], na.rm=TRUE)/sum(fk$flow); s_cl <- sum(fk$flow[fk$clean], na.rm=TRUE)/sum(fk$flow); s_af <- sum(fk$flow[fk$affordable], na.rm=TRUE)/sum(fk$flow)
    tagList(
      p(if (is.na(r$inc)) "Origin median income not available for this year, so the affordable set cannot be computed." else sprintf("Origin median household income: $%s. Of the metro's %s block groups, %s (%.0f%%) had a median home value within 3× that income, and %s (%.0f%%) were both affordable and cleaner than the metro median.",
          format(r$inc, big.mark=","), format(n_all, big.mark=","), format(n_aff, big.mark=","), 100*n_aff/n_all, format(n_ca, big.mark=","), 100*n_ca/n_all)),
      tags$ul(tags$li(sprintf("%.0f%% of its movers went somewhere cleaner than the metro median.", 100*s_cl)),
              tags$li(sprintf("%.0f%% went somewhere affordable by this rule; %.0f%% went somewhere both clean and affordable.", 100*s_af, 100*s_ca)),
              tags$li(if (n_ca/n_all > 0 && s_ca < n_ca/n_all) sprintf("Clean-and-affordable neighbourhoods were %.0f%% of the metro but received %.0f%% of these movers.", 100*n_ca/n_all, 100*s_ca) else sprintf("Clean-and-affordable neighbourhoods were %.0f%% of the metro and received %.0f%% of these movers.", 100*n_ca/n_all, 100*s_ca))),
      p(class="text-muted small", "Green outlines on the map mark the clean-and-affordable set."))
  })
  output$od_summary <- renderUI({
    req(HAS_OD, od_sel()); f <- od_flows(); o <- BGA[geoid_bg==od_sel()]; yr <- input$od_year
    opm <- o[[paste0("pm25_", yr)]]; omj <- MAJ_LAB[o[[paste0("maj_", yr)]]]
    if (!nrow(f)) return(p("No flows of one or more expected movers recorded for this block group in this year."))
    tot <- sum(f$flow); fk <- f[!is.na(pm)]; dpm <- weighted.mean(fk$pm, fk$flow); within <- sum(f$flow[f$cbsa==o$cbsa], na.rm=TRUE)/tot
    cleaner <- sum(fk$flow[fk$pm < opm])/sum(fk$flow); same <- sum(fk$flow[fk$mj==omj], na.rm=TRUE)/sum(fk$flow)
    tagList(
      p(strong(sprintf("Block group %s", od_sel())), sprintf(" — %s, PM2.5 %.1f µg/m³ in %s.", omj, opm, yr)),
      tags$ul(
        tags$li(sprintf("%s expected movers left in %s, spread over %s destination block groups (the largest %d account for %.0f%%).", format(round(tot), big.mark=","), yr, format(nrow(f), big.mark=","), min(nrow(f), 50), 100*sum(f$flow[1:min(nrow(f),50)])/tot)),
        tags$li(sprintf("Their destinations averaged %.1f µg/m³ — %s than the origin by %.2f.", dpm, if (dpm < opm) "cleaner" else "dirtier", abs(dpm-opm))),
        tags$li(sprintf("%.0f%% reached cleaner air than they left; %.0f%% stayed within the same metro; %.0f%% went to another %s neighbourhood.", 100*cleaner, 100*within, 100*same, tolower(sub(" \\(NH\\)", "", omj))))),
      p(class="text-muted small", "Expected movers are MIGRATE's model-based estimates, spread over many destinations and often fractional; totals are sums of those expectations."))
  })
  output$od_hist <- renderPlotly({
    req(HAS_OD, od_sel()); f <- od_flows()[!is.na(pm)]; o <- BGA[geoid_bg==od_sel()]; opm <- o[[paste0("pm25_", input$od_year)]]; req(nrow(f))
    g <- ggplot(f, aes(pm, weight=flow)) + geom_histogram(bins=30, fill="#3B6FA0", alpha=0.85) + geom_vline(xintercept=opm, colour="#B2403A", linewidth=1.2, linetype=2) +
      annotate("text", x=opm, y=Inf, label="origin", vjust=1.5, hjust=-0.1, colour="#B2403A") + labs(x=sprintf("Destination PM2.5, %s (µg/m³)", input$od_year), y="Movers") + theme_app()
    ggplotly(g) |> hleg()
  })

  tb <- reactive({ req(input$tb_file); fread(file.path("data", input$tb_file)) })
  output$tb_dt <- renderDT(datatable(tb(), options=list(pageLength=15, scrollX=TRUE), rownames=FALSE) |> formatRound(names(tb())[sapply(tb(), is.numeric)], 4))
  output$tb_dl <- downloadHandler(filename=function() input$tb_file, content=function(f) fwrite(tb(), f))
}

shinyApp(ui, server)
