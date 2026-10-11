# ============================================================
# Residential mobility and environmental-opportunity inequality
# Interactive results explorer (Shiny) for the paper:
#   deSouza, P.N. "Residential Mobility and the Reproduction of
#   Environmental-Opportunity Inequality through Housing Markets:
#   evidence from national migration flows and multi-pollutant
#   exposure, 2011-2019"
#
# All numbers come from the summary tables written by run_main.R
# (data/ next to this file; copied from pm25_avoidance_lag3_olag1/).
# Deploy: see README.md in this folder.
# ============================================================
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(bsicons); library(data.table); library(ggplot2); library(plotly); library(DT); library(leaflet)
})

# ---------- data ----------
rd <- function(f) { p <- file.path("data", f); if (file.exists(p)) fread(p) else NULL }
D <- list(
  trends   = rd("multipollutant_disparity_trends.csv"),
  decomp   = rd("multipollutant_decomposition.csv"),
  expdec   = rd("exposure_decomposition.csv"),
  cond     = rd("conditional_models_all_fe.csv"),
  thresh   = rd("race_threshold_sensitivity.csv"),
  select   = rd("selection_robustness.csv"),
  meerr    = rd("measurement_error_robustness.csv"),
  cross    = rd("cross_pollutant_estimates.csv"),
  gelbach  = rd("gelbach_decomposition.csv"),
  opp      = rd("opportunity_by_race.csv"),
  oppsens  = rd("opportunity_sensitivity.csv"),
  oppmp    = rd("multipollutant_opportunity.csv"),
  segq     = rd("gap_by_segregation_quartile.csv"),
  segmp    = rd("multipollutant_segregation.csv"),
  sortq    = rd("sorting_by_pm25_quintile.csv"),
  corr     = rd("pollutant_correlations.csv"),
  r1_band  = rd("r1_gap_by_distance_band.csv"),
  r1_dest  = rd("r1_black_gap_by_dest_type.csv"),
  r1_odfe  = rd("r1_gap_od_type_fe.csv"),
  r1_mat   = rd("r1_od_type_matrix.csv"),
  r1_dist  = rd("r1_distance_distribution_by_group.csv"),
  mlist    = rd("metro_list.csv"),
  mbg      = rd("metro_bg.csv"),
  mdest    = rd("metro_dest.csv"),
  msank    = rd("metro_sankey.csv")
)
if (!is.null(D$mbg)) { D$mbg[, geoid_bg := as.character(geoid_bg)]; D$mbg[, cbsa := as.character(cbsa)] }
if (!is.null(D$mdest)) { D$mdest[, cbsa := as.character(cbsa)] }
if (!is.null(D$msank)) { D$msank[, cbsa := as.character(cbsa)] }
if (!is.null(D$mlist)) { D$mlist[, cbsa := as.character(cbsa)] }
HAS_METRO <- !is.null(D$mbg) && !is.null(D$mlist) && !is.null(D$msank)
POLL <- c("Total PM2.5"="pm25", "NO2"="no2", "Ozone (warm-season MDA8)"="ozone", "Wildfire smoke PM2.5"="wildfire_pm25",
          "Black carbon"="bc", "Organic matter"="om", "Sulfate"="so4", "Nitrate"="no3")
UNIT <- c(pm25="µg/m³", no2="ppb", ozone="ppb", wildfire_pm25="µg/m³", bc="µg/m³", om="µg/m³", so4="µg/m³", no3="µg/m³")
plab <- function(p) names(POLL)[match(p, POLL)]
GRP_COL <- c(Black="#B2403A", Hispanic="#D98E2B", White="#3B6FA0", "No majority"="#7A7A7A", Gentrifying="#5B9A5B")
theme_app <- function() theme_minimal(base_size=13) + theme(panel.grid.minor=element_blank(), legend.position="bottom", plot.title.position="plot")
ci <- function(dt, ecol="est", scol="se") { dt[, lo := get(ecol) - 1.96*get(scol)]; dt[, hi := get(ecol) + 1.96*get(scol)]; dt }
fmt <- function(x, d=3) formatC(x, digits=d, format="f")

# headline numbers (from the tables, so they stay in sync with the pipeline)
h <- list()
if (!is.null(D$cond)) { b <- D$cond[term=="grpMajority Black" & fe=="cbsa"]; h$black_gap <- b$Estimate; h$black_lo <- b$Estimate-1.96*b$`Std. Error`; h$black_hi <- b$Estimate+1.96*b$`Std. Error` }
if (!is.null(D$trends)) { t <- D$trends[pollutant=="pm25"]; h$b11 <- t[year==2011]$Black; h$b19 <- t[year==2019]$Black; h$w11 <- t[year==2011]$White; h$w19 <- t[year==2019]$White; h$r11 <- t[year==2011]$ratio_BW; h$r19 <- t[year==2019]$ratio_BW }
if (!is.null(D$decomp)) { d <- D$decomp[pollutant=="pm25" & gap=="Black - White"]; h$cleanup <- d$cleanup; h$sorting <- d$sorting; h$total <- d$total_change }
if (!is.null(D$gelbach)) h$gelbach <- D$gelbach[group=="Black"]$pct_total
if (!is.null(D$opp)) { h$opp_b <- D$opp[o_race=="Majority Black"]$opp; h$opp_w <- D$opp[o_race=="Majority white (NH)"]$opp }
if (!is.null(D$cross)) { n <- D$cross[mechanism=="Conditional avoidance gap" & group=="Black" & pollutant=="no2"]; h$no2_gap <- n$est }

# ---------- UI ----------
ui <- page_navbar(
  title = "Mobility & environmental opportunity, 2011–2019",
  theme = bs_theme(version=5, bootswatch="flatly", primary="#3B6FA0", base_font=font_google("Source Sans 3"), heading_font=font_google("Source Serif 4")),
  fillable = FALSE,

  nav_panel("Overview",
    layout_columns(col_widths=c(8,4),
      card(card_header("What the study asks"),
        p("Air pollution in the United States fell sharply between 2011 and 2019. This study asks whether ",
          strong("moving house"), " helped close the racial and ethnic gap in exposure, or reproduced it. It links ",
          "national block-group migration flows (MIGRATE) to eight pollutant surfaces and measures, for the households leaving each neighbourhood, ",
          "the flow-weighted difference between the pollution at their destinations and at their origin — the ", em("exposure reallocation"), "."),
        p("Three results carry the paper:"),
        tags$ol(
          tags$li(strong("Cleanup, not migration, drove the decline."), " In-place reductions account for essentially all of each group's exposure change; population redistribution contributed only hundredths of a µg/m³."),
          tags$li(strong("Mobility distributed those gains unequally."), " Conditional on starting pollution, streams leaving majority-Black neighbourhoods reached destinations ",
                  if (!is.null(h$black_gap)) sprintf("%s µg/m³ dirtier (95%% CI %s–%s) ", fmt(h$black_gap), fmt(h$black_lo), fmt(h$black_hi)) else "measurably dirtier ",
                  "than streams leaving comparable majority-white neighbourhoods."),
          tags$li(strong("Destination prices do not explain it."), " Home values and rents at the destination account for ",
                  if (!is.null(h$gelbach)) sprintf("%.2f%%", h$gelbach) else "well under 1%", " of the gap; the gap is largest for traffic and combustion pollutants, absent for wildfire smoke, and reversed for ozone.")
        ),
        p(class="text-muted", "These are observational associations at the neighbourhood-stream level. The analysis documents a disparity and characterises it; it does not identify the mechanism that produces it. Use the tabs above to explore each result; every chart is drawn from the output tables of the public analysis code.")
      ),
      layout_columns(col_widths=12,
        value_box(title="Black-majority PM2.5, 2011 → 2019", value=if (!is.null(h$b11)) sprintf("%.1f → %.1f µg/m³", h$b11, h$b19) else "—",
                  showcase=bsicons::bs_icon("arrow-down-right"), theme="primary",
                  p(if (!is.null(h$r11)) sprintf("Relative to white-majority: %.0f%% above in 2011, %.0f%% above in 2019", 100*(h$r11-1), 100*(h$r19-1)) else "")),
        value_box(title="Conditional avoidance gap (Black vs white origins)", value=if (!is.null(h$black_gap)) sprintf("+%s µg/m³", fmt(h$black_gap)) else "—",
                  showcase=bsicons::bs_icon("signpost-split"), theme="danger", p("Dirtier destinations reached, net of origin pollution, SES and metro")),
        value_box(title="Share of gap explained by destination price", value=if (!is.null(h$gelbach)) sprintf("%.2f%%", h$gelbach) else "—",
                  showcase=bsicons::bs_icon("house"), theme="secondary", p("Gelbach decomposition: home value + gross rent")),
        value_box(title="Clean & affordable neighbourhoods in reach", value=if (!is.null(h$opp_b)) sprintf("%.0f%% vs %.0f%%", 100*h$opp_b, 100*h$opp_w) else "—",
                  showcase=bsicons::bs_icon("map"), theme="info", p("Typical majority-Black vs majority-white origin, share of metro"))
      )
    )
  ),

  nav_panel("Exposure trends",
    layout_sidebar(
      sidebar = sidebar(width=300,
        selectInput("tr_poll", "Pollutant", POLL, selected="pm25"),
        radioButtons("tr_view", "Show", c("Group means"="level", "Relative disparity (ratio to white)"="ratio", "Absolute gap"="gap")),
        helpText("Population-weighted mean exposure of residents of majority-Black, majority-Hispanic and majority-white block groups, by year. Absolute gaps narrowed for most pollutants; relative gaps did not.")
      ),
      card(full_screen=TRUE, card_header(textOutput("tr_title")), plotlyOutput("tr_plot", height="480px"))
    )
  ),

  nav_panel("Cleanup vs. sorting",
    layout_sidebar(
      sidebar = sidebar(width=300,
        selectInput("dc_poll", "Pollutant", POLL, selected="pm25"),
        helpText("Shift-share decomposition of the 2011–2019 change in each group's exposure, and in the Black–white and Hispanic–white gaps, into in-place cleanup (pollution falling where the group already lived) and sorting (the group redistributing across neighbourhoods). The two terms sum exactly to the total.")
      ),
      layout_columns(col_widths=c(6,6),
        card(card_header("Change in each group's mean exposure"), plotlyOutput("dc_groups", height="420px")),
        card(card_header("Change in the exposure gap"), plotlyOutput("dc_gaps", height="420px"))
      ),
      card(card_header("Why a small sorting term does not mean mobility was unimportant"),
        p("The decomposition asks how much redistribution shifted a group's ", em("average"), " exposure. Because some moves led to cleaner air and others to dirtier air, those changes largely cancel when averaged over a group. It does not ask whether the households who moved fared equally well — that is the conditional avoidance gap on the next tab, which compares movers who left equally polluted neighbourhoods."))
    )
  ),

  nav_panel("Conditional avoidance gap",
    layout_sidebar(
      sidebar = sidebar(width=320,
        radioButtons("cg_fe", "Fixed effects", c("Metropolitan area (main)"="cbsa", "County"="county", "National, year only"="none")),
        helpText("Exposure reallocation regressed on origin group, a natural spline in origin pollution (which absorbs mean reversion), log income, poverty, education, and year effects; weighted by out-flow; SEs clustered on metro. A positive estimate means streams from that origin type reached dirtier air than streams from majority-white origins that started equally polluted.")
      ),
      card(full_screen=TRUE, card_header("Conditional PM2.5 avoidance gap by origin group (µg/m³, 95% CI)"), plotlyOutput("cg_plot", height="400px")),
      layout_columns(col_widths=c(4,4,4),
        card(card_header("Majority threshold"), plotlyOutput("cg_thresh", height="300px")),
        card(card_header("Selection into moving"), plotlyOutput("cg_select", height="300px")),
        card(card_header("Measurement-error weighting"), plotlyOutput("cg_meerr", height="300px"))
      )
    )
  ),

  nav_panel("Eight pollutants",
    layout_sidebar(
      sidebar = sidebar(width=320,
        radioButtons("mp_group", "Origin group", c("Black", "Hispanic")),
        checkboxInput("mp_sd", "Express in national SD units (comparable across pollutants)", TRUE),
        helpText("Three mechanisms estimated separately for each pollutant. Avoidance gap: movers from minority origins reach dirtier air. Sorting: destinations drawing more minority in-migration are more polluted. Unequal cleanup: neighbourhoods with higher minority shares cleaned up more slowly from the same starting level. The disparity concentrates in traffic and combustion pollutants, vanishes for wildfire smoke, and reverses for ozone.")
      ),
      card(full_screen=TRUE, card_header("Three mechanisms across pollutants (point estimate, 95% CI)"), plotlyOutput("mp_plot", height="520px"))
    )
  ),

  nav_panel("Affordability & segregation",
    layout_columns(col_widths=c(6,6),
      card(card_header("Does destination price explain the gap? (Gelbach decomposition)"), plotlyOutput("af_gelbach", height="320px"),
        p(class="text-muted small", "Share of the conditional gap attributed to destination median home value and gross rent, on the common sample where both are observed. Descriptive, not causal mediation: destination price is chosen after the move.")),
      card(card_header("Clean and affordable neighbourhoods within reach"), plotlyOutput("af_opp", height="320px"),
        p(class="text-muted small", "Share of a metro's block groups that are both below the metro median PM2.5 and have a median home value within 3× the origin's median income. A price-based heuristic for reach, not any household's true choice set."))
    ),
    layout_columns(col_widths=c(6,6),
      card(card_header("The gap grows with metropolitan segregation"), plotlyOutput("af_seg", height="320px"),
        p(class="text-muted small", "Black avoidance gap by quartile of metropolitan Black–white dissimilarity. Consistent with discrimination and steering, but other mechanisms that operate more strongly in segregated markets could produce the same pattern.")),
      card(card_header("Opportunity-set definition sensitivity"), DTOutput("af_sens"))
    )
  ),

  nav_panel("Metro explorer",
    if (!HAS_METRO) card(p("Metro-level map data has not been built yet. Run scripts/build_app_metro_data.R (after run_main.R); it writes metro_*.csv into shiny_app/data/."))
    else layout_sidebar(
      sidebar = sidebar(width=330,
        selectInput("me_cbsa", "Metropolitan area", setNames(D$mlist$cbsa, D$mlist$name)),
        radioButtons("me_layer", "Map layer",
          c("Exposure change of out-movers, by origin block group"="re",
            "Where movers from majority-Black origins went"="dest_black",
            "Where movers from majority-white origins went"="dest_white",
            "Where movers from majority-Hispanic origins went"="dest_hisp",
            "PM2.5 in 2019"="pm25")),
        sliderInput("me_min", "Hide block groups with fewer out-moves than", min=0, max=200, value=20, step=10),
        helpText("Within-metro moves only, pooled 2011–2019. Each dot is a block-group centroid. In the origin layer, blue means the households leaving that block group reached cleaner air on average and red dirtier; in the destination layers, dot size is the number of arrivals from that origin group and colour is the destination's PM2.5.")
      ),
      card(full_screen=TRUE, card_header(textOutput("me_title")), leafletOutput("me_map", height="560px")),
      layout_columns(col_widths=c(7,5),
        card(card_header("Where each origin group's movers land on the metro's pollution gradient"), plotlyOutput("me_sankey", height="420px"),
          p(class="text-muted small", "Destination block groups ranked into quintiles of 2019 PM2.5 within the metro (Q1 cleanest, Q5 dirtiest). Hover a ribbon for the number of moves and the mean change in PM2.5 those movers experienced.")),
        card(card_header("This metro in numbers"), DTOutput("me_table"))
      )
    )
  ),

  nav_panel("Urban form & distance",
    if (is.null(D$r1_band) && is.null(D$r1_dest)) card(p("The urban-form and move-distance analyses (added in revision) have not been run yet. Run scripts/reviewer1_urban_form_distance.R and copy its r1_*.csv outputs into the app's data folder."))
    else tagList(
      layout_columns(col_widths=c(6,6),
        card(card_header("Black gap by move distance"), plotlyOutput("r1_band", height="340px"),
          p(class="text-muted small", "Conditional gap re-estimated within distance bands. Short moves are the least likely to be driven by a change of workplace.")),
        card(card_header("Black gap by destination type"), plotlyOutput("r1_dest", height="340px"),
          p(class="text-muted small", "Urban core = principal city; inner suburb = same county or within 15 km of its centre; periphery = rest of the metro."))
      ),
      layout_columns(col_widths=c(7,5),
        card(card_header("Where moves from the urban core go, by origin group"), plotlyOutput("r1_mat", height="380px")),
        card(card_header("Gap within origin–destination type"), DTOutput("r1_odfe"), DTOutput("r1_dist"))
      )
    )
  ),

  nav_panel("Data & methods",
    layout_columns(col_widths=c(4,8),
      card(card_header("Sources"),
        tags$ul(
          tags$li("Migration flows: MIGRATE block-group origin–destination database (Agostini et al., 2025)"),
          tags$li("PM2.5 and components with uncertainty: WashU SatPM V6.NA.01 (Shen et al., 2026)"),
          tags$li("NO2: Anenberg et al. (2022), ~1 km"),
          tags$li("Ozone: Liu et al. (2026), warm-season MDA8, 0.01°"),
          tags$li("Wildfire smoke PM2.5: Childs et al. (2022), 10 km daily"),
          tags$li("Neighbourhood covariates: ACS 5-year, NHGIS; HMDA"),
          tags$li("Individual movers: IPUMS USA ACS 2012–2019")
        ),
        p("Analysis code: ", a("github.com/pdez90/pm25_residentialmobility", href="https://github.com/pdez90/pm25_residentialmobility", target="_blank")),
        p(class="text-muted small", "Contrasts are defined by the racial and ethnic composition of the origin neighbourhood, not the race of individual movers; MIGRATE is not disaggregated by race. Results describe 2011–2019.")
      ),
      card(card_header("Browse the output tables"),
        selectInput("tb_file", "Table", choices=sort(list.files("data", pattern="\\.csv$")), width="60%"),
        DTOutput("tb_dt"), downloadButton("tb_dl", "Download CSV"))
    )
  )
)

# ---------- server ----------
server <- function(input, output, session) {

  # Exposure trends
  output$tr_title <- renderText(sprintf("%s, residents of majority-Black, majority-Hispanic and majority-white block groups", plab(input$tr_poll)))
  output$tr_plot <- renderPlotly({
    req(D$trends); t <- D$trends[pollutant==input$tr_poll]; u <- UNIT[[input$tr_poll]]
    if (input$tr_view=="level") {
      m <- melt(t[, .(year, Black, Hispanic, White)], id.vars="year", variable.name="group", value.name="v")
      g <- ggplot(m, aes(year, v, colour=group)) + geom_line(linewidth=1.1) + geom_point(size=2) + scale_colour_manual(values=GRP_COL) +
        labs(x=NULL, y=sprintf("Population-weighted mean (%s)", u), colour=NULL) + theme_app()
    } else if (input$tr_view=="ratio") {
      m <- melt(t[, .(year, `Black / White`=ratio_BW, `Hispanic / White`=ratio_HW)], id.vars="year", variable.name="group", value.name="v")
      g <- ggplot(m, aes(year, v, colour=group)) + geom_hline(yintercept=1, linetype=2, colour="grey50") + geom_line(linewidth=1.1) + geom_point(size=2) +
        scale_colour_manual(values=c(`Black / White`=GRP_COL[["Black"]], `Hispanic / White`=GRP_COL[["Hispanic"]])) + labs(x=NULL, y="Ratio of group means", colour=NULL) + theme_app()
    } else {
      m <- melt(t[, .(year, `Black − White`=gap_BW, `Hispanic − White`=gap_HW)], id.vars="year", variable.name="group", value.name="v")
      g <- ggplot(m, aes(year, v, colour=group)) + geom_hline(yintercept=0, linetype=2, colour="grey50") + geom_line(linewidth=1.1) + geom_point(size=2) +
        scale_colour_manual(values=c(`Black − White`=GRP_COL[["Black"]], `Hispanic − White`=GRP_COL[["Hispanic"]])) + labs(x=NULL, y=sprintf("Gap (%s)", u), colour=NULL) + theme_app()
    }
    ggplotly(g) |> layout(legend=list(orientation="h", y=-0.15))
  })

  # Cleanup vs sorting
  output$dc_groups <- renderPlotly({
    req(D$decomp)
    # group-level change is available for PM2.5 in exposure_decomposition; for other pollutants reconstruct from trends
    if (input$dc_poll=="pm25" && !is.null(D$expdec)) {
      m <- melt(D$expdec[, .(group, cleanup, sorting)], id.vars="group", variable.name="term", value.name="v")
    } else {
      t <- D$trends[pollutant==input$dc_poll]; tot <- t[year==2019, .(Black, Hispanic, White)] - t[year==2011, .(Black, Hispanic, White)]
      m <- data.table(group=c("Black","Hispanic","White"), term="total change", v=as.numeric(tot))
    }
    g <- ggplot(m, aes(group, v, fill=term)) + geom_col(position="stack") + geom_hline(yintercept=0, colour="grey40") +
      scale_fill_manual(values=c(cleanup="#3B6FA0", sorting="#D98E2B", `total change`="#7A7A7A")) +
      labs(x=NULL, y=sprintf("Change 2011–2019 (%s)", UNIT[[input$dc_poll]]), fill=NULL) + theme_app()
    ggplotly(g)
  })
  output$dc_gaps <- renderPlotly({
    req(D$decomp); d <- D$decomp[pollutant==input$dc_poll]
    m <- melt(d[, .(gap, cleanup, sorting)], id.vars="gap", variable.name="term", value.name="v")
    g <- ggplot(m, aes(gap, v, fill=term)) + geom_col(position="stack") + geom_hline(yintercept=0, colour="grey40") +
      geom_point(data=d, aes(gap, total_change), inherit.aes=FALSE, shape=21, size=3, fill="white") +
      scale_fill_manual(values=c(cleanup="#3B6FA0", sorting="#D98E2B")) +
      labs(x=NULL, y=sprintf("Change in gap 2011–2019 (%s); dot = total", UNIT[[input$dc_poll]]), fill=NULL) + theme_app()
    ggplotly(g)
  })

  # Conditional gap
  output$cg_plot <- renderPlotly({
    req(D$cond); d <- copy(D$cond)[fe==input$cg_fe & dimension=="o_race_majority"]
    d[, group := sub("^grp", "", term)]; d[, `:=`(lo=Estimate-1.96*`Std. Error`, hi=Estimate+1.96*`Std. Error`)]
    d[, group := factor(group, levels=rev(c("Majority Black","Majority Hispanic","No majority")))]
    g <- ggplot(d, aes(Estimate, group)) + geom_vline(xintercept=0, linetype=2, colour="grey50") +
      geom_errorbar(aes(xmin=lo, xmax=hi), width=0.2) + geom_point(size=3, colour="#B2403A") +
      labs(x="Gap relative to majority-white origins (µg/m³)", y=NULL) + theme_app()
    ggplotly(g)
  })
  output$cg_thresh <- renderPlotly({
    req(D$thresh); d <- D$thresh[, .(spec, est=black_gap, se=black_se)]; ci(d)
    g <- ggplot(d, aes(est, spec)) + geom_vline(xintercept=0, linetype=2, colour="grey50") + geom_errorbar(aes(xmin=lo, xmax=hi), width=0.2) + geom_point(size=2.5, colour="#B2403A") +
      labs(x="Black gap (µg/m³)", y=NULL) + theme_app(); ggplotly(g)
  })
  output$cg_select <- renderPlotly({
    req(D$select); d <- D$select[group=="Black"]; ci(d)
    d[, spec := factor(spec, levels=c("all","within_cbsa","short_dist","plus_dist"), labels=c("All moves","Within metro","≤ 50 km","+ distance control"))]
    g <- ggplot(d, aes(est, spec)) + geom_vline(xintercept=0, linetype=2, colour="grey50") + geom_errorbar(aes(xmin=lo, xmax=hi), width=0.2) + geom_point(size=2.5, colour="#B2403A") +
      labs(x="Black gap (µg/m³)", y=NULL) + theme_app(); ggplotly(g)
  })
  output$cg_meerr <- renderPlotly({
    req(D$meerr); d <- D$meerr[, .(spec, est=Black, se=Black_se)]; ci(d)
    g <- ggplot(d, aes(est, spec)) + geom_vline(xintercept=0, linetype=2, colour="grey50") + geom_errorbar(aes(xmin=lo, xmax=hi), width=0.2) + geom_point(size=2.5, colour="#B2403A") +
      labs(x="Black gap (µg/m³)", y=NULL) + theme_app(); ggplotly(g)
  })

  # Eight pollutants
  output$mp_plot <- renderPlotly({
    req(D$cross); d <- copy(D$cross)[group==input$mp_group]
    if (input$mp_sd) { d[, `:=`(e=est_sd, s=se_sd)]; xl <- "Estimate (national SD units)" } else { d[, `:=`(e=est, s=se)]; xl <- "Estimate (native units)" }
    d[, `:=`(lo=e-1.96*s, hi=e+1.96*s)]
    d[, pollutant_lab := factor(pollutant_lab, levels=rev(unique(D$cross[order(match(pollutant, POLL))]$pollutant_lab)))]
    d[, sig := fifelse(p < 0.05, "p < 0.05", "n.s.")]
    g <- ggplot(d, aes(e, pollutant_lab, colour=sig)) + geom_vline(xintercept=0, linetype=2, colour="grey50") +
      geom_errorbar(aes(xmin=lo, xmax=hi), width=0.25) + geom_point(size=2.8) + facet_wrap(~mechanism, scales="free_x") +
      scale_colour_manual(values=c(`p < 0.05`="#B2403A", n.s.="grey60")) + labs(x=xl, y=NULL, colour=NULL) + theme_app()
    ggplotly(g) |> layout(legend=list(orientation="h", y=-0.12))
  })

  # Affordability & segregation
  output$af_gelbach <- renderPlotly({
    req(D$gelbach); m <- melt(D$gelbach[, .(group, `Home value`=pct_homeval, `Gross rent`=pct_rent)], id.vars="group", variable.name="cov", value.name="pct")
    g <- ggplot(m, aes(group, pct, fill=cov)) + geom_col() + scale_fill_manual(values=c(`Home value`="#3B6FA0", `Gross rent`="#D98E2B")) +
      labs(x=NULL, y="% of conditional gap explained", fill=NULL) + theme_app(); ggplotly(g)
  })
  output$af_opp <- renderPlotly({
    req(D$opp); d <- copy(D$opp); d[, o_race := factor(o_race, levels=c("Majority white (NH)","No majority","Majority Hispanic","Majority Black"))]
    g <- ggplot(d, aes(o_race, 100*opp, fill=o_race)) + geom_col(show.legend=FALSE) +
      scale_fill_manual(values=c(`Majority white (NH)`=GRP_COL[["White"]], `No majority`=GRP_COL[["No majority"]], `Majority Hispanic`=GRP_COL[["Hispanic"]], `Majority Black`=GRP_COL[["Black"]])) +
      labs(x=NULL, y="% of metro block groups clean AND affordable") + theme_app() + theme(legend.position="none"); ggplotly(g)
  })
  output$af_seg <- renderPlotly({
    req(D$segq); d <- D$segq[, .(seg_q, est=black_gap, se)]; ci(d)
    g <- ggplot(d, aes(seg_q, est)) + geom_hline(yintercept=0, linetype=2, colour="grey50") + geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2) + geom_point(size=3, colour="#B2403A") +
      labs(x="Metro Black–white dissimilarity quartile", y="Black gap (µg/m³)") + theme_app(); ggplotly(g)
  })
  output$af_sens <- renderDT({ req(D$oppsens); datatable(D$oppsens, options=list(dom="t", pageLength=10), rownames=FALSE) |> formatPercentage(c("white","black","hisp"), 1) })

  # Urban form & distance (revision)
  output$r1_band <- renderPlotly({
    req(D$r1_band); d <- D$r1_band[group=="Majority Black"]; d[, band := factor(band, levels=unique(band))]
    g <- ggplot(d, aes(band, est)) + geom_hline(yintercept=0, linetype=2, colour="grey50") + geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2) + geom_point(size=3, colour="#B2403A") +
      labs(x="Move distance", y="Black gap (µg/m³)") + theme_app(); ggplotly(g)
  })
  output$r1_dest <- renderPlotly({
    req(D$r1_dest); d <- D$r1_dest; d[, dest_type := factor(dest_type, levels=unique(dest_type))]
    g <- ggplot(d, aes(dest_type, est)) + geom_hline(yintercept=0, linetype=2, colour="grey50") + geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2) + geom_point(size=3, colour="#B2403A") +
      labs(x="Destination type", y="Black gap (µg/m³)") + theme_app() + theme(axis.text.x=element_text(angle=20, hjust=1)); ggplotly(g)
  })
  output$r1_mat <- renderPlotly({
    req(D$r1_mat); d <- D$r1_mat[o_uform=="Urban core" & o_race %in% c("Majority white (NH)","Majority Black","Majority Hispanic") & !d_uform %in% "Metro (unmatched)"]
    g <- ggplot(d, aes(d_uform, o_race, fill=100*share_of_origin_flow)) + geom_tile() +
      geom_text(aes(label=sprintf("%.0f%%\n%+.2f", 100*share_of_origin_flow, re_pm25)), size=3.2) +
      scale_fill_gradient(low="white", high="#3B6FA0", name="% of out-flow") + labs(x="Destination type", y=NULL, caption="Cell text: share of out-flow; mean PM2.5 reallocation (µg/m³)") +
      theme_app() + theme(axis.text.x=element_text(angle=20, hjust=1)); ggplotly(g)
  })
  output$r1_odfe <- renderDT({ req(D$r1_odfe); datatable(D$r1_odfe[, .(spec, group, est=round(est,4), lo=round(lo,4), hi=round(hi,4), p=signif(p,2))], options=list(dom="t"), rownames=FALSE, caption="Conditional gap with origin×destination-type fixed effects") })
  output$r1_dist <- renderDT({ req(D$r1_dist); datatable(D$r1_dist[, lapply(.SD, function(v) if (is.numeric(v)) round(v,2) else v)], options=list(dom="t"), rownames=FALSE, caption="Flow-weighted move distance (km) and shares, by origin group") })


  # Metro explorer
  me_bg <- reactive({ req(HAS_METRO, input$me_cbsa); D$mbg[cbsa==input$me_cbsa & out_moves >= input$me_min] })
  me_name <- reactive({ req(HAS_METRO, input$me_cbsa); D$mlist[cbsa==input$me_cbsa]$name })
  output$me_title <- renderText({ req(HAS_METRO); switch(input$me_layer,
    re=sprintf("%s — mean PM2.5 change experienced by households leaving each block group (µg/m³)", me_name()),
    pm25=sprintf("%s — block-group PM2.5, 2019 (µg/m³)", me_name()),
    sprintf("%s — destinations of within-metro movers from %s origins", me_name(),
            c(dest_black="majority-Black", dest_white="majority-white", dest_hisp="majority-Hispanic")[[input$me_layer]])) })
  output$me_map <- renderLeaflet({
    req(HAS_METRO); b <- me_bg(); req(nrow(b) > 0)
    m <- leaflet(options=leafletOptions(preferCanvas=TRUE)) |> addProviderTiles(providers$CartoDB.Positron)
    if (input$me_layer == "re") {
      lim <- quantile(abs(b$re_pm25), 0.95, na.rm=TRUE); b[, v := pmax(pmin(re_pm25, lim), -lim)]
      pal <- colorNumeric(c("#2166AC","#F7F7F7","#B2182B"), domain=c(-lim, lim))
      m <- m |> addCircleMarkers(data=b, lng=~lon, lat=~lat, radius=~pmin(2 + sqrt(out_moves)/6, 9), stroke=FALSE, fillOpacity=0.75, fillColor=~pal(v),
             popup=~sprintf("<b>%s</b><br>%s<br>Out-moves 2011–19: %s<br>Mean PM2.5 change of out-movers: %+.2f µg/m³<br>PM2.5 2019: %.1f µg/m³", geoid_bg, race_majority, format(round(out_moves)), re_pm25, pm25_2019)) |>
        addLegend("bottomright", pal=pal, values=c(-lim, lim), title="Δ PM2.5 of out-movers<br>(blue = cleaner)", opacity=0.9)
    } else if (input$me_layer == "pm25") {
      pal <- colorNumeric("viridis", domain=b$pm25_2019, reverse=TRUE)
      m <- m |> addCircleMarkers(data=b, lng=~lon, lat=~lat, radius=4, stroke=FALSE, fillOpacity=0.75, fillColor=~pal(pm25_2019),
             popup=~sprintf("<b>%s</b><br>%s<br>PM2.5 2019: %.1f µg/m³<br>NO2 2019: %.1f ppb", geoid_bg, race_majority, pm25_2019, no2_2019)) |>
        addLegend("bottomright", pal=pal, values=b$pm25_2019, title="PM2.5 2019 (µg/m³)", opacity=0.9)
    } else {
      grp <- c(dest_black="Majority Black", dest_white="Majority white (NH)", dest_hisp="Majority Hispanic")[[input$me_layer]]
      d <- D$mdest[cbsa==input$me_cbsa & o_race==grp]; req(nrow(d) > 0)
      pal <- colorNumeric("viridis", domain=d$dest_pm25, reverse=TRUE)
      m <- m |> addCircleMarkers(data=d, lng=~lon, lat=~lat, radius=~pmin(1.5 + sqrt(in_moves)/3, 14), stroke=FALSE, fillOpacity=0.7, fillColor=~pal(dest_pm25),
             popup=~sprintf("<b>Destination %s</b><br>Arrivals from %s origins: %s<br>Destination PM2.5 2019: %.1f µg/m³<br>Mean PM2.5 of their origins: %.1f µg/m³", dest_bg, grp, format(round(in_moves)), dest_pm25, origin_pm25)) |>
        addLegend("bottomright", pal=pal, values=d$dest_pm25, title="Destination PM2.5<br>2019 (µg/m³)", opacity=0.9)
    }
    m
  })
  output$me_sankey <- renderPlotly({
    req(HAS_METRO, input$me_cbsa); s <- D$msank[cbsa==input$me_cbsa & o_race %in% c("Majority white (NH)","Majority Black","Majority Hispanic","No majority")]
    req(nrow(s) > 0)
    og <- c("Majority white (NH)","Majority Black","Majority Hispanic","No majority"); dq <- paste0("Q", 1:5)
    nodes <- c(og, paste0(dq, c(" (cleanest)","","","", " (dirtiest)")))
    s[, `:=`(src = match(o_race, og) - 1L, tgt = length(og) + match(dest_q, dq) - 1L)]
    s[, share := moves / sum(moves), by=o_race]
    lc <- c(`Majority white (NH)`="rgba(59,111,160,0.45)", `Majority Black`="rgba(178,64,58,0.45)", `Majority Hispanic`="rgba(217,142,43,0.45)", `No majority`="rgba(122,122,122,0.45)")
    plot_ly(type="sankey", orientation="h", arrangement="snap",
      node=list(label=nodes, pad=18, thickness=18, color=c(GRP_COL[["White"]], GRP_COL[["Black"]], GRP_COL[["Hispanic"]], GRP_COL[["No majority"]], rev(viridisLite::viridis(5))),
                line=list(width=0.5, color="grey")),
      link=list(source=s$src, target=s$tgt, value=s$moves, color=lc[s$o_race],
                customdata=sprintf("%s → %s<br>%s moves (%.0f%% of this group's within-metro moves)<br>mean Δ PM2.5: %+.2f µg/m³", s$o_race, s$dest_q, format(round(s$moves), big.mark=","), 100*s$share, s$re_pm25),
                hovertemplate="%{customdata}<extra></extra>")) |>
      layout(font=list(size=12), margin=list(l=10, r=10, t=10, b=10))
  })
  output$me_table <- renderDT({
    req(HAS_METRO, input$me_cbsa); b <- D$mbg[cbsa==input$me_cbsa & out_moves > 0]
    t <- b[, .(`Block groups`=.N, `Out-moves`=format(round(sum(out_moves)), big.mark=","), `Mean PM2.5 2019`=round(weighted.mean(pm25_2019, out_moves, na.rm=TRUE), 2),
              `Mean Δ PM2.5 of out-movers`=round(weighted.mean(re_pm25, out_moves, na.rm=TRUE), 3)), by=.(`Origin majority`=race_majority)][order(`Origin majority`)]
    datatable(t, options=list(dom="t"), rownames=FALSE)
  })

  # Data browser
  tb <- reactive({ req(input$tb_file); fread(file.path("data", input$tb_file)) })
  output$tb_dt <- renderDT(datatable(tb(), options=list(pageLength=15, scrollX=TRUE), rownames=FALSE) |> formatRound(names(tb())[sapply(tb(), is.numeric)], 4))
  output$tb_dl <- downloadHandler(filename=function() input$tb_file, content=function(f) fwrite(tb(), f))
}

shinyApp(ui, server)
