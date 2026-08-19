# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# SPECIFICATION LADDER for the conditional avoidance gap
# ------------------------------------------------------------
# Shows how the Black and Hispanic avoidance-gap coefficient evolves as
# controls are added, to demonstrate the result is not an artefact of any
# single control and to make the mean-reversion logic visible:
#   (1) Raw (no controls)                    -> gap is null/negative
#   (2) + origin-PM2.5 spline                -> gap emerges (mean reversion removed)
#   (3) + income
#   (4) + poverty
#   (5) + education
#   (6) Main (+ CBSA and year fixed effects) -> stable
# Reads origin_year_affordability.csv (stage 7). Fast.
# Output: pm25_avoidance_lag3_olag1/ spec_ladder.csv + figures/figS_spec_ladder.png
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(fixest); library(splines); library(ggplot2); library(stringr) })
base_dir <- migrate_data_dir()
out_dir  <- file.path(base_dir,"pm25_avoidance_lag3_olag1"); fig_dir <- file.path(out_dir,"figures")
oyf <- file.path(out_dir,"origin_year_affordability.csv")
stopifnot("run affordability_robustness.R first" = file.exists(oyf))
oy <- fread(oyf, colClasses=c(origin_bg="character"))
oy[, grp := relevel(factor(o_race), ref="Majority white (NH)")]
oy <- oy[out_flow>0]

specs <- list(
  "Raw"                    = c("re_pm25 ~ grp", ""),
  "+ origin pollution"     = c("re_pm25 ~ grp + ns(orig_pm25,4)", ""),
  "+ income"               = c("re_pm25 ~ grp + ns(orig_pm25,4) + log(o_inc0)", ""),
  "+ poverty"              = c("re_pm25 ~ grp + ns(orig_pm25,4) + log(o_inc0) + o_pov", ""),
  "+ education"            = c("re_pm25 ~ grp + ns(orig_pm25,4) + log(o_inc0) + o_pov + o_bach", ""),
  "Main (+ CBSA, year FE)" = c("re_pm25 ~ grp + ns(orig_pm25,4) + log(o_inc0) + o_pov + o_bach", "| o_cbsa + year"))
res <- rbindlist(lapply(names(specs), function(nm){
  s <- specs[[nm]]
  fit <- feols(as.formula(paste(s[1], s[2])), oy, weights=~out_flow, cluster=~o_cbsa)
  ct <- as.data.table(coeftable(fit), keep.rownames="term")[grepl("grpMajority (Black|Hispanic)", term)]
  ct[, .(spec=nm, race=sub("grpMajority ","",term), est=Estimate, se=`Std. Error`)]
}))
res[, spec := factor(spec, levels=names(specs))]
fwrite(res, file.path(out_dir,"spec_ladder.csv"))
cat("\n=== SPECIFICATION LADDER (avoidance gap as controls are added) ===\n")
print(dcast(res, spec~race, value.var="est")[order(spec)])

p <- ggplot(res, aes(spec, est, color=race, group=race)) +
  geom_hline(yintercept=0, linewidth=.3, color="grey60") +
  geom_line(linewidth=.8) +
  geom_pointrange(aes(ymin=est-1.96*se, ymax=est+1.96*se), size=.4, position=position_dodge(.15)) +
  scale_color_manual(values=c(Black="#1b7837", Hispanic="#762a83"), name=NULL) +
  labs(title="The avoidance gap emerges once mean reversion is removed and is then stable",
       subtitle="Conditional gap (ug/m3) vs comparable white-origin movers, adding controls left to right",
       x=NULL, y="Avoidance gap (ug/m3)") +
  theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"), legend.position="bottom",
        axis.text.x=element_text(angle=20, hjust=1))
ggsave(file.path(fig_dir,"figS_spec_ladder.png"), p, width=8.5, height=4.6, dpi=200)
cat("\nDONE. Specification ladder in ", out_dir, "\n", sep="")
