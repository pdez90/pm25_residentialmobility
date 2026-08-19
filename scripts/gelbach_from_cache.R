# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# Quick Gelbach decomposition from the cached origin-year panel.
# Reproduces a0 (base) / a2 (full) / auxiliary regressions on the SAME common
# sample and computes the order-invariant per-covariate attribution
#   b_base - b_full = sum_k delta_k * gamma_k   (Gelbach 2016).
# Runs in seconds -- no MIGRATE re-aggregation. Self-checks the identity.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(fixest); library(splines) })
out_dir <- file.path(migrate_data_dir(), "pm25_avoidance_lag3_olag1")
oy <- fread(file.path(out_dir,"origin_year_affordability.csv"))
oy[, grp := relevel(factor(o_race), ref="Majority white (NH)")]
ctrl <- "+ ns(orig_pm25,4) + log(o_inc0) + o_pov + o_bach"
bcoef <- function(fit,t) tryCatch(coeftable(fit)[t,"Estimate"], error=function(e) NA_real_)

oy_aff <- oy[dest_homeval>0 & dest_rent>0]                                        # common sample
a0 <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," | o_cbsa + year")), oy_aff, weights=~out_flow, cluster=~o_cbsa)
a2 <- feols(as.formula(paste0("re_pm25 ~ grp",ctrl," + log(dest_homeval) + log(dest_rent) | o_cbsa + year")), oy_aff, weights=~out_flow, cluster=~o_cbsa)

X3   <- c("log(dest_homeval)","log(dest_rent)")
tmap <- c(`grpMajority Black`="Black", `grpMajority Hispanic`="Hispanic", `grpNo majority`="No majority")
gam  <- sapply(X3, function(v) coeftable(a2)[v,"Estimate"])
aux  <- lapply(X3, function(v) feols(as.formula(paste0(v," ~ grp",ctrl," | o_cbsa + year")), oy_aff, weights=~out_flow, cluster=~o_cbsa))
names(aux) <- X3

gelb <- rbindlist(lapply(names(tmap), function(t){
  b0 <- bcoef(a0,t); b2 <- bcoef(a2,t)
  contr <- sapply(X3, function(v) bcoef(aux[[v]], t)) * gam
  data.table(group=tmap[[t]], gap_base=b0, gap_full=b2, total_explained=b0-b2,
             contrib_homeval=contr[["log(dest_homeval)"]], contrib_rent=contr[["log(dest_rent)"]],
             pct_homeval=round(100*contr[["log(dest_homeval)"]]/b0,2),
             pct_rent   =round(100*contr[["log(dest_rent)"]]/b0,2),
             pct_total  =round(100*sum(contr)/b0,2),
             identity_gap=(b0-b2)-sum(contr)) }), fill=TRUE)
fwrite(gelb, file.path(out_dir,"gelbach_decomposition.csv"))
print(gelb)
maxgap <- max(abs(gelb$identity_gap), na.rm=TRUE)
cat(sprintf("\nadditivity check  max|(b_base-b_full) - sum_k delta_k*gamma_k| = %.3e  (should be ~0)\n", maxgap))
if (!is.finite(maxgap) || maxgap > 1e-5) stop("Gelbach identity violated -- see message in affordability_robustness.R.")
cat("common sample:", nrow(oy_aff), "of", nrow(oy), "origin-years\n")
