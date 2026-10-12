# Small extras for the revision (Reviewer 3 #8 and Reviewer 2 #6), ~3 min.
#   (1) number of CBSA clusters in the main model (and county clusters, and observations)
#   (2) unweighted estimates of the main model for Table S4 / Table 1
#   (3) share of block groups and population below 1 km2 (Table S25 / section 3.4 pointer)
# Reads the cached origin-year panel written by idea2_pm25_avoidance_v3.R.
# Output: pm25_avoidance_lag3_olag1/revision_extras.csv and revision_unweighted.csv
if (!exists(".migrate_paths_loaded")) source(if (file.exists("scripts/_paths.R")) "scripts/_paths.R" else "_paths.R")
suppressPackageStartupMessages({ library(data.table); library(fixest); library(splines) })
base_dir <- migrate_data_dir(); out_dir <- file.path(base_dir, "pm25_avoidance_lag3_olag1")
oy <- fread(file.path(out_dir, "origin_year_pm25_reallocation.csv"), colClasses=c(origin_bg="character", o_state_county="character"))
if (!"o_county" %in% names(oy)) oy[, o_county := stringr::str_pad(o_state_county, 5, pad="0")]
if (!"o_cbsa" %in% names(oy)) {
  cxw <- fread(file.path(base_dir, "pollution_bg", "county_cbsa_xwalk.csv"), colClasses=c(county="character", cbsa="character"))
  oy <- merge(oy, unique(cxw, by="county")[, .(o_county=county, o_cbsa=cbsa)], by="o_county", all.x=TRUE); oy[is.na(o_cbsa), o_cbsa := paste0("nonmetro_", o_county)]
}
d <- oy[!is.na(o_race_majority) & out_flow > 0 & !is.na(o_med_hh_income) & o_med_hh_income > 0 & !is.na(orig_pm25) & !is.na(re_pm25)]
d[, grp := relevel(factor(o_race_majority), ref="Majority white (NH)")]
f <- re_pm25 ~ grp + ns(orig_pm25,4) + log(o_med_hh_income) + o_poverty_pct + o_bach_pct | o_cbsa + year
m_w <- feols(f, d, weights=~out_flow, cluster=~o_cbsa)
m_u <- feols(f, d, cluster=~o_cbsa)
ct <- function(m, w) { x <- as.data.table(coeftable(m), keep.rownames="term")[grepl("^grp", term)]; setnames(x, c("term","est","se","t","p")); x[, `:=`(weighting=w, lo=est-1.96*se, hi=est+1.96*se)][] }
unw <- rbind(ct(m_w, "Flow-weighted (main)"), ct(m_u, "Unweighted"))
fwrite(unw, file.path(out_dir, "revision_unweighted.csv"))
cat("\n=== Main model, flow-weighted vs unweighted (CBSA + year FE, CBSA-clustered) ===\n"); print(unw[, .(weighting, term, est=round(est,4), lo=round(lo,4), hi=round(hi,4), p=signif(p,2))])
ex <- data.table(item=c("n_obs_origin_years","n_origin_bg","n_cbsa_clusters","n_cbsa_clusters_metro","n_cbsa_clusters_nonmetro","n_county_clusters","n_years"),
                 value=c(nobs(m_w), uniqueN(d$origin_bg), uniqueN(d$o_cbsa), sum(!grepl("^nonmetro", unique(d$o_cbsa))), sum(grepl("^nonmetro", unique(d$o_cbsa))), uniqueN(d$o_county), uniqueN(d$year)))
af <- file.path(out_dir, "blockgroup_area_summary.csv")
if (file.exists(af)) { a <- fread(af); ex <- rbind(ex, data.table(item=c("pct_bg_below_1km2","pct_pop_below_1km2"), value=c(a$pct_bg[1], a$pct_pop[1]))) }
fwrite(ex, file.path(out_dir, "revision_extras.csv")); cat("\n"); print(ex)
cat("\nWrote", file.path(out_dir, "revision_extras.csv"), "\n")
