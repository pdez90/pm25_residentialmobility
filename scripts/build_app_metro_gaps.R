# Metro-specific conditional avoidance gaps for the Shiny metro explorer (~2-3 min).
# For each of the app's metros, re-estimates the main model WITHIN that metro:
#   re_pm25 ~ grp + ns(orig_pm25,4) + log(income) + poverty + bach | year, weights = out_flow,
# SEs clustered on origin county. Same specification as the paper, restricted to one CBSA,
# so the estimate is the metro's own adjusted gap (descriptive; small metros are noisy).
# Reads the cached origin-year panel written by idea2_pm25_avoidance.R.
# Output: shiny_app/data/metro_gap.csv
if (!exists(".migrate_paths_loaded")) source(if (file.exists("scripts/_paths.R")) "scripts/_paths.R" else "_paths.R")
suppressPackageStartupMessages({ library(data.table); library(stringr); library(fixest); library(splines) })
base_dir <- migrate_data_dir(); out_dir <- file.path(base_dir,"pm25_avoidance_lag3_olag1"); app_dir <- file.path(repo_root,"shiny_app","data")
top <- fread(file.path(app_dir,"metro_list.csv"), colClasses=c(cbsa="character"))
oy <- fread(file.path(out_dir,"origin_year_pm25_reallocation.csv"), colClasses=c(origin_bg="character", o_state_county="character"))
oy[, county := str_pad(o_state_county, 5, pad="0")]
cxw <- fread(file.path(base_dir,"pollution_bg","county_cbsa_xwalk.csv"), colClasses=c(county="character",cbsa="character"))
oy <- merge(oy, unique(cxw, by="county"), by="county")
oy <- oy[cbsa %in% top$cbsa & !is.na(o_race_majority) & out_flow>0 & !is.na(o_med_hh_income) & o_med_hh_income>0 & !is.na(orig_pm25) & !is.na(re_pm25)]
oy[, grp := relevel(factor(o_race_majority), ref="Majority white (NH)")]
f <- re_pm25 ~ grp + ns(orig_pm25,4) + log(o_med_hh_income) + o_poverty_pct + o_bach_pct | year
res <- rbindlist(lapply(top$cbsa, function(cb) {
  d <- oy[cbsa==cb]; if (nrow(d) < 200) return(NULL)
  fit <- tryCatch(feols(f, d, weights=~out_flow, cluster=~county), error=function(e) NULL); if (is.null(fit)) return(NULL)
  ct <- as.data.table(coeftable(fit), keep.rownames="term")[term %in% c("grpMajority Black","grpMajority Hispanic","grpNo majority")]
  setnames(ct, c("term","est","se","t","p"))
  ct[, `:=`(cbsa=cb, group=sub("^grp","",term), n_origin_years=nrow(d),
            n_group=sapply(term, function(t) d[o_race_majority==sub("^grp","",t), .N]),
            lo=est-1.96*se, hi=est+1.96*se)]
  ct[, .(cbsa, group, est=round(est,4), se=round(se,4), lo=round(lo,4), hi=round(hi,4), p=signif(p,3), n_origin_years, n_group)]
}))
fwrite(res, file.path(app_dir,"metro_gap.csv"))
cat("wrote", nrow(res), "rows for", uniqueN(res$cbsa), "metros to", file.path(app_dir,"metro_gap.csv"), "\n")
print(merge(res[group=="Majority Black"], top[, .(cbsa, name)], by="cbsa")[order(-est), .(name, est, lo, hi, n_group)])
