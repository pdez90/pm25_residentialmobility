#!/usr/bin/env bash
# ============================================================
# Fetch / stage the raw inputs for the paper pipeline.
# Run natively on the Mac from the repository root:
#
#   bash setup/get_inputs.sh
#
# Everything lands under $MIGRATE_DATA_DIR (default: ./data).
# Box-hosted products (ACAG) and Dropbox-hosted products (smoke) do not
# allow scripted downloads; this script opens them in the browser and
# tells you where to drop the files. Zenodo / NHGIS / Census parts are
# automated where possible.
# ============================================================
set -euo pipefail
ROOT="${MIGRATE_DATA_DIR:-$(pwd)/data}"
mkdir -p "$ROOT"/raw/{acag_pm25/{PM25,SO4,NO3,BC,OM},acag_pm25_uncertainty,no2,ozone/annual,smoke/10km_grid_wgs84,epa_greenbook} \
         "$ROOT"/{acs_gentrification,nhgis,ipums,pollution_bg}

echo "Data root: $ROOT"

# ---- 1. MIGRATE flows: unzip whatever is still zipped ---------------
MIG_ZIPS="${MIGRATE_ZIP_DIR:-$HOME/Downloads/Migration/MIGRATE}"
for y in 2011 2012 2013 2014 2015 2016 2017 2018 2019; do
  if [ ! -s "$ROOT/$y.csv" ] && [ -f "$MIG_ZIPS/$y.csv.zip" ]; then
    echo "unzipping MIGRATE $y ..."
    unzip -o -q -j "$MIG_ZIPS/$y.csv.zip" "$y.csv" -d "$ROOT"
  fi
done
ls -lh "$ROOT"/20??.csv 2>/dev/null || echo "!! no MIGRATE csvs yet (set MIGRATE_ZIP_DIR)"

# ---- 2. ACAG V6.NA.01 PM2.5, components, uncertainty (Box, manual) --
cat <<EOF

== ACAG SatPM V6.NA.01 (manual, Box) ==
  Total + components (annual NetCDF):
    https://wustl.box.com/v/ACAG-V6NA01-CNNPM25
  Uncertainty:
    https://wustl.box.com/v/ACAG-V6NA01-CNNPM25-Uncertainty
  Download the ANNUAL files for 2006-2019 and place them as:
    $ROOT/raw/acag_pm25/PM25/   V6NA01.CNNPM25.PM25.YYYY01-YYYY12.nc   (14 files)
    $ROOT/raw/acag_pm25/SO4/    ...SO4.YYYY01-YYYY12.nc
    $ROOT/raw/acag_pm25/NO3/    ...NO3.YYYY01-YYYY12.nc
    $ROOT/raw/acag_pm25/BC/     ...BC.YYYY01-YYYY12.nc
    $ROOT/raw/acag_pm25/OM/     ...OM.YYYY01-YYYY12.nc
    $ROOT/raw/acag_pm25_uncertainty/  V6NA01.CNNPM25.AbsoluteUncertainty.YYYY01-YYYY12.nc
  (the scripts match on '.<VAR>.YYYY01-YYYY12.nc' so the exact prefix
   does not matter, but keep only ONE file per year per folder.)
EOF
open "https://wustl.box.com/v/ACAG-V6NA01-CNNPM25" 2>/dev/null || true

# ---- 3. Surface NO2 (Anenberg et al. 2022; GWU Box, manual) ---------
cat <<EOF

== Surface NO2 1990-2020 ~1 km (manual, Box) ==
  figshare record: https://figshare.com/articles/dataset/Global_surface_NO2_concentrations_1990-2020/12968114
  files are hosted at: https://gwu.box.com/s/4wey4183hb0nlp3dg4wl02chca3nqnl9
  Download years 2006-2019 into $ROOT/raw/no2/ (one file per year; the
  filename must contain the 4-digit year).
  NOTE: extract_pollution_multi_bg.R expects GeoTIFFs (NO2_PATTERN "\\.tif$").
  If the current release is NetCDF, either convert with
    gdal_translate -of GTiff in.nc out_YYYY.tif
  or change NO2_PATTERN to "\\.nc$" in scripts/extract_pollution_multi_bg.R.
EOF

# ---- 4. Warm-season MDA8 ozone (Liu et al. 2026) --------------------
cat <<EOF

== Ozone WSMDA8 (manual) ==
  Liu, Chu, Deziel, Chen 2026, ES&T, doi:10.1021/acs.est.5c16412
  Follow the paper's Data Availability statement (Yale / Kai Chen lab);
  annual NetCDFs with variable 'WSMDA8', one per year, go in
    $ROOT/raw/ozone/annual/
EOF

# ---- 5. Wildfire smoke PM2.5 (Childs et al. 2022, v2.0 2006-2023) ---
cat <<EOF

== Smoke PM2.5 daily 10 km (manual, Dropbox) ==
  v2.0 beta (2006-2023): https://www.dropbox.com/scl/fo/91k0aq80vp57qixkm508q/AKQSIJ5C1kDMQLz8oh02UAA?rlkey=nutebc9pn2vsupr0p9ks4k73u&dl=0
  (v1 2006-2020 fallback: https://www.dropbox.com/sh/9mbcxy65crd5cex/AAAk2_1QE6A5rAS7M5We-iapa?dl=0
   or Harvard Dataverse doi:10.7910/DVN/DJVMTV)
  Place:
    $ROOT/raw/smoke/smokePM2pt5_predictions_daily_10km_20060101-20231231.rds
    $ROOT/raw/smoke/10km_grid_wgs84/10km_grid_wgs84.{shp,shx,dbf,prj}
  (if you end up with the v1 file, set
     export SMOKE_PM25_RDS=\$ROOT/raw/smoke/<actual name>.rds )
EOF

# ---- 6. EPA Green Book (run_full.R only) ----------------------------
if [ ! -f "$ROOT/raw/epa_greenbook/phistory.xls" ]; then
  curl -fsSL -o "$ROOT/raw/epa_greenbook/phistory.xls" \
    "https://www3.epa.gov/airquality/greenbook/downld/phistory.xls" \
    && echo "Green Book downloaded" || echo "!! Green Book: download phistory.xls manually from https://www.epa.gov/green-book"
fi

# ---- 7. API keys for the automated parts -----------------------------
cat <<EOF

== API keys (put in ~/.Renviron, then restart R) ==
  CENSUS_API_KEY=...   https://api.census.gov/data/key_signup.html   (tidycensus, race 2013-19, upstream panel)
  IPUMS_API_KEY=...    https://account.ipums.org/api_keys            (NHGIS 2011-12 race; IPUMS USA microdata)

Then, from R at the repo root:
  source("upstream/nhgis_download.R")                      # -> data/nhgis/*.zip
  source("upstream/blockgroup_gentrification_earliest_acs.R")
  source("upstream/hmda_capital_reinvestment.R")           # needs HMDA LAR 2012-2019
  -> data/acs_gentrification/blockgroup_gentrification_2009_2019_lag3_with_HMDA.csv

HMDA LAR:
  2012-2017 national LAR: https://www.consumerfinance.gov/data-research/hmda/historic-data/
  2018-2019: https://ffiec.cfpb.gov/data-browser/  (API; see hmda_capital_reinvestment.R)
EOF
