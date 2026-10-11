#!/usr/bin/env bash
# ============================================================
# Download the raw exposure surfaces that CAN be fetched by URL
# (verified 11 Oct 2026). Run natively from the repository root:
#
#   bash setup/download_rasters.sh            # everything below
#   bash setup/download_rasters.sh no2        # or one of: no2 ozone smoke
#
# Resumable: curl -C - skips/continues files already present.
# Total: NO2 ~14 x ~0.6-1 GB, ozone 14 x 64 MB, smoke ~0.6 GB.
#
# NOT covered here (no scriptable URL):
#   * ACAG SatPM V6.NA.01 PM2.5 + components + uncertainty: the Box links on
#     satpm.org (ACAG-V6NA01-CNNPM25, -Uncertainty) were dead on 11 Oct 2026.
#     Email support@satpm.org for a working link; see REBUILD_CHECKLIST.md.
# ============================================================
set -euo pipefail
ROOT="${MIGRATE_DATA_DIR:-$(pwd)/data}"
WHAT="${1:-all}"
dl() {  # dl <url> <dest>
  local url="$1" dest="$2"
  mkdir -p "$(dirname "$dest")"
  if [ -s "$dest" ]; then echo "  have  $(basename "$dest")"; return; fi
  echo "  fetch $(basename "$dest")"
  curl -fL --retry 5 --retry-delay 10 -C - -o "$dest" "$url"
}

# ---- NO2: Anenberg et al. 2022 global surface NO2, ~1 km, annual NetCDF ----
# Shared Box folder: https://gwu.box.com/s/4wey4183hb0nlp3dg4wl02chca3nqnl9
# Box direct-download form: index.php?rm=box_download_shared_file&shared_name=<share>&file_id=f_<id>
if [ "$WHAT" = all ] || [ "$WHAT" = no2 ]; then
  echo "== NO2 -> $ROOT/raw/no2"
  SHARE=4wey4183hb0nlp3dg4wl02chca3nqnl9
  for pair in 2006:902361035267 2007:902358900739 2008:902361187350 2009:902360534114 \
              2010:902364695903 2011:902361777446 2012:902361656382 2013:902354878555 \
              2014:902360064205 2015:902366684527 2016:902361029352 2017:902360383547 \
              2018:902360689807 2019:902359970546; do
    y=${pair%%:*}; id=${pair##*:}
    dl "https://gwu.app.box.com/index.php?rm=box_download_shared_file&shared_name=$SHARE&file_id=f_$id" \
       "$ROOT/raw/no2/${y}_final_1km.nc"
  done
fi

# ---- Ozone: Liu, Chu, Deziel & Chen 2026, Yale Dataverse doi:10.60600/YU/M1WT9R ----
# annual_metrics files carry AMDA8 and WSMDA8 (warm-season May-Oct mean MDA8, ppb)
if [ "$WHAT" = all ] || [ "$WHAT" = ozone ]; then
  echo "== Ozone -> $ROOT/raw/ozone/annual"
  for pair in 2006:43268 2007:43253 2008:43271 2009:43265 2010:43274 2011:43293 2012:43277 \
              2013:43254 2014:43269 2015:43283 2016:43266 2017:43288 2018:43270 2019:43278; do
    y=${pair%%:*}; id=${pair##*:}
    dl "https://dataverse.yale.edu/api/access/datafile/$id" \
       "$ROOT/raw/ozone/annual/CONUS_O3_MDA8_p01_annual_metrics_${y}.nc"
  done
  dl "https://dataverse.yale.edu/api/access/datafile/43689" "$ROOT/raw/ozone/README.txt"
fi

# ---- Wildfire smoke: Childs et al. 2022, ECHO Lab v2.0 (2006-2023), Dropbox ----
if [ "$WHAT" = all ] || [ "$WHAT" = smoke ]; then
  echo "== Smoke -> $ROOT/raw/smoke"
  RK=nutebc9pn2vsupr0p9ks4k73u
  dl "https://www.dropbox.com/scl/fo/91k0aq80vp57qixkm508q/AHRwfIPVxa4s6ko_YrjDuf8/10km_grid/smokePM2pt5_predictions_daily_10km_20060101-20231231.rds?rlkey=$RK&dl=1" \
     "$ROOT/raw/smoke/smokePM2pt5_predictions_daily_10km_20060101-20231231.rds"
  dl "https://www.dropbox.com/scl/fo/91k0aq80vp57qixkm508q/AOzBwciodiQAl0XgIoBo-LQ/10km_grid/README.txt?rlkey=$RK&dl=1" \
     "$ROOT/raw/smoke/README_smoke.txt"
  if [ ! -f "$ROOT/raw/smoke/10km_grid_wgs84/10km_grid_wgs84.shp" ]; then
    dl "https://www.dropbox.com/scl/fo/91k0aq80vp57qixkm508q/AOr-fAml-EKKgLaMC1UM9Rw/10km_grid/10km_grid_wgs84?rlkey=$RK&dl=1" \
       "$ROOT/raw/smoke/10km_grid_wgs84.zip"
    unzip -o -q "$ROOT/raw/smoke/10km_grid_wgs84.zip" -d "$ROOT/raw/smoke/10km_grid_wgs84"
    # Dropbox zips may nest a folder; flatten if so
    f=$(find "$ROOT/raw/smoke/10km_grid_wgs84" -name '10km_grid_wgs84.shp' | head -1)
    [ -n "$f" ] && [ "$(dirname "$f")" != "$ROOT/raw/smoke/10km_grid_wgs84" ] && mv "$(dirname "$f")"/* "$ROOT/raw/smoke/10km_grid_wgs84/"
  fi
fi

echo
echo "== sanity check =="
for d in no2 ozone/annual smoke; do printf "%-14s %s files\n" "$d" "$(find "$ROOT/raw/$d" -maxdepth 1 -type f | wc -l | tr -d ' ')"; done
command -v ncdump >/dev/null && ncdump -h "$ROOT/raw/ozone/annual/CONUS_O3_MDA8_p01_annual_metrics_2019.nc" | grep -E "WSMDA8|AMDA8" | head -2 || true
