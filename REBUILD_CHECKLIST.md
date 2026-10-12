# Rebuild checklist (Oct 2026)

Original data and code were lost; code was recovered from GitHub
(github.com/pdez90/pm25_residentialmobility, commit 6c0c7ee). Data must be
re-downloaded. Everything below is read relative to `MIGRATE_DATA_DIR`,
which defaults to `./data` inside this folder, so the simplest layout is to
keep using that default.

Run from Terminal / RStudio on the Mac (R >= 4.5 recommended; package
versions pinned in `renv.lock`):

```bash
cd ~/data/pm25_residentialmobility
Rscript setup/install_pinned.R        # or R -e 'renv::restore()'
Rscript run_main.R                    # paper pipeline
# MIGRATE_RUN_VALIDATION=TRUE MIGRATE_RUN_QUASIEXPERIMENT=TRUE Rscript run_full.R
```

Expect the extraction stages (0a-0d) to be the slow part: national 1-km
rasters x 14 years x ~217k block groups. Plan for a machine with >= 16 GB RAM
and ~100 GB free disk; terra scratch goes under `data/pollution_bg/`.

## 0. Status (11 Oct 2026)

- MIGRATE 2011-2019: unzipped and verified in `data/`.
- EPA Green Book `phistory.xls`: downloaded.
- NO2, ozone, smoke: scriptable — run `bash setup/download_rasters.sh` (verified URLs, resumable).
  Ozone now comes from Yale Dataverse doi:10.60600/YU/M1WT9R (`annual/CONUS_O3_MDA8_p01_annual_metrics_YYYY.nc`, var `WSMDA8`).
  NO2 arrives as `YYYY_final_1km.nc`; `extract_pollution_multi_bg.R` was patched to accept `.nc`.
- ACAG V6.NA.01 concentrations: DONE. Working share (11 Oct 2026):
  https://wustl.app.box.com/s/qwvj00cvghiebo7wdf1voip42mcdvk89 (the `/v/ACAG-V6NA01-CNNPM25` vanity link is dead).
  `Annual.zip` from that share -> `data/raw/acag_pm25/<SP>/V6NA01.CNNPM25.<SP>.YYYY01-YYYY12.nc` for ALL EIGHT species
  (PM25, SO4, NO3, BC, OM, NH4, DUST, SS - `extract_pm25_components_bg.R` hard-fails if any folder is missing), 2006-2019.
  Source zips kept in `~/Downloads/ACAG_V6_NA/` (33 GB; only Annual.zip is needed - the monthly zips can be deleted).
- ACAG V6.NA.01 uncertainty: DONE. Share: https://wustl.app.box.com/s/njjw6nbgx8zeac9z4hv2ijxyu9ib8ckg
  (`Annual.zip` -> `data/raw/acag_pm25_uncertainty/PM25/V6NA01.CNNPM25.AbsoluteUncertainty.YYYY01-YYYY12.nc`, 2006-2019).
- Still manual: Census + IPUMS API keys, `upstream/` panel build, HMDA LAR, IPUMS microdata extract.

## 1. Inputs required by `run_main.R`

| # | What | Where it must go | Exact naming the code looks for | Source |
|---|------|------------------|---------------------------------|--------|
| 1 | MIGRATE block-group OD flows, 2011-2019 | `data/2011.csv` … `data/2019.csv` | `<year>.csv`, columns incl. origin/destination `geoid_bg` and flow count | Agostini et al., *Nat. Commun.* 2025 (MIGRATE dataset; Zenodo / OSF link in the paper) |
| 2 | WashU ACAG SatPM V6.NA.01 total PM2.5, 2006-2019 | `data/raw/acag_pm25/PM25/` | `*.PM25.<yyyy>01-<yyyy>12.nc` (annual NetCDF, subdataset `PM25`) | https://www.satpm.org/v6-na-01 → Box `wustl.box.com/v/ACAG-V6NA01-CNNPM25` |
| 3 | ACAG speciated components (all eight) | `data/raw/acag_pm25/{SO4,NO3,BC,OM,NH4,DUST,SS}/` | `*.SO4.<yyyy>01-<yyyy>12.nc` etc. (subds = component name) | same Box folder (compositional files) |
| 4 | ACAG PM2.5 absolute uncertainty | `data/raw/acag_pm25_uncertainty/` | `*AbsoluteUncertainty.<yyyy>01-<yyyy>12.nc` (subds `PM25`) | Box `wustl.box.com/v/ACAG-V6NA01-CNNPM25-Uncertainty` |
| 5 | Surface NO2 annual, ~1 km | `data/raw/no2/` | one `*.nc` or `*.tif` per year (no duplicates/tiles in the folder) | Anenberg et al. 2022, figshare 12968114 → GWU Box `gwu.box.com/s/4wey4183hb0nlp3dg4wl02chca3nqnl9` (NetCDF; accepted as .nc or .tif) |
| 6 | Warm-season MDA8 ozone, annual | `data/raw/ozone/annual/` | NetCDF with variable `WSMDA8`, one per year | Liu, Chu, Deziel, Chen 2026 *ES&T* doi:10.1021/acs.est.5c16412 — see paper's Data Availability (Yale) |
| 7 | Wildfire smoke PM2.5 daily 10 km | `data/raw/smoke/smokePM2pt5_predictions_daily_10km_20060101-20231231.rds` + `data/raw/smoke/10km_grid_wgs84/10km_grid_wgs84.shp` | exact names above | Childs et al. 2022 v2.0 (2006-2023) Dropbox, linked from stanfordecholab.com/wildfire_smoke; v1 on Harvard Dataverse doi:10.7910/DVN/DJVMTV |
| 8 | Neighborhood covariate panel | `data/acs_gentrification/blockgroup_gentrification_2009_2019_lag3_with_HMDA.csv` | exact name | built by `upstream/` (needs Census + IPUMS API keys) — see §2 |
| 9 | NHGIS B03002 block-group race tables, ACS5 2007-11 and 2008-12 | `data/nhgis/nhgis*_csv.zip` | zip containing `*_20115_blck_grp.csv` and `*_20125_blck_grp.csv` + codebooks | `upstream/nhgis_download.R` (IPUMS API key) or manual NHGIS extract |
| 10 | Race 2013-2019 | fetched automatically | via `tidycensus` (Census API key: `tidycensus::census_api_key(..., install=TRUE)`) | — |
| 11 | Block-group boundaries (2019 vintage, cb=TRUE) | fetched automatically | via `tigris`, cached in `data/pollution_bg/cache/bg_national.rds` | — |

Override any raw location without moving files:
`PM25_RAW_DIR PM25_UNC_DIR NO2_RAW_DIR OZONE_RAW_DIR SMOKE_PM25_RDS SMOKE_GRID_SHP EPA_GREENBOOK_XLS`.

## 2. Rebuilding the covariate panel (`upstream/`)

Run once, in order, from the repo root:

1. `upstream/nhgis_download.R` — set `IPUMS_API_KEY` in `~/.Renviron`.
2. `upstream/blockgroup_gentrification_earliest_acs.R` — needs Census API key.
3. `upstream/hmda_capital_reinvestment.R` — needs HMDA LAR files: CFPB historic national LAR 2012-2017 (files.consumerfinance.gov) and the FFIEC/CFPB data-browser API for 2018-2019.

## 3. Extra inputs for `run_full.R` only

| What | Where | Source |
|------|-------|--------|
| IPUMS USA ACS 1-yr 2012-2019 migration extract | `data/ipums/usa_acs_migration.csv` (built by `scripts/acs_ipums_pull.R`) | IPUMS API key |
| PUMA ↔ MIGPUMA crosswalk | `data/ipums/puma_migpuma_xwalk.csv` (built by `scripts/build_puma_migpuma_xwalk.R`) | IPUMS USA site |
| EPA Green Book history | `data/raw/epa_greenbook/phistory.xls` | https://www.epa.gov/green-book |

## 4. Verifying the rebuild against the published numbers

`run_all_pm25.R` writes `data/pm25_avoidance_lag3_olag1/logs/manifest_<ts>.csv`
(row counts + md5 per output) and `run_diagnostics.R` hard-fails on known
problems. After the run, compare the headline quantities against the
manuscript (Table 1, Fig. 2-4) and the SI tables; the shift-share Gelbach
additivity check (`gelbach_from_cache.R`) should pass exactly.

## 5. Why this was not run from the Claude session

The sandboxed shells available to Claude (cloud container, and the Linux VM
on this Mac) have no R installation and can only reach github.com, so none of
the data hosts above are reachable from them, and the VM has 3 GB RAM. The
pipeline must be run natively.
