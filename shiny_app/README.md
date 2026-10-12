# Results explorer (Shiny)

Interactive companion to the paper. Everything it shows is read from the summary
CSVs in `data/`, which are copies of the small output tables written by
`run_main.R` (and `scripts/reviewer1_urban_form_distance.R`) into
`pm25_avoidance_lag3_olag1/`. No raw data, no large panels, nothing fetched at
runtime.

## Run locally

```r
install.packages(c("shiny","bslib","bsicons","data.table","ggplot2","plotly","DT","leaflet","viridisLite","rsconnect"))
shiny::runApp("shiny_app")
```

## Refresh the data after a pipeline run

```bash
bash shiny_app/refresh_data.sh              # copies every output CSV under 50 KB into shiny_app/data/
Rscript scripts/reviewer1_urban_form_distance.R   # urban-form classes and principal cities, needed by the two scripts below
Rscript scripts/build_app_metro_data.R      # metro explorer: map + Sankey tables for the 30 largest metros (~20 min)
Rscript scripts/build_app_metro_pc.R        # principal-city markers (1 min)
Rscript scripts/build_app_metro_gaps.R      # metro-specific adjusted gaps (2-3 min)
Rscript scripts/build_od_parquet.R          # national origin-destination Parquet for 'Trace a neighbourhood' (~40 min)
OD_ATTR_ONLY=1 Rscript scripts/build_od_parquet.R   # only rebuild data/bg_attr.parquet (income, home value, rent; ~2 min)
```

`build_app_metro_data.R` also writes `metro_transition.csv` (origin PM2.5 quintile x
destination quintile, by origin group and year), which drives the transition matrix on
the metro page; the tab falls back to the Sankey only when that file is absent.

## Shareable links

The query string mirrors the current view (tab, metro, year, origin block group,
affordability rule, highlights), e.g.
`?tab=trace&origin=360850323001&year=2019&rule=4` or `?tab=metro&metro=35620&grp=black&year=2015&view=matrix`,
and each explorer sidebar has a "Link to this view" box.

## Publish to Posit Connect Cloud (free tier is enough)

Connect Cloud deploys straight from a public GitHub repository and needs a
`manifest.json` in the app folder describing the R packages.

1. Generate the manifest once (and again whenever `app.R`'s packages change):
   ```r
   rsconnect::writeManifest(appDir = "shiny_app")
   ```
2. Commit and push `shiny_app/` (including `manifest.json` and `data/`).
   `data/*.csv` is normally git-ignored in this repo; `shiny_app/data/` is
   explicitly un-ignored in `.gitignore` so the app's tables are tracked.
3. Go to https://connect.posit.cloud → **Publish** → **Shiny (R)** → connect
   your GitHub account → choose `pdez90/pm25_residentialmobility`, branch
   `main`, app path `shiny_app/app.R`.
4. Each later `git push` to `main` can redeploy automatically (enable in the
   content settings).

The result is a public URL you can put in the paper's Data Availability
statement next to the repository link.

## Alternative: shinyapps.io

```r
rsconnect::setAccountInfo(name="<account>", token="<token>", secret="<secret>")   # from shinyapps.io dashboard
rsconnect::deployApp("shiny_app", appName="pm25-residential-mobility")
```

## "Trace a neighbourhood" (national origin-destination flows)

The full flows (~23 GB CSV) are not bundled. `scripts/build_od_parquet.R` packs
them as one Parquet file per year (sorted by origin, ~1 GB each, every modelled flow kept) in
`app_od/`, plus a small `shiny_app/data/bg_attr.parquet`. Host the yearly files
anywhere that serves HTTP range requests and point the app at them:

```bash
gh release create od-v1 app_od/od_*.parquet --repo pdez90/pm25_residentialmobility \
   -t "Origin-destination flows (Parquet) for the results explorer"
```

Then in Connect Cloud → the app → Settings → Environment variables add
`OD_PARQUET_BASE = https://github.com/pdez90/pm25_residentialmobility/releases/download/od-v1`
and add `duckdb` and `arrow` to the manifest (`install.packages(c("duckdb","arrow"))`, then
`rsconnect::writeManifest("shiny_app")`). The tab stays hidden until the variable is set.
DuckDB reads only the row groups for the clicked origin (a few MB per click).
