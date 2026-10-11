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
Rscript scripts/build_app_metro_data.R      # metro explorer: map + Sankey tables for the 30 largest metros (~20 min)
```

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
