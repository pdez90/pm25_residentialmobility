# Residential mobility and environmental-opportunity inequality

Analysis code for *Residential Mobility Reproduces Environmental-Opportunity Inequality through Housing Markets: Evidence from national migration flows and multi-pollutant exposure, 2011–2019.*

Priyanka N. deSouza — Department of Urban and Regional Planning, University of Colorado Denver; CU Population Center, University of Colorado Boulder.

The study links national block-group origin–destination migration flows (MIGRATE, 2011–2019) to fine-resolution surfaces for eight pollutants and estimates an **exposure-reallocation estimand**: for the stream of households leaving each neighborhood, the flow-weighted difference between the pollution they moved to and the pollution they left.

---

## Quick start

```bash
git clone https://github.com/pdez90/pm25_residentialmobility.git
cd pm25_residentialmobility

export MIGRATE_DATA_DIR=/path/to/your/data    # see "Data layout" below
Rscript setup/install_pinned.R                # or: R -e 'renv::restore()'
Rscript run_main.R
```

`run_main.R` reproduces the paper pipeline and needs no credentials. `run_full.R` additionally runs the individual-level ACS replication (needs an IPUMS extract) and the exploratory nonattainment analysis (needs the EPA Green Book file); both fail loudly if their inputs are absent rather than silently skipping.

Run scripts from the repository root. There are **no absolute paths anywhere in this repository**: every location is resolved in `scripts/_paths.R` from `MIGRATE_DATA_DIR` and, where relevant, per-dataset environment variables.

## Data layout

`MIGRATE_DATA_DIR` is the working root. The pipeline reads and writes:

```
$MIGRATE_DATA_DIR/
  <year>.csv                  MIGRATE block-group flow files, 2011-2019
  acs_gentrification/         neighborhood covariate panel (see "Upstream inputs")
  pollution_bg/               derived block-group pollutant panels (created by the pipeline)
  pm25_avoidance_lag3_olag1/  all analysis outputs, figures, logs, manifests
  raw/                        raw exposure surfaces (or point to them individually)
    acag_pm25/                total PM2.5 and speciated components
    acag_pm25_uncertainty/    paired uncertainty product
    no2/                      surface NO2
    ozone/annual/             warm-season MDA8 ozone
    smoke/                    wildfire smoke PM2.5 predictions + 10 km grid shapefile
    epa_greenbook/            phistory.xls (run_full.R only)
```

Raw surfaces are large and often stored elsewhere. Override any of them individually:
`PM25_RAW_DIR`, `PM25_UNC_DIR`, `NO2_RAW_DIR`, `OZONE_RAW_DIR`, `SMOKE_PM25_RDS`, `SMOKE_GRID_SHP`, `EPA_GREENBOOK_XLS`.

## Data sources

All inputs are public; none are redistributed here.

| Input | Source |
|---|---|
| Block-group origin–destination migration flows | MIGRATE (Agostini et al., *Nat. Commun.* 2025) |
| Total PM2.5 and components (SO4, NO3, BC, OM) with uncertainty | Washington University ACAG V6.NA.01 (Shen et al.; van Donkelaar et al.) |
| Surface NO2 (~1 km) | Anenberg et al. 2022 / Kerr et al. 2021, figshare |
| Warm-season MDA8 ozone | Liu et al., *Environ. Sci. Technol.* |
| Daily wildfire smoke PM2.5 (10 km) | Childs et al., *Environ. Sci. Technol.* 2022 |
| ACS five-year block-group tables (income, poverty, education, housing cost, race B03002) | Census API via `tidycensus`; earlier years via NHGIS/IPUMS |
| HMDA loan records | CFPB historic national LAR (2012–2017); FFIEC/CFPB data browser API (2018–2019) |
| Individual movers (race, ethnicity, residence one year earlier) | IPUMS USA, ACS one-year samples 2012–2019 |
| PM2.5 nonattainment designations | EPA Green Book (`phistory.xls`) |

## Upstream inputs

The analysis reads a prepared block-group covariate panel,
`acs_gentrification/blockgroup_gentrification_2009_2019_lag3_with_HMDA.csv`.
The code that builds it from the raw public sources is included in `upstream/`, so the chain is reproducible end to end:

```
upstream/nhgis_download.R                        ACS block-group tables via the IPUMS/NHGIS API
upstream/blockgroup_gentrification_earliest_acs.R  neighborhood change + gentrification measures
upstream/hmda_capital_reinvestment.R             adds the HMDA capital-reinvestment measure
```

These three scripts are shared with a separate project on neighborhood change and are reproduced here for provenance. They require an IPUMS API key and a Census API key, and are not invoked by `run_main.R`; run them once to build the panel, or obtain the panel directly.

## Reproducibility

`scripts/run_all_pm25.R` is the engine behind both drivers. It sets a fixed seed (`20260709`, `L'Ecuyer-CMRG`), tees a timestamped log, asserts the expected outputs of every stage, runs the diagnostic suite, and writes a manifest with row counts and md5 checksums for every output table. Any failure stops the run with an explicit message.

Stages, in order: PM2.5 panel → speciated components → multi-pollutant panel → uncertainty panel → core avoidance → figures and specifications → sorting vs. unequal cleanup → cross-pollutant → diagnostic maps and Sankey → affordability and selection robustness → opportunity set and segregation → multi-pollutant robustness → [ACS replication] → [NAAQS exploratory] → diagnostics.

`scripts/run_diagnostics.R` hard-fails on duplicate block-group-year rows, coverage below threshold in any year, disagreement between the two PM2.5 sources, flow-join or race-join shortfalls, implausible pollutant ranges, missing outputs, and violation of the Gelbach additivity identity.

**Environment.** `renv.lock` and `setup/install_pinned.R` pin the package versions recorded in the run log that produced the published results (R 4.5.0). The lockfile was reconstructed from that log rather than captured by `renv::snapshot()`, so it carries versions and sources but not hashes; regenerating it with `renv::snapshot()` after a successful restore will add them.

`gelbach_from_cache.R` reproduces the affordability decomposition in seconds from the cached origin-year panel and self-checks that the per-covariate contributions sum exactly to the total attenuation.

## Where each result comes from

| Paper element | Script |
|---|---|
| Conceptual framework figure | `conceptual_figure.R` |
| Exposure trends and relative disparities | `idea2_pm25_avoidance_v3.R` |
| Shift-share decomposition (cleanup vs. sorting) | `idea2_pm25_avoidance_v3.R`; across pollutants, `multipollutant_disparity.R` |
| Conditional avoidance gap | `idea2_pm25_avoidance_v3.R` |
| Three mechanisms across eight pollutants | `idea2_pollution_crosscut.R` |
| Sorting and unequal cleanup, robustness | `pm25_sorting_vs_siting.R` |
| Gelbach affordability decomposition, selection and uncertainty checks | `affordability_robustness.R`, `gelbach_from_cache.R` |
| Opportunity set and segregation gradient | `opportunity_heterogeneity.R` |
| Individual-level replication | `acs_ipums_pull.R` → `build_puma_migpuma_xwalk.R` → `acs_micro_validation.R` |
| Specification ladder | `spec_ladder.R` |
| Causal-structure diagram | `dag_figure.R` |
| Block-group area vs. grid resolution | `blockgroup_area_summary.R` |
| Descriptive maps and migration Sankey diagrams | `descriptive_maps_si.R`, `national_maps_si.R`, `pollution_sankey.R` |
| Exploratory nonattainment analysis | `build_greenbook.R` → `naaqs_nonattainment.R` |

`idea2_pm25_avoidance.R` builds the origin-year reallocation panel and subgroup summaries that `idea2_pm25_avoidance_v3.R` and later stages consume; both are required, and the driver runs them in that order.

## Interpretation notes

Contrasts in the flow analysis are defined by the racial and ethnic composition of the **origin neighborhood**, not the race of individual movers: MIGRATE is not disaggregated by race. Statements from that analysis are of the form "streams leaving majority-Black neighborhoods reached…". Individual-level racial statements rest on the ACS microdata replication, at coarser migration-PUMA resolution.

Destination housing cost is a post-move, chosen quantity. The Gelbach decomposition is therefore a descriptive attribution, not causal mediation.

The uncertainty-weighted specification is a sensitivity check, not a formal measurement-error correction: the ACAG uncertainty product is not an ordinary independent sampling standard error.

The nonattainment analysis is exploratory and reported in the SI only; its pre-trends fail and it is not relied upon.

## License

MIT (see `LICENSE`).
