# Residential mobility and environmental-opportunity inequality

Analysis code for *Residential Mobility Reproduces Environmental-Opportunity Inequality through Housing Markets: Evidence from national migration flows and multi-pollutant exposure, 2011–2019.*

Priyanka N. deSouza — Department of Urban and Regional Planning, University of Colorado Denver; CU Population Center, University of Colorado Boulder.

The study links national block-group origin–destination migration flows (MIGRATE, 2011–2019) to fine-resolution surfaces for eight pollutants, and estimates an **exposure-reallocation estimand**: for the stream of households leaving each neighborhood, the flow-weighted difference between the pollution they moved to and the pollution they left.

---

## What this repository contains

R code only. It does **not** contain input data or model output; every input is public and listed below.

```
scripts/
  run_all_pm25.R              reproducible driver: runs the whole pipeline, stage by stage
  run_diagnostics.R           hard-fail verification of panels, joins, ranges, and outputs
  extract_*.R                 build block-group pollutant panels from the raw exposure surfaces
  rebuild_pollution_panel.R   re-merge cached extractions (fixes FIPS zero-padding)
  idea2_*.R                   core estimation: avoidance gap, sorting, unequal cleanup
  multipollutant_*.R          the same estimators across all eight pollutants
  affordability_robustness.R  Gelbach decomposition, selection and uncertainty robustness
  gelbach_from_cache.R        standalone Gelbach decomposition from the cached origin-year panel
  opportunity_heterogeneity.R environmental opportunity set and segregation interaction
  acs_*.R, build_puma_*.R     individual-level replication in ACS/IPUMS microdata
  build_greenbook.R,
  naaqs_nonattainment.R       exploratory Clean Air Act nonattainment analysis (SI only)
  *_maps_si.R, *_sankey.R,
  conceptual_figure.R,
  dag_figure.R, spec_ladder.R figures for the main text and SI
```

## Running the analysis

```r
Rscript scripts/run_all_pm25.R
```

The driver sets a fixed seed (`20260709`, `L'Ecuyer-CMRG`), tees a timestamped log, asserts the expected outputs of every stage, runs the full diagnostic suite, and writes a manifest with row counts and md5 checksums for every output table. Any failure stops the run with an explicit message; stages you do not want are switched off at the top of the file rather than silently skipped.

Stages, in order: PM2.5 panel → speciated components → multi-pollutant panel → uncertainty panel → core avoidance → figures and specifications → sorting vs. unequal cleanup → cross-pollutant → diagnostic maps and Sankey → affordability and selection robustness → opportunity set and segregation → multi-pollutant robustness → ACS individual-level replication → NAAQS exploratory analysis → diagnostics.

`gelbach_from_cache.R` reproduces the affordability decomposition in seconds from the cached origin-year panel, and self-checks the Gelbach additivity identity (contributions must sum to the total attenuation; verified to machine precision).

### Requirements

R 4.5.0 with `data.table`, `fixest`, `splines`, `ggplot2`, `ggrepel`, `sf`, `terra`, `exactextractr`, `tigris`, `tidycensus`, `stringr`, `readxl`, `ggalluvial`, `foreign`, `ipumsr`. The driver installs anything missing on first run.

## Data sources (all public, none redistributed here)

| Input | Source |
|---|---|
| Block-group origin–destination migration flows | MIGRATE (Agostini et al., *Nat. Commun.* 2025) |
| Total PM2.5 and speciated components (SO4, NO3, BC, OM), with uncertainty | Washington University ACAG V6.NA.01 (Shen et al.; van Donkelaar et al.) |
| Surface NO2 (~1 km) | Anenberg et al. 2022 / Kerr et al. 2021, figshare |
| Warm-season MDA8 ozone | Liu et al., *Environ. Sci. Technol.* |
| Daily wildfire smoke PM2.5 (10 km) | Childs et al., *Environ. Sci. Technol.* 2022 |
| Neighborhood race, income, poverty, education, housing cost | ACS five-year estimates; race from table B03002 |
| Individual movers (race, ethnicity, residence one year earlier) | IPUMS USA, ACS one-year samples 2012–2019 |
| PM2.5 nonattainment designations | EPA Green Book (`phistory.xls`) |

Two inputs must be obtained before the corresponding stages will run: an IPUMS extract (`scripts/acs_ipums_pull.R` builds it via the IPUMS API given an `IPUMS_API_KEY`) and the EPA Green Book history file. The block-group neighborhood covariate panel is assembled upstream of this repository from the ACS five-year series and HMDA.

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

## Notes

Contrasts in the flow analysis are defined by the racial and ethnic composition of the **origin neighborhood**, not the race of individual movers; MIGRATE is not disaggregated by race. The ACS microdata replication addresses this at the individual level, at coarser migration-PUMA resolution.

The nonattainment analysis is exploratory and reported in the SI only: its pre-trends fail, and it is not relied upon.

## License

MIT (see `LICENSE`).
