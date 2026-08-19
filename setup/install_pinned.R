# ============================================================
# Install the exact package versions used for the published results.
# Versions are those recorded by the driver in the run log for the
# analysis reported in the paper (R 4.5.0, macOS aarch64).
#
#   Rscript setup/install_pinned.R
#
# Preferred alternative, if you use renv:
#   renv::restore()        # uses renv.lock
# ============================================================
pkgs <- c(
  data.table    = "1.18.2.1",
  stringr       = "1.6.0",
  fixest        = "0.13.2",
  ggplot2       = "4.0.3",
  ggrepel       = "0.9.6",
  sf            = "1.0.21",
  terra         = "1.8.70",
  exactextractr = "0.10.0",
  tigris        = "2.2.1",
  tidycensus    = "1.7.1",
  readxl        = "1.4.5",
  ggalluvial    = "0.12.6",
  foreign       = "0.8.90"
)
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")
for (p in names(pkgs)) {
  cur <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  if (identical(cur, unname(pkgs[[p]]))) { message(sprintf("%-14s %s already installed", p, pkgs[[p]])); next }
  message(sprintf("%-14s installing %s (have %s)", p, pkgs[[p]], ifelse(is.na(cur), "none", cur)))
  remotes::install_version(p, version = pkgs[[p]], repos = "https://cloud.r-project.org", upgrade = "never")
}
message("\nNote: `splines` ships with R (4.5.0). `here` is optional but recommended for path resolution.")
