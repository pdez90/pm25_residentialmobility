# Copy the descriptive figures removed from the SI (Figures S22-S60 of the submission)
# into si_data_appendix/ so they are published with the repository.
# Runs descriptive_maps_si.R first if its outputs are missing (needs the IPUMS extract);
# the pollutant sankeys are written by run_main.R (pollution_sankey.R).
if (!exists(".migrate_paths_loaded")) source(if (file.exists("scripts/_paths.R")) "scripts/_paths.R" else "_paths.R")
base_dir <- migrate_data_dir(); out_dir <- file.path(base_dir, "pm25_avoidance_lag3_olag1"); fig_dir <- file.path(out_dir, "figures")
dest <- file.path(repo_root, "si_data_appendix"); dir.create(dest, showWarnings=FALSE)
dm <- file.path(fig_dir, "descriptive_maps")
if (!dir.exists(dm) || !length(list.files(dm, "\\.png$"))) { cat("descriptive maps missing -> running scripts/descriptive_maps_si.R\n"); source(file.path(repo_root, "scripts", "descriptive_maps_si.R")) }
src <- c(list.files(dm, "\\.png$", full.names=TRUE), list.files(fig_dir, "^sankey_.*\\.png$", full.names=TRUE), list.files(file.path(fig_dir, "maps"), "^map_.*\\.png$", full.names=TRUE))
ok <- file.copy(src, file.path(dest, basename(src)), overwrite=TRUE)
cat(sprintf("copied %d of %d figures to %s (%.0f MB)\n", sum(ok), length(src), dest, sum(file.size(file.path(dest, basename(src))))/1e6))
# append a file listing to the README
rd <- file.path(dest, "README.md"); txt <- readLines(rd); txt <- txt[seq_len(if (any(grepl("^## Files", txt))) which(grepl("^## Files", txt))[1]-1 else length(txt))]
writeLines(c(txt, "## Files", "", paste0("- `", sort(basename(src)), "`")), rd)
cat("README updated:", rd, "\n")
