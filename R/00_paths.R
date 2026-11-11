#### Edit these placeholders before running any script in this folder.
#### `repo` is the project root that contains Enet/, Endpoints/, Review/, Heatmap/, etc.

if (!exists("repo", inherits = FALSE)) {
  repo <- "[path]"
}

if (!exists("biomarker_rdata", inherits = FALSE)) {
  biomarker_rdata <- "[path]"
}

if (!exists("met_names_csv", inherits = FALSE)) {
  met_names_csv <- "[path]"
}
