# Load dataset paths. Set CONTEXT_PATHS to a JSON file, or use the default
# ~/Paths/context.json. See data/context-paths.example.json.
context_config_file <- Sys.getenv("CONTEXT_PATHS", unset = path.expand("~/Paths/context.json"))
if (!file.exists(context_config_file)) {
  stop("Path configuration not found: ", context_config_file,
       ". Copy data/context-paths.example.json, edit it, and set CONTEXT_PATHS.", call. = FALSE)
}
if (!requireNamespace("rjson", quietly = TRUE)) {
  stop("Install the configuration reader with install.packages('rjson').", call. = FALSE)
}
paths <- rjson::fromJSON(file = context_config_file)$paths
if (!is.list(paths) || !all(c("mirtcga", "scratch", "results") %in% names(paths))) {
  stop("The configuration must contain paths.mirtcga, paths.scratch and paths.results.", call. = FALSE)
}
for (context_path_key in c("mirtcga", "scratch", "results")) {
  if (!is.character(paths[[context_path_key]]) || length(paths[[context_path_key]]) != 1L ||
      is.na(paths[[context_path_key]]) || !nzchar(paths[[context_path_key]])) {
    stop("Invalid path configuration value: ", context_path_key, call. = FALSE)
  }
}
rm(context_config_file, context_path_key)
