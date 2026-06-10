###############################
### Find the "home" directory.
home <- Sys.getenv("HOME", unset = NA)
if (is.na(home)) {
  if(.Platform$OS.type == "windows") {
    home <- paste(Sys.getenv("HOMEDRIVE"),
                  Sys.getenv("HOMEPATH"), sep = "")
    if (!file.exists(home)) home <- NA
  }
}
if (is.na(home)) stop("Cannot find 'HOME' from environment variables.")

### Make sure there is a "Paths" subdirectory.
pathsdir <- file.path(home, "Paths")
if (!file.exists(pathsdir)) {
  stop("There is no 'Paths' folder inside '", home, "'\n",
       "Please create one.\n")
}

### Read the local JSON file.
jinfo <- file.path(home, "Paths", "context.json")
if (!file.exists(jinfo)) stop("Cannot locate file: '", jinfo, "'.\n", sep='')
library(rjson)
temp <- fromJSON(file = jinfo)
paths <- temp$paths

### Delete any temporary objects, just leving "paths".
rm(home, pathsdir, jinfo, temp)

