###############################
### Set the root directory, which is where we can find the TCGA
### data on the local machine.

home <- Sys.getenv("HOME", unset = NA)
if (is.na(home)) stop("Cannot find 'HOME' from environment variable s.")

jinfo <- file.path(home, "Paths", "tcga.json")
if (!file.exists(jinfo)) stop("Cannot locate file: '", jinfo, "'.\n", sep='')
library(rjson)
temp <- fromJSON(file = jinfo)
tcgapaths <- temp$paths

rm(home, jinfo, temp)
