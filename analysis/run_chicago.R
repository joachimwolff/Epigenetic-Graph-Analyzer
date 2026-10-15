# Standard CHiCAGO run (default settings): merge the chinput files of all libraries, call interactions and export
# those with score >= 5 (or --cutoff=<x>) as ibed and WashU text, plus the full CHiCAGO object. Replicates that are
# combined later by EGA's replicate rule (present in all replicates, called in at least one) are exported with a
# lower cutoff, e.g. --cutoff=1, so that sub-threshold scores of a pair remain visible.
#
# Usage: Rscript run_chicago.R [--cutoff=<x>] <design dir> <output prefix> <chinput> [<chinput> ...]

suppressPackageStartupMessages(library(Chicago))

args <- commandArgs(trailingOnly = TRUE)
cutoff <- 5
if (any(grepl("^--cutoff=", args))) {
  cutoff <- as.numeric(sub("^--cutoff=", "", args[grepl("^--cutoff=", args)][1]))
  args <- args[!grepl("^--cutoff=", args)]
}
design_dir <- args[1]
prefix <- args[2]
chinputs <- args[-(1:2)]

cd <- setExperiment(designDir = design_dir)
cd <- readAndMerge(files = chinputs, cd = cd)
cd <- chicagoPipeline(cd, outprefix = prefix, printMemory = FALSE)
exportResults(cd, outfileprefix = prefix, cutoff = cutoff, format = c("interBed", "washU_text"))
saveRDS(cd, paste0(prefix, ".Rds"))
