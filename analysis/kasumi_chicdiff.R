# Chicdiff (Cairns et al. 2019) on the Kasumi-1 capture Hi-C replicates called with the CHiCAGO pipeline of
# HiCExplorer v4 (kasumi_call_chicago.sh). Chicdiff needs, per replicate, the chinput counts and a CHiCAGO object
# with the per-pair background (Bmean, Tmean, score), bait factors s_j and bait pools tblb, other-end factors s_i and
# pools tlb, the distance bins (distbin, refBinMean) and the dispersion; it also accepts a data.table with an attribute
# "dispersion" in place of a chicagoData object. These are assembled from the per-replicate <name>_scores.txt and
# <name>_background.txt of HiCExplorer: s_j [BaitFactors], tblb [BaitPools], tlb [OtherEndPools], s_i of the pool
# [PoolFactors], dispersion = alpha [Dispersion] (NB size, CHiCAGO's theta), and refBinMean of every 20 kb distance bin
# up to maxLBrownEst from the fitted cubic [DistFun], so that Chicdiff's refit of the distance function returns the
# same cubic. The peak matrix holds every pair with score >= 5 in at least one replicate and the score of the pair
# in each replicate (0 where the pair is not in the score >= 1 export).
#
# Usage:
#   Rscript kasumi_chicdiff.R prepare <scores.txt> <background.txt> <out.Rds>
#   Rscript kasumi_chicdiff.R peaks <rmap> <baitmap> <out.txt> name=<score1.txt> [name=<score1.txt> ...]
#   Rscript kasumi_chicdiff.R run <designDir> <outprefix> <peakmatrix> <condA>:<name>=<chinput>,<Rds>[;...] <condB>:...

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
command <- args[1]

read_section <- function(lines, name) {
  start <- which(lines == paste0("[", name, "]"))
  if (length(start) != 1) stop("section ", name, " not found")
  ends <- which(startsWith(lines, "[") & seq_along(lines) > start)
  end <- if (length(ends)) min(ends) - 1 else length(lines)
  body <- lines[(start + 1):end]
  body <- body[nchar(body) > 0]
  fread(text = paste(body, collapse = "\n"), sep = "\t")
}

key_values <- function(lines, name) {
  start <- which(lines == paste0("[", name, "]"))
  ends <- which(startsWith(lines, "[") & seq_along(lines) > start)
  body <- lines[(start + 1):(min(c(ends, length(lines) + 1)) - 1)]
  body <- body[nchar(body) > 0]
  parts <- strsplit(body, "\t")
  setNames(as.numeric(sapply(parts, `[`, 2)), sapply(parts, `[`, 1))
}

if (command == "prepare") {
  lines <- readLines(args[3])
  settings <- key_values(lines, "Settings")
  fun <- key_values(lines, "DistFun")
  alpha <- key_values(lines, "Dispersion")[["alpha"]]
  bait_factors <- read_section(lines, "BaitFactors")
  bait_pools <- read_section(lines, "BaitPools")
  oe_pools <- read_section(lines, "OtherEndPools")
  pool_factors <- read_section(lines, "PoolFactors")
  x <- fread(args[2], select = c("baitID", "otherEndID", "N", "distSign", "Bmean", "Tmean", "score"))
  x[bait_factors, s_j := i.s_j, on = "baitID"]
  x[bait_pools, tblb := i.tblb, on = "baitID"]
  setnames(oe_pools, "pool", "tlb")
  oe_pools[pool_factors, s_i := i.s_i, on = c(tlb = "pool")]
  x[oe_pools, `:=`(tlb = i.tlb, s_i = i.s_i), on = "otherEndID"]
  # distance bins: refBinMean of bin b = f(midpoint_b) of the cubic fit; one row per bin carries it
  binsize <- settings[["binsize"]]
  n_bins <- floor(settings[["maxLBrownEst"]] / binsize)
  midpoint <- seq(from = round(binsize / 2), by = binsize, length.out = n_bins)
  l <- log(midpoint)
  ref <- exp(fun[["cubic0"]] + fun[["cubic1"]] * l + fun[["cubic2"]] * l^2 + fun[["cubic3"]] * l^3)
  if (abs(log(min(midpoint)) - fun[["obsMinLog"]]) > 1e-6 || abs(log(max(midpoint)) - fun[["obsMaxLog"]]) > 1e-6)
    stop("distance bins do not match the observed range of the fit")
  if (is.unsorted(rev(ref))) stop("refBinMean not decreasing with distance")
  x[, distbin := NA_integer_]
  x[, refBinMean := NA_real_]
  x[seq_len(n_bins), `:=`(distbin = seq_len(n_bins), refBinMean = ref)]
  setattr(x, "dispersion", alpha)
  message(sprintf("%s: %d pairs, %d baits with s_j, %d with tblb, other ends with tlb %.3f, dispersion %.4f",
                  args[2], nrow(x), sum(!is.na(bait_factors$s_j)), nrow(bait_pools), mean(!is.na(x$tlb)), alpha))
  saveRDS(x, args[4], compress = FALSE)
} else if (command == "peaks") {
  rmap <- fread(args[2], col.names = c("chr", "start", "end", "ID"))
  baitmap <- fread(args[3], col.names = c("chr", "start", "end", "ID", "name"))
  out <- args[4]
  samples <- args[-(1:4)]
  tables <- list()
  for (s in samples) {
    name <- sub("=.*", "", s)
    t <- fread(sub("^[^=]*=", "", s), select = c("bait_chr", "bait_start", "bait_end", "otherEnd_chr",
                                                 "otherEnd_start", "otherEnd_end", "score"))
    t[rmap, baitID := i.ID, on = c(bait_chr = "chr", bait_start = "start", bait_end = "end")]
    t[rmap, oeID := i.ID, on = c(otherEnd_chr = "chr", otherEnd_start = "start", otherEnd_end = "end")]
    if (anyNA(t$baitID) || anyNA(t$oeID)) stop(s, ": fragments not in the rmap")
    t <- t[, .(score = max(score)), by = .(baitID, oeID)]
    setnames(t, "score", name)
    tables[[name]] <- t
  }
  m <- Reduce(function(a, b) merge(a, b, by = c("baitID", "oeID"), all = TRUE), tables)
  for (name in names(tables)) set(m, which(is.na(m[[name]])), name, 0)
  m <- m[do.call(pmax, m[, names(tables), with = FALSE]) >= 5]
  m[rmap, `:=`(baitChr = i.chr, baitStart = i.start, baitEnd = i.end), on = c(baitID = "ID")]
  m[rmap, `:=`(oeChr = i.chr, oeStart = i.start, oeEnd = i.end), on = c(oeID = "ID")]
  m[baitmap, baitName := i.name, on = c(baitID = "ID")]
  m[, oeName := "."]
  m[, dist := ifelse(baitChr == oeChr, round((oeStart + oeEnd) / 2 - (baitStart + baitEnd) / 2), NA)]
  setcolorder(m, c("baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr", "oeStart", "oeEnd", "oeID",
                   "oeName", "dist", names(tables)))
  message(sprintf("peak matrix: %d pairs with score >= 5 in at least one of %d replicates", nrow(m), length(tables)))
  fwrite(m, out, sep = "\t")
} else if (command == "run") {
  suppressPackageStartupMessages({library(Chicago); library(Chicdiff)})
  design_dir <- normalizePath(args[2])
  peaks <- normalizePath(args[4])
  # Chicdiff writes its plots to "./<outprefix>...": run in the output directory with a relative prefix
  setwd(dirname(args[3]))
  prefix <- basename(args[3])
  chicago_data <- list()
  count_data <- list()
  for (spec in args[-(1:4)]) {
    condition <- sub(":.*", "", spec)
    entries <- strsplit(sub("^[^:]*:", "", spec), ";")[[1]]
    names_ <- sub("=.*", "", entries)
    files <- strsplit(sub("^[^=]*=", "", entries), ",")
    count_data[[condition]] <- setNames(sapply(files, `[`, 1), names_)
    chicago_data[[condition]] <- setNames(sapply(files, `[`, 2), names_)
  }
  settings <- setChicdiffExperiment(designDir = design_dir, chicagoData = chicago_data, countData = count_data,
                                    peakfiles = peaks, outprefix = prefix, settings = list(printMemory = TRUE))
  start <- Sys.time()
  output <- chicdiffPipeline(settings)
  message("Chicdiff pipeline: ", format(Sys.time() - start))
  fwrite(as.data.table(output), paste0(prefix, "_results.tsv.gz"), sep = "\t")
} else {
  stop("unknown command ", command)
}
