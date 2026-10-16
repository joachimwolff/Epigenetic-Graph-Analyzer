# Distance-matched randomisations (ChAseR `randomize(dist.match = TRUE)`) of the step-1 network written as edge lists,
# so that their assortativity can be recomputed per chromosome (analysis/block_bootstrap_results.py). The ChAseR
# assortativity of every randomisation is written alongside as a check of the reimplementation.
#
# Usage: Rscript chaser_random_networks.R <dir with network_step1.txt and features.txt> <output dir> [randomisations=10] [seed=1]

suppressPackageStartupMessages(library(chaser))

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
out <- args[2]
n_random <- if (length(args) > 2) as.integer(args[3]) else 10
seed <- if (length(args) > 3) as.integer(args[4]) else 1
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)
net <- chaser::make_chromnet(file.path(dir, "network_step1.txt"))
net <- chaser::load_features(net, file.path(dir, "features.txt"), type = "features_on_nodes", missingv = 0)
random <- chaser::randomize(net, nrandom = n_random, dist.match = TRUE)
rows <- list()
for (i in seq_along(random)) {
  edges <- random[[i]]$edgesdf[, 1:6]
  write.table(edges, file.path(out, sprintf("random_step1_%02d.txt", i)), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  values <- chaser::chas(random[[i]])
  rows[[i]] <- data.frame(randomisation = i, mark = names(values), assortativity = as.numeric(values))
}
write.table(do.call(rbind, rows), file.path(out, "random_step1_chas.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
