# Chromatin assortativity (ChAseR, Pancaldi et al. 2016; Madrid-Mencia et al. 2020) of fragment pairs at graph
# distance 1 to 5, and of distance-matched randomised networks of the direct interactions.
#
# Usage: Rscript chaser_assortativity.R <dir with network_step*.txt and features.txt> <output.tsv> [randomisations=10]
#        [seed=1]

suppressPackageStartupMessages(library(chaser))

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
output <- args[2]
n_random <- if (length(args) > 2) as.integer(args[3]) else 10
# the randomisations were unseeded, so their assortativity changed between runs
set.seed(if (length(args) > 3) as.integer(args[4]) else 1)
features <- file.path(dir, "features.txt")
rows <- list()
for (step in 1:5) {
  net <- chaser::make_chromnet(file.path(dir, sprintf("network_step%d.txt", step)))
  net <- chaser::load_features(net, features, type = "features_on_nodes", missingv = 0)
  values <- chaser::chas(net)
  rows[[length(rows) + 1]] <- data.frame(network = sprintf("step %d", step), mark = names(values),
                                         assortativity = as.numeric(values))
  if (step == 1) {
    random <- chaser::randomize(net, nrandom = n_random, dist.match = TRUE)
    random_values <- sapply(random, function(r) chaser::chas(r))
    rows[[length(rows) + 1]] <- data.frame(network = "step 1, distance-matched random", mark = rownames(random_values),
                                           assortativity = rowMeans(random_values),
                                           sd = apply(random_values, 1, sd))
  }
}
result <- do.call(rbind, lapply(rows, function(r) { if (!"sd" %in% names(r)) r$sd <- NA; r }))
write.table(result, output, sep = "\t", quote = FALSE, row.names = FALSE)
print(result)
