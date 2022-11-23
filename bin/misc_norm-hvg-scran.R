# Normalization
lsf <- scater::librarySizeFactors(ds)
nn <- scran::buildSNNGraph(ds)
clust <- igraph::cluster_louvain(nn)
table(clust$membership)
sf <- scran::calculateSumFactors(ds, cluster = factor(clust$membership))

plot(lsf, sf, col = clust$membership, log = 'xy')

raw <- ds # TODO: remove
ds <- scuttle::normalizeCounts(ds, size.factors = sf)

# Feature Selection
data <- as.data.frame(scran::modelGeneVar(ds, block = obs$sample))
hvg <- scran::getTopHVGs(data, n = 5000)
data$variable <- rownames(data) %in% hvg

ggplot2::ggplot(data, ggplot2::aes(mean, total, col = variable)) +
  ggplot2::geom_point(shape = 21) +
  ggplot2::geom_smooth(col = "blue") +
  ggplot2::theme_classic(20) +
  ggplot2::scale_x_continuous() +
  ggplot2::scale_y_continuous() +
  ggplot2::coord_fixed(ratio = .5)