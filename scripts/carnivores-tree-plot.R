library(dplyr)
library(ggplot2)


# altTreesPath <- "../treeflow/examples/demo-out/carnivores-alt-trees.nexus"
# baseTreesPath <- "../treeflow/examples/demo-out/carnivores-base-trees.nexus"

print("Loading base trees...")
baseTreesPath <- snakemake@input[["base_tree_samples"]]
baseTrees <- treeio::read.beast(baseTreesPath)

print("Loading alt trees...")
altTreesPath <- snakemake@input[["alt_tree_samples"]]
altTrees <- treeio::read.beast(altTreesPath)

print("Done")

addHeights <- function(tree) {
    treeDf <- tibble::as_tibble(tree)
    depths <- ape::node.depth.edgelength(ape::as.phylo(tree))
    withHeights <- mutate(treeDf, height = max(depths) - depths) %>%
        mutate(parentHeight = height[parent])
    withHeights
}

getHeightDf <- function(trees) {
    bind_rows(lapply(trees, addHeights), .id = "index")
}

heightsDf <- tibble::as_tibble(bind_rows(Base = getHeightDf(baseTrees), `Kappa variation` = getHeightDf(altTrees), .id = "model"))

loP <- 0.025
upP <- 0.975
summaryFuncs <- list(
    lo = purrr::partial(quantile, probs = loP),
    up = purrr::partial(quantile, probs = upP),
    mean = mean
)

# Summarise each internal node's own height. Summarising the children's
# `parentHeight` instead emits every internal node twice (once per child), which
# duplicated every point and error bar in the plot. In this representation tips
# are the nodes that never appear in the `parent` column; every other node
# (including the root, which is its own parent) is internal.
internalNodes <- unique(heightsDf$parent)

summaryDf <- heightsDf %>%
    filter(node %in% internalNodes) %>%
    group_by(node, model) %>%
    summarise(
        lo = quantile(height, probs = loP),
        up = quantile(height, probs = upP),
        Age = mean(height)
    )
maxAge <- max(summaryDf$up)
plotDf <- tidyr::pivot_wider(
    summaryDf,
    names_from = model,
    values_from = c(lo, up, Age),
    names_sep = " "
) %>%
    dplyr::rename(
        `Age in base model` = `Age Base`,
        `Age in kappa variation model` = `Age Kappa variation`
    )

# One row per internal node: 62 taxa => 61 internal nodes.
print(paste("Internal nodes per model:", nrow(plotDf)))
stopifnot(nrow(plotDf) == length(internalNodes))

fig <- ggplot(plotDf, aes(x = `Age in base model`, y = `Age in kappa variation model`)) +
    geom_errorbar(aes(ymin = `lo Kappa variation`, ymax = `up Kappa variation`, color = "95% posterior quantile interval"), alpha = 0.5) +
    geom_errorbarh(aes(xmin = `lo Base`, xmax = `up Base`, color = "95% posterior quantile interval"), alpha = 0.5) +
    geom_point(aes(color = "Posterior mean")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    expand_limits(x = maxAge, y = maxAge) +
    scale_color_manual(name = "", values = c("black", "navy")) +
    guides(colour = guide_legend(override.aes = list(
        linetype = c("solid", "blank"),
        shape = c(NA, 19)
    ))) +
    theme(legend.position = "bottom")


ggplot2::ggsave(snakemake@output[[1]], fig, width = 4.5, height = 5)
