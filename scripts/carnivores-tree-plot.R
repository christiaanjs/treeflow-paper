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

heightsDf <- bind_rows(Base = getHeightDf(baseTrees), `Kappa variation` = getHeightDf(altTrees), .id = "model")

loP <- 0.025
upP <- 0.975
summaryFuncs <- list(
    lo = purrr::partial(quantile, probs = loP),
    up = purrr::partial(quantile, probs = upP),
    mean = mean
)

summaryDf <- heightsDf %>%
    filter(parentHeight != height) %>%
    group_by(node, model) %>%
    summarise(
        lo = quantile(parentHeight, probs = loP),
        up = quantile(parentHeight, probs = upP),
        Age = mean(parentHeight)
    )

plotDf <- tidyr::pivot_wider(
    summaryDf,
    names_from = model,
    values_from = c(lo, up, Age),
    names_sep = " "
)

fig <- ggplot(plotDf, aes(x = `Age Base`, y = `Age Kappa variation`)) +
    geom_point(color = "navy") +
    geom_errorbar(aes(ymin = `lo Kappa variation`, ymax = `up Kappa variation`), alpha = 0.5) +
    geom_errorbarh(aes(xmin = `lo Base`, xmax = `up Base`), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted")


ggplot2::ggsave(snakemake@output[[1]], fig, width = 8, height = 6)
