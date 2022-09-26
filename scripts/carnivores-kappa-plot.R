library(dplyr)
library(ggplot2)

# treesPath <- "../treeflow/examples/demo-out/carnivores-alt-trees.nexus"
print("Loading trees...")
treesPath <- snakemake@input[["tree_samples"]]
trees <- treeio::read.beast(treesPath)
print("Trees loaded")

addHeights <- function(tree) {
    treeDf <- tibble::as_tibble(tree)
    depths <- ape::node.depth.edgelength(ape::as.phylo(tree))
    withHeights <- mutate(treeDf, height = max(depths) - depths) %>%
        mutate(parentHeight = height[parent])
    withHeights
}

withHeightsDf <- bind_rows(lapply(trees, addHeights), .id = "index")
withMid <- mutate(withHeightsDf, midHeight = parentHeight - height)

loP <- 0.025
upP <- 0.975
summaryFuncs <- list(
    lo = purrr::partial(quantile, probs = loP),
    up = purrr::partial(quantile, probs = upP),
    mean = mean
)

plotDf <- withMid %>%
    group_by(node) %>%
    filter(!is.na(kappa)) %>%
    summarise(across(c(midHeight, kappa), summaryFuncs)) %>%
    rename(Age = midHeight_mean, Kappa = kappa_mean)

fig <- ggplot(plotDf, aes(x = Age, y = Kappa)) +
    geom_errorbar(aes(ymin = kappa_lo, ymax = kappa_up), alpha = 0.5) +
    geom_errorbarh(aes(xmin = midHeight_lo, xmax = midHeight_up), alpha = 0.5) +
    geom_point(color = "navy")

ggplot2::ggsave(snakemake@output[[1]], fig, width = 8, height = 6)
