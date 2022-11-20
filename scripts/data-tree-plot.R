library(dplyr)
library(ggplot2)

print("Loading beast trees...")
# beastTreesPath <- "out/dengue/beast.trees"
beastTreesPath <- snakemake@input[["beast_tree_samples"]]
beastTrees <- treeio::read.beast(beastTreesPath)

print("Loading variational trees...")
# variationalTreesPath <- "out/dengue/variational-tree-samples.nexus"
variationalTreesPath <- snakemake@input[["vi_tree_samples"]]
variationalTrees <- treeio::read.beast(variationalTreesPath)

print("Done")

addHeights <- function(tree) {
    treeDf <- tibble::as_tibble(tree)
    depths <- ape::node.depth.edgelength(ape::as.phylo(tree))
    withHeights <- mutate(treeDf, height = max(depths) - depths)
    withHeights
}

getHeightDf <- function(trees) {
    bind_rows(lapply(trees, addHeights), .id = "index")
}

heightsDf <- bind_rows(`BEAST 2` = getHeightDf(beastTrees), `TreeFlow VI` = getHeightDf(variationalTrees), .id = "method")



summaryDf <- heightsDf %>%
    filter(is.na(label)) %>%
    group_by(node, method) %>%
    summarise(
        `Standard Deviation` = sd(height),
        Mean = mean(height)
    )
longForm <- tidyr::pivot_longer(
    summaryDf,
    c(Mean, `Standard Deviation`),
    names_to = "Statistic"
)
limits <- longForm %>%
    group_by(Statistic) %>%
    summarise(min = min(value), max = max(value)) %>%
    tidyr::pivot_longer(c(min, max))
plotDf <- longForm %>%
    tidyr::pivot_wider(
        names_from = method,
        values_from = value
    )

fig <- ggplot(plotDf, aes(x = `BEAST 2`, y = `TreeFlow VI`)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    geom_point() +
    geom_blank(data = limits, aes(x = value, y = value)) +
    facet_wrap(~Statistic, scales = "free")

# outputFile <- "manuscript/figures/flu-tree-plot.png"
outputFile <- snakemake@output[[1]]
ggplot2::ggsave(outputFile, fig, width = 7, height = (4 * 11.7 / 12.5))
