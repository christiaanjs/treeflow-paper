library(dplyr)
library(ggplot2)

addHeights <- function(tree) {
    treeDf <- tibble::as_tibble(tree)
    depths <- ape::node.depth.edgelength(ape::as.phylo(tree))
    withHeights <- mutate(treeDf, height = max(depths) - depths)
    withHeights
}

getHeightDf <- function(trees) {
    bind_rows(lapply(trees, addHeights), .id = "index")
}

print("Loading beast trees...")
beastTreesPath <- snakemake@input[["beast_tree_samples"]]
beastTrees <- treeio::read.beast(beastTreesPath)
beastHeights <- getHeightDf(beastTrees)

# Supports either a single `vi_tree_samples` input (the generic per-dataset
# comparison, one VI method) or one or more `vi_tree_samples_<approx>` inputs,
# each compared against the same BEAST 2 trees. Only the inputs the invoking
# rule actually declares are used.
viMethodInputs <- list(
    `Treeflow VI (root full rank)` = "vi_tree_samples_root_full_rank",
    `Treeflow VI (full rank)` = "vi_tree_samples_full_rank",
    `Treeflow VI (mean field)` = "vi_tree_samples_mean_field",
    `Treeflow VI` = "vi_tree_samples"
)

viHeightsByMethod <- Filter(Negate(is.null), lapply(viMethodInputs, function(inputKey) {
    path <- snakemake@input[[inputKey]]
    if (is.null(path)) {
        return(NULL)
    }
    print(paste("Loading variational trees:", inputKey))
    getHeightDf(treeio::read.beast(path))
}))
stopifnot(length(viHeightsByMethod) > 0)

print("Done")

heightsDf <- tibble::as_tibble(bind_rows(
    c(list(`BEAST 2` = beastHeights), viHeightsByMethod),
    .id = "method"
))

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
    ) %>%
    tidyr::pivot_longer(
        cols = tidyselect::any_of(names(viHeightsByMethod)),
        names_to = "Approximation",
        values_to = "Treeflow VI"
    )

fig <- ggplot(plotDf, aes(x = `BEAST 2`, y = `Treeflow VI`, colour = Approximation)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    geom_point(alpha = 0.6) +
    geom_blank(data = limits, aes(x = value, y = value, colour = NULL)) +
    ylab("TreeFlow VI") +
    facet_wrap(~Statistic, scales = "free")

if (length(viHeightsByMethod) == 1) {
    # Nothing to distinguish -- a one-entry "Approximation" legend is just noise
    fig <- fig + ggplot2::guides(colour = "none")
}

outputFile <- snakemake@output[[1]]
ggplot2::ggsave(outputFile, fig, width = 8, height = (4 * 11.7 / 12.5))
