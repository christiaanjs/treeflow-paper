# Summary statistics comparing per-node divergence time estimates between
# TreeFlow's variational inference and BEAST 2 MCMC, written as YAML for
# substitution into the manuscript text.
#
# The figure these accompany (scripts/data-tree-plot.R) shows the per-node
# posterior mean and standard deviation against each other; these are the
# numbers quoted alongside it, computed from the same inputs so the two cannot
# disagree.

library(dplyr)

addHeights <- function(tree) {
    treeDf <- tibble::as_tibble(tree)
    depths <- ape::node.depth.edgelength(ape::as.phylo(tree))
    mutate(treeDf, height = max(depths) - depths)
}

nodeSummary <- function(path) {
    trees <- treeio::read.beast(path)
    bind_rows(lapply(trees, addHeights), .id = "index") %>%
        # Internal nodes only: tips have no height uncertainty (the sampling
        # times are fixed), so including them would inflate the agreement.
        filter(is.na(label)) %>%
        group_by(node) %>%
        summarise(sd = sd(height), mean = mean(height))
}

beast <- nodeSummary(snakemake@input[["beast_tree_samples"]])
vi <- nodeSummary(snakemake@input[["vi_tree_samples"]])
joined <- inner_join(beast, vi, by = "node", suffix = c("_beast", "_vi"))

sdRatio <- joined$sd_vi / joined$sd_beast
quartiles <- quantile(joined$sd_beast, c(0.25, 0.75))
mostPrecise <- joined$sd_beast <= quartiles[1]
leastPrecise <- joined$sd_beast > quartiles[2]

stats <- list(
    flu_internal_node_count = format(nrow(joined), big.mark = ","),
    flu_tree_mean_correlation = sprintf("%.2f", cor(joined$mean_beast, joined$mean_vi)),
    flu_tree_mean_slope = sprintf("%.2f", coef(lm(mean_vi ~ mean_beast, joined))[2]),
    flu_tree_mean_relative_error = sprintf(
        "%.0f\\%%",
        100 * median(abs(joined$mean_vi - joined$mean_beast) / joined$mean_beast)
    ),
    flu_tree_sd_correlation = sprintf("%.2f", cor(joined$sd_beast, joined$sd_vi)),
    flu_tree_sd_ratio_low = sprintf("%.1f", median(sdRatio[mostPrecise])),
    flu_tree_sd_ratio_high = sprintf("%.1f", median(sdRatio[leastPrecise]))
)

# Written directly rather than through the `yaml` package, which is not among
# this project's R dependencies. Single-quoted YAML scalars are used because
# they do not process backslash escapes, and one of the values carries a LaTeX
# `\%`, which is not a valid escape inside a double-quoted scalar. No value
# contains a single quote.
stopifnot(!any(grepl("'", unlist(stats), fixed = TRUE)))
writeLines(
    sprintf("%s: '%s'", names(stats), unlist(stats)),
    snakemake@output[[1]]
)
