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

# Batch means SE: split chain into blocks, compute mean per block, take SD across blocks
n_blocks <- 10
batchMeanSE <- function(x) {
    n <- length(x)
    block_size <- floor(n / n_blocks)
    if (block_size < 2) return(sd(x) / sqrt(n))
    block_means <- sapply(seq_len(n_blocks), function(i) {
        start <- (i - 1) * block_size + 1
        end <- min(i * block_size, n)
        mean(x[start:end])
    })
    sd(block_means) / sqrt(n_blocks)
}

plotDf <- withMid %>%
    group_by(node) %>%
    filter(!is.na(kappa)) %>%
    summarise(
        across(c(midHeight, kappa), summaryFuncs),
        kappa_mcse = batchMeanSE(kappa),
        midHeight_mcse = batchMeanSE(midHeight)
    ) %>%
    rename(Age = midHeight_mean, Kappa = kappa_mean) %>%
    mutate(
        kappa_mc_lo = Kappa - 2 * kappa_mcse,
        kappa_mc_up = Kappa + 2 * kappa_mcse,
        age_mc_lo = Age - 2 * midHeight_mcse,
        age_mc_up = Age + 2 * midHeight_mcse
    )

fig <- ggplot(plotDf, aes(x = Age, y = Kappa)) +
    geom_errorbar(aes(ymin = kappa_lo, ymax = kappa_up, color = "95% posterior quantile interval"), alpha = 0.5) +
    geom_errorbarh(aes(xmin = midHeight_lo, xmax = midHeight_up, color = "95% posterior quantile interval"), alpha = 0.5) +
    geom_errorbar(aes(ymin = kappa_mc_lo, ymax = kappa_mc_up, color = "Monte Carlo SE (\u00b12)"), width = 0, linewidth = 0.8) +
    geom_errorbarh(aes(xmin = age_mc_lo, xmax = age_mc_up, color = "Monte Carlo SE (\u00b12)"), height = 0, linewidth = 0.8) +
    geom_point(aes(color = "Posterior mean")) +
    scale_color_manual(
        name = "",
        values = c("95% posterior quantile interval" = "black",
                   "Monte Carlo SE (\u00b12)" = "firebrick",
                   "Posterior mean" = "navy")
    ) +
    guides(colour = guide_legend(override.aes = list(
        linetype = c("solid", "solid", "blank"),
        shape = c(NA, NA, 19),
        linewidth = c(0.5, 0.8, NA)
    ))) +
    theme(legend.position = "bottom")

ggplot2::ggsave(snakemake@output[[1]], fig, width = 4.5, height = 5)
