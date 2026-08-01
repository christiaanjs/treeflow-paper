# Supplementary figure: how the choice of variational approximation family
# affects the fitted posterior marginals.
#
# One run of each of the three families TreeFlow provides for fixed-topology
# models -- mean_field (every coordinate independent), full_rank (a full
# covariance over every coordinate, including all node heights), and
# root_full_rank (a full covariance over the model parameters and the tree's
# root height, with the remaining node-height ratios independent) -- against
# the same BEAST 2 reference posterior. Single runs, since this figure is about
# the systematic differences between the families rather than their Monte Carlo
# error; the main-text figures show the inter-run band for root_full_rank.

library(magrittr)

pythonExecutable <- snakemake@params[["python_executable"]]
reticulate::use_python(pythonExecutable)
pythonModule <- reticulate::import("treeflow_pipeline.manuscript")

approxLabels <- c(
    mean_field = "TreeFlow VI (mean field)",
    full_rank = "TreeFlow VI (full rank)",
    root_full_rank = "TreeFlow VI (root full rank)"
)

approxes <- unlist(snakemake@params[["approxes"]])
samplePaths <- unlist(snakemake@input[["vi_samples"]])
stopifnot(length(approxes) == length(samplePaths))

readBeastTrace <- function(filename, columns, burnIn = 0.1) {
    raw <- readr::read_tsv(filename, comment = "#")
    burnedIn <- dplyr::filter(raw, dplyr::row_number() > nrow(raw) * burnIn)
    renamed <- dplyr::rename(burnedIn, tree_height = `tree.height`, tree_length = `tree.treeLength`) %>%
        dplyr::rename_with(
            ~ paste("frequencies", as.numeric(stringr::str_sub(.x, start = -1)) - 1, sep = "_"),
            tidyselect::starts_with("frequencies")
        )
    dplyr::select(renamed, tidyselect::all_of(columns))
}

viByMethod <- stats::setNames(
    lapply(samplePaths, readr::read_csv),
    approxLabels[approxes]
)

viColumns <- colnames(viByMethod[[1]])

dfs <- c(
    list(`Beast 2` = readBeastTrace(snakemake@input[["beast_samples"]], viColumns)),
    viByMethod
)

stacked <- dplyr::bind_rows(dfs, .id = "Method")
renamed <- pythonModule$rename_marginal_df(stacked)
pivoted <- tidyr::pivot_longer(renamed, !Method, names_to = "variable", values_to = "Value")

# BEAST 2 is the reference the approximations are being judged against, so give
# it a visually distinct (dashed, black) line rather than another colour.
methodLevels <- c("Beast 2", unname(approxLabels[approxes]))
pivoted$Method <- factor(pivoted$Method, levels = methodLevels)

fig <- ggplot2::ggplot(pivoted) +
    ggplot2::geom_density(
        ggplot2::aes(Value, colour = Method, linetype = Method)
    ) +
    ggplot2::scale_linetype_manual(
        values = c("dashed", rep("solid", length(approxes)))
    ) +
    ggplot2::scale_colour_manual(
        values = c("black", scales::hue_pal()(length(approxes)))
    ) +
    ggplot2::scale_y_continuous(name = "Density") +
    ggplot2::scale_x_continuous(n.breaks = 4) +
    ggplot2::facet_wrap(~variable, scales = "free") +
    ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(snakemake@output[[1]], fig, width = 9, height = 7.5)
