# Regenerates the phylogenetic likelihood benchmark figure (Figure 5).
#
# Layout: one panel per benchmarked task (model x computation), with the three
# software implementations (TreeFlow, bito/BEAGLE, JAX) shown as separate series.
# This makes the difference in scaling (the slope on log-log axes) between the
# implementations directly comparable within each task, addressing the reviewer
# suggestion to split the benchmark results by task rather than by software.
#
# Previously this plot faceted by software with the four tasks as coloured lines
# (via treeflowbenchmarksr::comparisonPlot); it is now self-contained and reads
# plot-data.csv directly.

library(readr)
library(dplyr)
library(ggplot2)

plotDataPath <- snakemake@input[["plot_data"]]
outFile <- snakemake@output[["plot"]]

df <- readr::read_csv(plotDataPath, show_col_types = FALSE)

# Display labels matching treeflow_pipeline.manuscript. The inlined
# treeflow/experiments benchmark adds the native C++-op series (treeflow_native).
# For bito/BEAGLE only the direct benchmarkable (beagle_bito_direct, which drives
# the bito instance without the TensorFlow wrapper) is shown; the wrapper-overhead
# variant (beagle_bito) is left to the benchmark's own exploratory plots.
methodLabels <- c(
    treeflow = "TreeFlow",
    treeflow_native = "TreeFlow (native)",
    beagle_bito_direct = "bito/BEAGLE",
    jax = "JAX",
    jax_jit = "JAX (JIT)"
)
modelLabels <- c(jc = "JC", full = "GTR/Weibull")
computationLabels <- c(
    likelihood_time = "Likelihood",
    phylo_gradients_time = "Gradients"
)

methodOrdering <- c("TreeFlow", "TreeFlow (native)", "bito/BEAGLE", "JAX", "JAX (JIT)")
# Panels grouped by model (rows) and computation (columns) when nrow = 2.
taskOrdering <- c(
    "Likelihood, JC",
    "Gradients, JC",
    "Likelihood, GTR/Weibull",
    "Gradients, GTR/Weibull"
)

plotDf <- df %>%
    mutate(
        Method = factor(methodLabels[method], levels = methodOrdering),
        Model = modelLabels[model],
        Computation = computationLabels[computation],
        Task = factor(paste(Computation, Model, sep = ", "), levels = taskOrdering)
    ) %>%
    # Drop any method not in methodLabels (e.g. jax_jit, which the manuscript
    # export excludes) so it never becomes a stray NA series.
    filter(!is.na(Method)) %>%
    group_by(Task, Method, taxon_count) %>%
    summarise(time = mean(time), .groups = "drop")

plot <- ggplot(plotDf, aes(x = taxon_count, y = time, colour = Method)) +
    geom_line() +
    geom_point(size = 1) +
    facet_wrap(~Task, nrow = 2) +
    scale_x_continuous(trans = "log2") +
    scale_y_continuous(trans = "log10") +
    labs(x = "Taxon count", y = "Time (s)", colour = "Software")

ggplot2::ggsave(filename = outFile, plot = plot, width = 8, height = 6)
