# Multi-run marginals figure, shared by the carnivores and flu (H3N2) datasets.
#
# Compares BEAST 2 MCMC against TreeFlow VI, each run independently 4 times
# (see the multi_run_variational_fit rule in workflow/data.smk) so a Monte
# Carlo error band can be shown for each method: a bootstrap band for BEAST 2's
# MCMC samples, and an inter-run band (the spread of the per-run kernel density
# estimates) for the VI approximation. The `vi_samples_*` inputs are each the
# pooled samples from all 4 runs of one approximation, annotated with a `run`
# column; the main-text figures pass only root_full_rank, but the script takes
# any subset of the approximation families.

library(magrittr)

pythonExecutable <- snakemake@params[["python_executable"]]
reticulate::use_python(pythonExecutable)
pythonModule <- reticulate::import("treeflow_pipeline.manuscript")

# Named list mapping a display label to the VI approximation's snakemake@input
# key, e.g. list(`TreeFlow VI (full rank)` = "vi_samples_full_rank", ...).
viMethodInputs <- list(
    `TreeFlow VI` = "vi_samples_root_full_rank",
    `TreeFlow VI (full rank)` = "vi_samples_full_rank",
    `TreeFlow VI (mean field)` = "vi_samples_mean_field"
)

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

# Monte Carlo error band for an MCMC estimate: bootstrap the samples and take
# quantiles of the resulting kernel density estimates.
bootstrapDensityBands <- function(x, n_boot = 200, n_grid = 512, ci = 0.95) {
    rng <- range(x)
    padding <- 0.1 * diff(rng)
    grid <- seq(rng[1] - padding, rng[2] + padding, length.out = n_grid)
    bw <- stats::bw.nrd0(x)

    boot_densities <- matrix(NA, nrow = n_boot, ncol = n_grid)
    n <- length(x)
    for (i in seq_len(n_boot)) {
        boot_sample <- sample(x, n, replace = TRUE)
        d <- stats::density(boot_sample, bw = bw, from = grid[1], to = grid[n_grid], n = n_grid)
        boot_densities[i, ] <- d$y
    }

    alpha <- 1 - ci
    data.frame(
        Value = grid,
        ymin = apply(boot_densities, 2, quantile, probs = alpha / 2),
        ymax = apply(boot_densities, 2, quantile, probs = 1 - alpha / 2)
    )
}

# Monte Carlo error band for a VI estimate: the range of the per-run kernel
# density estimates across the independent runs.
interRunDensityBands <- function(x, run, n_grid = 512) {
    rng <- range(x)
    padding <- 0.1 * diff(rng)
    grid <- seq(rng[1] - padding, rng[2] + padding, length.out = n_grid)
    bw <- stats::bw.nrd0(x)

    run_densities <- sapply(split(x, run), function(xr) {
        stats::density(xr, bw = bw, from = grid[1], to = grid[n_grid], n = n_grid)$y
    })
    data.frame(
        Value = grid,
        ymin = apply(run_densities, 1, min),
        ymax = apply(run_densities, 1, max)
    )
}

# Load each available VI approximation's pooled samples (all four runs,
# annotated with `run`); a dataset may not have every approximation available.
viRawByMethod <- Filter(Negate(is.null), lapply(viMethodInputs, function(inputKey) {
    path <- snakemake@input[[inputKey]]
    if (is.null(path)) {
        return(NULL)
    }
    raw <- readr::read_csv(path)
    stopifnot("run" %in% colnames(raw))
    raw
}))
stopifnot(length(viRawByMethod) > 0)

viColumns <- setdiff(colnames(viRawByMethod[[1]]), "run")

dfs <- c(
    list(`BEAST 2` = readBeastTrace(snakemake@input[["beast_samples"]], viColumns)),
    lapply(viRawByMethod, function(raw) dplyr::select(raw, tidyselect::all_of(viColumns)))
)

stacked <- dplyr::bind_rows(dfs, .id = "Method")
renamed <- pythonModule$rename_marginal_df(stacked)
pivoted <- tidyr::pivot_longer(renamed, !Method, names_to = "variable", values_to = "Value")

variables <- setdiff(colnames(renamed), "Method")

# BEAST 2 Monte Carlo error band (bootstrap of the MCMC samples)
beastData <- dplyr::filter(renamed, Method == "BEAST 2")
beastRibbon <- do.call(rbind, lapply(variables, function(v) {
    bands <- bootstrapDensityBands(beastData[[v]])
    bands$variable <- v
    bands$Method <- "BEAST 2"
    bands
}))

# TreeFlow VI Monte Carlo error bands (variability between independent runs),
# one per approximation.
viRibbons <- do.call(rbind, lapply(names(viRawByMethod), function(methodLabel) {
    viRenamed <- pythonModule$rename_marginal_df(viRawByMethod[[methodLabel]])
    do.call(rbind, lapply(variables, function(v) {
        bands <- interRunDensityBands(viRenamed[[v]], viRenamed$run)
        bands$variable <- v
        bands$Method <- methodLabel
        bands
    }))
}))

ribbonDf <- rbind(beastRibbon, viRibbons)

postFig <- ggplot2::ggplot(pivoted) +
    ggplot2::geom_ribbon(
        data = ribbonDf,
        ggplot2::aes(x = Value, ymin = ymin, ymax = ymax, fill = Method),
        alpha = 0.3
    ) +
    ggplot2::geom_density(ggplot2::aes(Value, colour = Method)) +
    ggplot2::scale_y_continuous(name = "Density") +
    ggplot2::scale_x_continuous(n.breaks = 4) +
    ggplot2::facet_wrap(~variable, scales = "free")

ggplot2::ggsave(snakemake@output[[1]], postFig, width = 9, height = 7)
