library(magrittr)

pythonExecutable <- snakemake@params[["python_executable"]]
reticulate::use_python(pythonExecutable)
pythonModule <- reticulate::import("treeflow_pipeline.manuscript")

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

readViTrace <- function(filename) {
    raw <- readr::read_csv(filename)
    raw
}

# Compute bootstrap density bands for a numeric vector
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

viTrace <- readViTrace(snakemake@input["vi_samples"])

dfs <- list(
    `Beast 2` = readBeastTrace(snakemake@input["beast_samples"], colnames(viTrace)),
    `Treeflow VI` = viTrace
)

stacked <- dplyr::bind_rows(dfs, .id = "Method")
renamed <- pythonModule$rename_marginal_df(stacked)
pivoted <- tidyr::pivot_longer(renamed, !Method, names_to = "variable", values_to = "Value")

# Compute bootstrap bands for Beast 2 densities
beastData <- dplyr::filter(renamed, Method == "Beast 2")
variables <- setdiff(colnames(renamed), "Method")
ribbonDf <- do.call(rbind, lapply(variables, function(v) {
    bands <- bootstrapDensityBands(beastData[[v]])
    bands$variable <- v
    bands
}))

postFig <- ggplot2::ggplot(pivoted) +
    ggplot2::geom_ribbon(
        data = ribbonDf,
        ggplot2::aes(x = Value, ymin = ymin, ymax = ymax),
        fill = "grey70", alpha = 0.4
    ) +
    ggplot2::geom_density(ggplot2::aes(Value, colour = Method)) +
    ggplot2::scale_y_continuous(name = "Density") +
    ggplot2::scale_x_continuous(n.breaks = 4) +
    ggplot2::facet_wrap(~variable, scales = "free")

ggplot2::ggsave(snakemake@output[[1]], postFig, width = 8, height = 6)
