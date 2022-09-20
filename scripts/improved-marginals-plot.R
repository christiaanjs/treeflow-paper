library(magrittr)

pythonExecutable <- snakemake@params[["python_executable"]]
reticulate::use_python(pythonExecutable)
pythonModule <- reticulate::import("treeflow_pipeline.manuscript")

readBeastTrace <- function(filename, columns, burnIn = 0.1) {
    raw <- readr::read_tsv(filename, comment = "#")
    burnedIn <- dplyr::filter(raw, dplyr::row_number() > nrow(raw) * burnIn)
    # dplyr::transmute(burnedIn, Kappa=kappa, `Birth rate`=birthRate, `Gamma shape`=gammaShape, `Tree height`=TreeHeight, )
    renamed <- dplyr::rename(burnedIn, tree_height = `tree.height`) %>%
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

viTrace <- readViTrace(snakemake@input["vi_samples"])

dfs <- list(
    `Beast 2` = readBeastTrace(snakemake@input["beast_samples"], colnames(viTrace)),
    `Treeflow VI` = viTrace
)

stacked <- dplyr::bind_rows(dfs, .id = "Method")
renamed <- pythonModule$rename_marginal_df(stacked)
pivoted <- tidyr::pivot_longer(renamed, !Method, names_to = "variable", values_to = "Value")
postFig <- ggplot2::ggplot(pivoted) +
    ggplot2::geom_density(ggplot2::aes(Value, colour = Method)) +
    ggplot2::scale_y_continuous(name = "Density") +
    ggplot2::facet_wrap(~variable, scales = "free")

ggplot2::ggsave(snakemake@output[[1]], postFig, width = 8, height = 6)
