print(snakemake)
library(magrittr)

readBeastTrace <- function(filename, columns, burnIn = 0.1) {
  raw <- readr::read_tsv(filename, comment = "#")
  burnedIn <- dplyr::filter(raw, dplyr::row_number() > nrow(raw) * burnIn)
  # dplyr::transmute(burnedIn, Kappa=kappa, `Birth rate`=birthRate, `Gamma shape`=gammaShape, `Tree height`=TreeHeight, )
  renamed <- dplyr::rename(burnedIn, tree_height = `tree.height`, tree_length = `tree.treeLength`) %>%
    dplyr::rename_with(
      ~ paste("frequencies", as.numeric(stringr::str_sub(.x, start = -1)) - 1, sep = "_"),
      tidyselect::starts_with("frequencies")
    )
  dplyr::select(renamed, tidyselect::all_of(columns))
}

readViTrace <- function(filename) {
  raw <- readr::read_csv(filename)
  # dplyr::transmute(raw, Kappa=kappa, `Birth rate`=birth_rate, `Gamma shape`=site_gamma_shape, `Tree height`=tree_height)
  raw
}

viTrace <- readViTrace(snakemake@input["vi_samples"])

includeMl <- snakemake@config[["include_ml_plot"]]
includeEmp <- snakemake@config[["include_empirical_plot"]]
includePointEsts <- includeMl || includeEmp

mlPivoted <- if (includeMl) {
  mlTrace <- readViTrace(snakemake@input["ml_variables"])
  dplyr::mutate(mlTrace, Method = "Treeflow ML") %>%
    tidyr::pivot_longer(!Method, names_to = "variable")
} else {
  NULL
}


dfs <- list(
  `Beast 2` = readBeastTrace(snakemake@input["beast_samples"], colnames(viTrace)),
  `Treeflow VI` = viTrace
)

stacked <- dplyr::bind_rows(dfs, .id = "Method")
pivoted <- tidyr::pivot_longer(stacked, !Method, names_to = "variable")
pointPivoted <- if (includeEmp && ("frequencies_0" %in% colnames(viTrace))) {
  freqsDf <- readr::read_csv(snakemake@input["empirical_frequencies"])
  freqsPivoted <- dplyr::mutate(freqsDf, Method = "Empirical") %>%
    tidyr::pivot_longer(!Method, names_to = "variable")
  dplyr::bind_rows(mlPivoted, freqsPivoted)
} else {
  mlPivoted
}

alpha <- 0.05
fl <- purrr::partial(signif, digits = 3)

postSummaryTableNarrow <- dplyr::group_by(pivoted, Method, variable) %>%
  dplyr::summarise(
    mean = mean(value),
    lower = quantile(value, alpha / 2.0),
    upper = quantile(value, 1.0 - alpha / 2.0)
  ) %>%
  dplyr::mutate(summaryString = glue::glue("{fl(mean)} ({fl(lower)}, {fl(upper)})")) %>%
  dplyr::select(Method, variable, summaryString)

summaryTableNarrow <- if (includePointEsts) {
  pointSummaryTableNarrow <- dplyr::transmute(pointPivoted, Method, variable, summaryString = as.character(fl(value)))
  dplyr::bind_rows(postSummaryTableNarrow, pointSummaryTableNarrow)
} else {
  postSummaryTableNarrow
}

summaryTable <- tidyr::pivot_wider(summaryTableNarrow, names_from = Method, values_from = summaryString)

print(summaryTable)

tableGrob <- gridExtra::tableGrob(summaryTable)
# pairPlotData <- dplyr::select(stacked, !tidyselect::starts_with("frequencies"))
# pairPlot <- GGally::ggpairs(pairPlotData, ggplot2::aes(color = Method)) %>% GGally::ggmatrix_gtable()
postFig <- ggplot2::ggplot(pivoted) +
  ggplot2::geom_density(ggplot2::aes(value, colour = Method))
unfacetedFig <- if (includePointEsts) {
  postFig + ggplot2::geom_vline(ggplot2::aes(xintercept = value, colour = Method), data = pointPivoted)
} else {
  postFig
}
fig <- unfacetedFig + ggplot2::facet_wrap(~variable, scales = "free")

# composite <- gridExtra::grid.arrange(fig, pairPlot, tableGrob, heights = c(3, 5, 3), nrow = 3)
composite <- gridExtra::grid.arrange(fig, tableGrob, heights = c(3, 2), nrow = 2)

ggplot2::ggsave(snakemake@output[[1]], composite, width = 8, height = 8)
