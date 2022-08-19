print(snakemake)
library(magrittr)

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
  # dplyr::transmute(raw, Kappa=kappa, `Birth rate`=birth_rate, `Gamma shape`=site_gamma_shape, `Tree height`=tree_height)
  raw
}

viTrace <- readViTrace(snakemake@input["vi_samples"])
mlTrace <- readViTrace(snakemake@input["ml_variables"])
dfs <- list(
  `Beast 2` = readBeastTrace(snakemake@input["beast_samples"], colnames(viTrace)),
  `Treeflow VI` = viTrace
)
stacked <- dplyr::bind_rows(dfs, .id = "Method")
pivoted <- tidyr::pivot_longer(stacked, !Method, names_to = "variable")
mlPivoted <- dplyr::mutate(mlTrace, Method = "Treeflow ML") %>%
  tidyr::pivot_longer(!Method, names_to = "variable")

alpha <- 0.05
fl <- purrr::partial(signif, digits = 3)

mlSummaryTable <- dplyr::transmute(mlPivoted, Method, variable, summaryString = as.character(value))
summaryTable <- dplyr::group_by(pivoted, Method, variable) %>%
  dplyr::summarise(
    mean = mean(value),
    lower = quantile(value, alpha / 2.0),
    upper = quantile(value, 1.0 - alpha / 2.0)
  ) %>%
  dplyr::mutate(summaryString = glue::glue("{fl(mean)} ({fl(lower)}, {fl(upper)})")) %>%
  dplyr::select(Method, variable, summaryString) %>%
  dplyr::bind_rows(mlSummaryTable) %>%
  tidyr::pivot_wider(names_from = Method, values_from = summaryString)

print(summaryTable)

tableGrob <- gridExtra::tableGrob(summaryTable)
# pairPlot <- GGally::ggpairs(stacked, ggplot2::aes(color = Method)) %>% GGally::ggmatrix_gtable()
fig <- ggplot2::ggplot(pivoted) +
  ggplot2::geom_density(ggplot2::aes(value, colour = Method)) +
  ggplot2::geom_vline(ggplot2::aes(xintercept = value, colour = Method), data = mlPivoted) +
  ggplot2::facet_wrap(~variable, scales = "free")

# composite <- gridExtra::grid.arrange(fig, pairPlot, tableGrob, heights = c(3, 5, 2), nrow = 3)
composite <- gridExtra::grid.arrange(fig, tableGrob, heights = c(3, 2), nrow = 2)

ggplot2::ggsave(snakemake@output[[1]], composite, width = 8, height = 8)
