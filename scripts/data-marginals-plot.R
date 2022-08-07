print(snakemake)
library(magrittr)

readBeastTrace <- function(filename, columns, burnIn = 0.1) {
  raw <- readr::read_tsv(filename, comment = "#")
  burnedIn <- dplyr::filter(raw, dplyr::row_number() > nrow(raw) * burnIn)
  # dplyr::transmute(burnedIn, Kappa=kappa, `Birth rate`=birthRate, `Gamma shape`=gammaShape, `Tree height`=TreeHeight, )
  renamed <- dplyr::rename(burnedIn, tree_height = `tree.height`)
  dplyr::select(renamed, tidyselect::all_of(columns))
}

readViTrace <- function(filename) {
  raw <- readr::read_csv(filename)
  # dplyr::transmute(raw, Kappa=kappa, `Birth rate`=birth_rate, `Gamma shape`=site_gamma_shape, `Tree height`=tree_height)
}

viTrace <- readViTrace(snakemake@input["vi_samples"])
dfs <- list(
  `Beast 2` = readBeastTrace(snakemake@input["beast_samples"], colnames(viTrace)),
  `Treeflow VI` = viTrace
)

stacked <- dplyr::bind_rows(dfs, .id = "Method")
pivoted <- tidyr::pivot_longer(stacked, !Method, names_to = "variable")

alpha <- 0.05
fl <- purrr::partial(signif, digits = 3)

summaryTable <- dplyr::group_by(pivoted, Method, variable) %>%
  dplyr::summarise(
    mean = mean(value),
    lower = quantile(value, alpha / 2.0),
    upper = quantile(value, 1.0 - alpha / 2.0)
  ) %>%
  dplyr::mutate(summaryString = glue::glue("{fl(mean)} ({fl(lower)}, {fl(upper)})")) %>%
  dplyr::select(Method, variable, summaryString) %>%
  tidyr::pivot_wider(names_from = variable, values_from = summaryString)

print(summaryTable)

tableGrob <- gridExtra::tableGrob(summaryTable)

fig <- ggplot2::ggplot(pivoted) +
  ggplot2::geom_density(ggplot2::aes(value, colour = Method)) +
  ggplot2::facet_wrap(~variable, scales = "free")

composite <- gridExtra::grid.arrange(fig, tableGrob, heights = c(3, 2), nrow = 2)

ggplot2::ggsave(snakemake@output[[1]], composite, width = 8, height = 6)
