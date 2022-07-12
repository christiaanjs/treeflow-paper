readBeastTrace <- function(filename, burnIn=0.1){
  raw <- readr::read_tsv(filename, comment="#")
  burnedIn <- dplyr::filter(raw, dplyr::row_number() > nrow(raw) * burnIn)
  dplyr::transmute(burnedIn, Kappa=kappa, `Birth rate`=birthRate, `Gamma shape`=gammaShape, `Tree height`=TreeHeight, )
}

readViTrace <- function(filename){
  raw <- readr::read_csv(filename)
  dplyr::transmute(raw, Kappa=kappa, `Birth rate`=birth_rate, `Gamma shape`=site_gamma_shape, `Tree height`=tree_height)
}

dfs <- list(
  `Beast 2`=readBeastTrace(snakemake@params["beast_trace"]),
  `Treeflow VI`=readViTrace(snakemake@params["vi_trace"])
)

stacked <- dplyr::bind_rows(dfs, .id="Method")
pivoted <- tidyr::pivot_longer(stacked, !Method, names_to="variable")

fig <- ggplot2::ggplot(pivoted) +
  ggplot2::geom_density(ggplot2::aes(value, colour=Method)) +
  ggplot2::facet_wrap(~variable, scales="free")

ggplot2::ggsave(snakemake@output[[1]], fig, width=8, height=6)
