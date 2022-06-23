setClass("Snakemake", representation(input="list"))
snakemake <- new("Snakemake", input=list(
  beast_trace="../out/carnivores-beast2.log",
  vi_trace="../../treeflow/examples/demo-out/carnivores-base-samples.csv"
))

readBeastTrace <- function(filename, burnIn=0.1){
  raw <- readr::read_tsv(filename, comment="#")
  burnedIn <- dplyr::filter(raw, dplyr::row_number() > nrow(raw) * burnIn)
  dplyr::transmute(burnedIn, Kappa=kappa, `Tree height`=TreeHeight, `Birth rate`=birthRate, `Gamma shape`=gammaShape)
}

dfs <- list(
  beast=readBeastTrace(snakemake@input["beast_trace"]),
  vi=readr::read_csv(snakemake@input["vi_trace"])
)
