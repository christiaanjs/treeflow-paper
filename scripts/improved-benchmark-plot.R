# plotDataPath <- "/home/christiaan/uni/treeflow-benchmarks/out/plot-data.csv"
# outFile <- "manuscript/figures/benchmark-log-scale-plot.png"
# pythonExecutable <- "/home/christiaan/.pyenv/versions/3.8.13/envs/treeflow/bin/python3.8"


plotDataPath <- snakemake@input[["plot_data"]]
outFile <- snakemake@output[["plot"]]
pythonExecutable <- snakemake@params[["python_executable"]]
reticulate::use_python(pythonExecutable)
pythonModule <- reticulate::import("treeflow_pipeline.manuscript")
renameFunc <- pythonModule$get_benchmark_colname
df <- readr::read_csv(plotDataPath, show_col_types = FALSE)
remappedDf <- pythonModule$remap_and_sort_benchmark_df(df)
ggsaveArgs <- list(width = 8, height = 3)
treeflowbenchmarksr::comparisonPlot(
    remappedDf,
    outFile = outFile,
    ggsaveArgs = ggsaveArgs,
    renameFunc = renameFunc,
    colorFunc = purrr::partial(paste, sep = ", ")
)
