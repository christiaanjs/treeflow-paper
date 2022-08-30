library(dplyr)
library(ggplot2)

plotDataFilename <- snakemake@input[[1]]
outFilename <- snakemake@output[[1]]

df <- readr::read_csv(plotDataFilename)
onlyLast <- group_by(df, across(c(-var_index, -value))) %>%
    arrange(desc(var_index)) %>%
    filter(row_number() == 1)
withPlotColumns <- mutate(onlyLast,
    var = paste(model_var_name, var_type),
    loss = var_type == "loss"
)
fig <- ggplot(withPlotColumns) +
    geom_line(aes(
        x = index,
        y = value,
        color = model_var_name,
        linetype = var_type
    )) +
    facet_grid(rows = vars(loss), cols = vars(method), scales = "free")
ggsave(outFilename, fig)
