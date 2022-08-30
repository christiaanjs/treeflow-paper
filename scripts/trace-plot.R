# plotDataFilename <- snakemake@input["vi_samples"]

library(dplyr)
library(ggplot2)

plotDataFilename <- "out/dengue_coal_easy/trace-plot-data.csv"
df <- readr::read_csv(plotDataFilename)
onlyLast <- group_by(df, across(c(-var_index))) %>%
    arrange(desc(var_index)) %>%
    filter(row_number() == 1)
withPlotColumns <- mutate(onlyLast,
    var = paste(model_var_name, var_type),
    loss = var_type == "loss"
)
fig <- ggplot(withPlotColumns) +
    geom_line(aes(x = index, y = value, color = var)) +
    facet_grid(rows = vars(loss), cols = vars(method), scales = "free_y")
