import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import os
import jinja2
import treeflow_pipeline.results
from treeflow_pipeline.util import text_input
import treeflow_pipeline.model

DPI = 300
MAX_WIDTH_PIXELS = 2250
MAX_WIDTH = MAX_WIDTH_PIXELS/DPI
MAX_HEIGHT_PIXELS = 2625
MAX_HEIGHT = MAX_HEIGHT_PIXELS / DPI

param_name_mapping = {
    "pop_size": "Pop size",
    "clock_rate": "Clock rate",
    "rate_sd": "Rate prior scale",
    "rate_stats.mean": "Rate mean",
    "rate_stats.coefficientOfVariation": "Rate CoV",
    "tree.height": "Tree height",
    "tree.treeLength": "Tree length",
    "height": "Node heights",
    "rate": "Rates"
}

method_name_mapping = {
    "beast": "MCMC",
    "variational-samples-mean_field": "Variational (mean field)",
    "variational-samples-scaled": "Variational (scaled)"
}

np_func = lambda f: (lambda x: f(x).numpy())

def get_stats(log, sim_trace, name):
    run_filter = log[f"{name}.ESS"] > treeflow_pipeline.results.NUMERICAL_ISSUE_N
    trace_filtered = log[run_filter]
    return [sim_trace[name][run_filter]] + [trace_filtered[f"{name}.{stat}"] for stat in ["mean", "95%HPDlo", "95%HPDup"]]

BLUE = "#0384fc"
RED = "#ff3d51"

def coverage_plot(log_file_dict, sim_trace_file, model, stats, output, prior_scale=False):
    MARKER_WIDTH = 3
    TICK_LABEL_SIZE = 6

    params = model.free_params()
    stats_params = list(params.keys()) + stats
    sim_trace = pd.read_table(sim_trace_file)
    nrows = len(stats_params)
    ncols = len(log_file_dict)
    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(MAX_WIDTH, MAX_HEIGHT),
        dpi=DPI,
        constrained_layout=True
    )

    for i, name in enumerate(stats_params):
        axs[i, 0].set_ylabel(param_name_mapping[name])

    for j, name in enumerate(log_file_dict.keys()):
        ax = axs[0, j]
        ax.set_xlabel(method_name_mapping[name])
        ax.xaxis.set_label_position('top') 

    for j, (method, method_log_file) in enumerate(log_file_dict.items()):
        log = pd.read_table(method_log_file)
        for i, name in enumerate(stats_params):

            ax = axs[i, j]            
            if name in params and prior_scale:
                prior = treeflow_pipeline.model.get_dist(params[name])
                scale_functions = (np_func(prior.cdf), np_func(prior.quantile))
                ax.set_xscale("function", functions=scale_functions)
                ax.set_yscale("function", functions=scale_functions)
                
            true, mean, lower, upper = get_stats(log, sim_trace, name)
            covered = (lower <= true) & (true <= upper)

            ax.set_xlim(0, max(true))
            ax.set_ylim(0, max(upper[covered]))

            ax.vlines(true, lower, upper, color=np.where(covered, BLUE, RED), linewidths=MARKER_WIDTH, alpha=0.5)
            ax.scatter(true, mean, marker="_", color="black", s=MARKER_WIDTH ** 2, linewidth=1.0)
            ax.plot([0, max(true)], [0, max(true)], color="black", linestyle="--", linewidth=1.0)
            ax.tick_params(axis='both', which='major', labelsize=TICK_LABEL_SIZE)
    
    fig.savefig(output)

def coverage_table(coverage_table_file, stats, output):
    df = pd.read_csv(coverage_table_file)
    methods = df.method
    t = df.drop("method", axis=1)[stats].T
    t.columns = methods
    t = t[list(method_name_mapping.keys())]
    renamed = t.rename(columns=method_name_mapping, index=param_name_mapping)
    renamed.columns.name = "Method"
    renamed.index.name = "Statistic"
    (renamed * 100).to_latex(output, float_format="%.0f%%")


latex_jinja_env = jinja2.Environment(
    block_start_string='\BLOCK{',
    block_end_string='}',
    variable_start_string='\VAR{',
    variable_end_string='}',
    comment_start_string='\#{',
    comment_end_string='}',
    line_statement_prefix='%%',
    line_comment_prefix='%#',
    trim_blocks=True,
    autoescape=False,
    loader=jinja2.FileSystemLoader(".")
)

def build_manuscript(template_file, content_template_file,
                     figures_dict, tables_dict,
                     submission=False):
    if submission:
        figures = lambda x, args: None
    else:
        figures = lambda x, args: f"\\includegraphics{args}{{{str(os.path.splitext(figures_dict[x])[0])}}}"
    tables = { key: text_input(filename) for key, filename in tables_dict.items() }
    template = latex_jinja_env.get_template(content_template_file)
    return template.render(template=template_file, tables=tables, figures=figures)