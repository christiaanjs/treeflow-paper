import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import treeflow_pipeline.model
import treeflow_pipeline.results

DPI = 300
FULL_WIDTH = 2250
MAX_WIDTH = 2250/DPI

param_name_mapping = {
    "pop_size": "Population size",
    "clock_rate": "Clock rate",
    "rate_sd": "Branch rate prior scale",
    "rate_stats.mean": "Branch rate mean",
    "rate_stats.coefficientOfVariation": "Branch rate coeff of var",
    "tree.height": "Tree height",
    "tree.treeLength": "Tree length"
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

    params = model.free_params()
    stats_params = list(params.keys()) + stats
    sim_trace = pd.read_table(sim_trace_file)
    nrows = len(stats_params)
    ncols = len(log_file_dict)
    width = MAX_WIDTH
    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(width, (width/ncols)*nrows),
        dpi=DPI,
        constrained_layout=True
    )

    for i, name in enumerate(stats_params):
        axs[i, 0].set_ylabel(param_name_mapping[name])

    for j, name in enumerate(log_file_dict.keys()):
        axs[-1, j].set_xlabel(method_name_mapping[name])

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
            ax.tick_params(axis='both', which='major', labelsize=MARKER_WIDTH)
    
    fig.savefig(output)

