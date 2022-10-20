from functools import reduce
import pathlib
import typing as tp
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("AGG")
import matplotlib.pyplot as plt
import os
import jinja2
from treeflow import parse_newick
import treeflow_pipeline.results
from treeflow_pipeline.util import text_input
import treeflow_pipeline.model

DPI = 300
MAX_WIDTH_PIXELS = 2250
MAX_WIDTH = MAX_WIDTH_PIXELS / DPI
MAX_HEIGHT_PIXELS = 2625
MAX_HEIGHT = MAX_HEIGHT_PIXELS / DPI

param_name_mapping = {
    "pop_size": "Pop size",
    "clock_rate": "Clock rate",
    "rate_sd": "Rate prior scale",
    "rate_stats.mean": "Rate mean",
    "rate_stats.coefficientOfVariation": "Rate CoV",
    "tree_height": "Tree height",
    "tree_length": "Tree length",
    "height": "Node heights",
    "rate": "Rates",
    "site_gamma_shape": "Site rate shape",
    "birth_rate": "Birth rate",
    **{f"frequencies_{i}": f"Frequencies ({char})" for i, char in enumerate("ACGT")},
}

method_name_mapping = {
    "beast": "MCMC",
    "variational-samples-mean_field": "Variational (mean field)",
    "variational-samples-scaled": "Variational (scaled)",
}

np_func = lambda f: (lambda x: f(x).numpy())


def get_stats(log, sim_trace, name):
    run_filter = log[f"{name}.ESS"] > treeflow_pipeline.results.NUMERICAL_ISSUE_N
    trace_filtered = log[run_filter]
    return [sim_trace[name][run_filter]] + [
        trace_filtered[f"{name}.{stat}"] for stat in ["mean", "95%HPDlo", "95%HPDup"]
    ]


BLUE = "#0384fc"
RED = "#ff3d51"


def coverage_plot(
    log_file_dict, sim_trace_file, model, stats, output, prior_scale=False
):
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
        constrained_layout=True,
    )

    for i, name in enumerate(stats_params):
        axs[i, 0].set_ylabel(param_name_mapping[name])

    for j, name in enumerate(log_file_dict.keys()):
        ax = axs[0, j]
        ax.set_xlabel(method_name_mapping[name])
        ax.xaxis.set_label_position("top")

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

            ax.vlines(
                true,
                lower,
                upper,
                color=np.where(covered, BLUE, RED),
                linewidths=MARKER_WIDTH,
                alpha=0.5,
            )
            ax.scatter(
                true,
                mean,
                marker="_",
                color="black",
                s=MARKER_WIDTH**2,
                linewidth=1.0,
            )
            ax.plot(
                [0, max(true)],
                [0, max(true)],
                color="black",
                linestyle="--",
                linewidth=1.0,
            )
            ax.tick_params(axis="both", which="major", labelsize=TICK_LABEL_SIZE)

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


colname_mapping = {
    "method": "Method",
    "computation": "Computation",
    "model": "Model",
    "slope": "log-log slope",
}
inverse_colname_mapping = {value: key for key, value in colname_mapping.items()}

model_mapping = {"jc": "JC", "full": "GTR/Weibull"}
method_mapping = {
    "treeflow": "TreeFlow",
    "jax": "JAX",
    "beagle_bito_direct": "bito/BEAGLE",
}
computation_mapping = {
    "likelihood_time": "Likelihood",
    "phylo_gradients_time": "Gradients",
}
value_mappings = dict(
    method=method_mapping,
    model=model_mapping,
    computation=computation_mapping,
)
index_columns = ["model", "computation", "method"]

computation_ordering = ["Likelihood", "Gradients"]
model_ordering = ["JC", "GTR/Weibull"]
method_ordering = ["TreeFlow", "bito/BEAGLE", "JAX"]

orderings_dict = dict(
    Method=method_ordering, Computation=computation_ordering, Model=model_ordering
)
sort_key_dict = {
    colname: {value: i for i, value in enumerate(ordering)}
    for colname, ordering in orderings_dict.items()
}

benchmark_plot_colname_mapping = dict(
    colname_mapping, colour="Task", taxon_count="Taxon count", time="Time (s)"
)


def build_col_mapping_df(mapping_dict, ordering, colname, reverse_order=False):
    return pd.DataFrame(
        list(mapping_dict.items()), columns=[colname, colname_mapping[colname]]
    ).merge(
        pd.DataFrame(
            list(enumerate(ordering[::-1] if reverse_order else ordering)),
            columns=[f"{colname}_index", colname_mapping[colname]],
        )
    )


def get_benchmark_colname(colname: tp.Union[tp.List[str], str]):
    if isinstance(colname, str):
        return benchmark_plot_colname_mapping.get(colname, colname)
    else:
        return [get_benchmark_colname(x) for x in colname]


def remap_and_sort_benchmark_df(
    df: pd.DataFrame, reverse_order=False
):  # columns method, computation, model, taxon_count, seed
    mapping_dfs = [
        build_col_mapping_df(
            value_mappings[colname],
            orderings_dict[colname_mapping[colname]],
            colname,
            reverse_order=reverse_order,
        )
        for colname in index_columns
    ]
    merged = reduce(pd.DataFrame.merge, mapping_dfs, df)
    sorted = merged.sort_values([f"{colname}_index" for colname in index_columns])

    index_columns_set = set(index_columns)
    renamed = sorted[
        [
            (get_benchmark_colname(col) if (col in index_columns_set) else col)
            for col in df.columns
        ]
    ].rename(columns=inverse_colname_mapping)
    return renamed


def rename_marginal_df(df: pd.DataFrame):
    return df.rename(columns=param_name_mapping)


def benchmark_summary_table(
    plot_data_path, fit_table_path, output_path, taxon_counts=(512,)
):
    fit_table = pd.read_csv(fit_table_path)
    plot_data = pd.read_csv(plot_data_path)

    time_summaries = (
        plot_data[plot_data.taxon_count.isin(taxon_counts)]
        .groupby(index_columns + ["taxon_count"])["time"]
        .mean()
        .to_frame()
        .reset_index()
    )
    times_pivoted = time_summaries.pivot(index=index_columns, columns="taxon_count")

    times_pivoted_renamed = times_pivoted.set_axis(
        times_pivoted.columns.to_flat_index().map(
            lambda var_name: f"Mean time ({var_name[1]} taxa)"
        ),
        axis=1,
    ).reset_index()
    merged = times_pivoted_renamed.merge(fit_table[list(colname_mapping.keys())])
    improved = merged.replace(
        dict(
            method=method_mapping, model=model_mapping, computation=computation_mapping
        )
    ).rename(columns=colname_mapping)
    indexed = improved.set_index(
        [colname_mapping[x] for x in index_columns]
    ).sort_index(key=lambda x: x.to_series().replace(sort_key_dict[x.name]))

    model_index = index_columns.index("model")
    comp_index = index_columns.index("computation")

    def style_func(row):
        if (
            row.name[model_index] == model_mapping["full"]
            and row.name[comp_index] == computation_mapping["phylo_gradients_time"]
        ):
            return ["font-weight: bold;" for x in row]
        else:
            return ["" for x in row]

    styled = indexed.style.apply(style_func, axis=1).format(precision=2)
    styled.to_latex(output_path, convert_css=True, multirow_align="t", hrules=True)


latex_jinja_env = jinja2.Environment(
    block_start_string="\BLOCK{",
    block_end_string="}",
    variable_start_string="\VAR{",
    variable_end_string="}",
    comment_start_string="\#{",
    comment_end_string="}",
    line_statement_prefix="%%",
    line_comment_prefix="%#",
    trim_blocks=True,
    autoescape=False,
    loader=jinja2.FileSystemLoader("."),
)


def get_treeflow_manuscript_vars(treeflow_benchmarks_config, flu_dataset, out_dir):
    out_dir = pathlib.Path(out_dir)
    flu_dir = out_dir / flu_dataset
    flu_tree = parse_newick(str(flu_dir / "topology.nwk"))
    return dict(
        min_sequence_count=min(treeflow_benchmarks_config["full_taxon_counts"]),
        max_sequence_count=max(treeflow_benchmarks_config["full_taxon_counts"]),
        sequence_length=treeflow_benchmarks_config["sequence_length"],
        replicate_count=treeflow_benchmarks_config["replicates"],
        sample_count=treeflow_benchmarks_config["sample_count"],
        flu_yaml_file=str(flu_dir / "model.yaml"),
        flu_taxon_count=flu_tree.taxon_count,
    )


def build_manuscript(
    template_file,
    content_template_file,
    figures_dict,
    tables_dict,
    vars,
    submission=False,
):
    if submission:
        figures = lambda x, args: None
    else:
        figures = (
            lambda x, args: f"\\includegraphics{args}{{{str(os.path.splitext(figures_dict[x])[0])}}}"
        )
    tables = {key: text_input(filename) for key, filename in tables_dict.items()}
    template = latex_jinja_env.get_template(content_template_file)
    return template.render(
        template=template_file,
        tables=tables,
        figures=figures,
        submission=submission,
        **vars,
    )
