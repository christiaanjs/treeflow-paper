import typing as tp
import dendropy
import numpy as np
from treeflow.tree.io import parse_newick, write_tensor_trees
from treeflow.evolution.seqio import parse_fasta
from treeflow.model.ml import MLResults
from treeflow.vi.util import VIResults
from treeflow.model.io import flatten_samples_to_dict, calculate_tree_stats
import treeflow_pipeline.model
import pandas as pd
import io
import arviz
from treeflow_pipeline.util import pickle_input, beast_log_input


def construct_precedes_map(
    taxon_order,
):  # [x, y] = True if x precedes y TODO: Do something that isn't quadratic
    precedes = {}
    for i in range(len(taxon_order)):
        for j in range(i):
            precedes[taxon_order[i], taxon_order[j]] = False
        for j in range(i + 1, len(taxon_order)):
            precedes[taxon_order[i], taxon_order[j]] = True
    return precedes


def construct_leaf_descendant_map(tree):  # [x] = y only if y is a leaf descendant of x
    leaf_descendants = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            leaf_descendants[node] = node
        else:
            leaf_descendants[node] = leaf_descendants[node.child_nodes()[0]]
    return leaf_descendants


def postorder_node_iter(tree, taxon_order):
    if taxon_order is None:
        return tree.postorder_node_iter()
    else:
        precedes = construct_precedes_map(taxon_order)
        leaf_desc_map = construct_leaf_descendant_map(tree)
    stack = [tree.seed_node]
    reverse = []
    while len(stack) > 0:
        reverse.append(stack.pop())
        if not reverse[-1].is_leaf():
            left, right = reverse[-1].child_nodes()
            if precedes[
                leaf_desc_map[left].taxon.label, leaf_desc_map[right].taxon.label
            ]:
                stack.append(left)
                stack.append(right)
            else:
                stack.append(right)
                stack.append(left)

    return reverse[::-1]


def parse_beast_trees(tree_file, format="nexus", metadata_keys=[], taxon_order=None):
    trees = dendropy.TreeList.get(
        path=tree_file,
        schema=format,
        rooting="default-rooted",
        preserve_underscores=True,
    )
    taxon_count = len(trees.taxon_namespace)
    branch_lengths = np.zeros((len(trees), 2 * taxon_count - 2))
    metadata = {
        key: np.zeros((len(trees), 2 * taxon_count - 2)) for key in metadata_keys
    }

    label_length = max([len(taxon.label) for taxon in trees.taxon_namespace])
    label_dtype = ">U{0}".format(label_length)
    taxon_names = np.empty((len(trees), taxon_count), dtype=label_dtype)

    for i, tree in enumerate(trees):
        leaf_index = 0
        node_index = taxon_count
        for node in list(postorder_node_iter(tree, taxon_order))[:-1]:
            if node.is_leaf():
                branch_lengths[i, leaf_index] = node.edge.length
                taxon_names[i, leaf_index] = node.taxon.label
                for key in metadata_keys:
                    metadata[key][i, leaf_index] = float(node.annotations[key].value)
                leaf_index += 1
            else:
                branch_lengths[i, node_index] = node.edge.length
                for key in metadata_keys:
                    metadata[key][i, node_index] = float(node.annotations[key].value)
                node_index += 1

    return dict(
        branch_lengths=branch_lengths, taxon_names=taxon_names, metadata=metadata
    )


def get_heights(branch_lengths, parent_indices, preorder_indices):
    root_distances = np.zeros(branch_lengths.shape[:-1] + (len(preorder_indices),))
    for i in preorder_indices[1:]:
        root_distances[..., i] = (
            root_distances[..., parent_indices[i]] + branch_lengths[..., i]
        )
    return np.max(root_distances, axis=-1, keepdims=True) - root_distances


def remove_burn_in(x, burn_in):
    return x[np.arange(len(x)) > int(len(x) * burn_in)]


def process_beast_results(tree_file, trace_file, topology_file, beast_config, model):
    burn_in = beast_config["burn_in"]
    trace = pd.read_table(trace_file, comment="#")

    trees = dendropy.TreeList.get(
        path=tree_file,
        schema="nexus",
        rooting="default-rooted",
        preserve_underscores=True,
    )
    tree, taxon_names = treeflow.tree_processing.parse_newick(topology_file)
    topology = treeflow.tree_processing.update_topology_dict(tree["topology"])

    relaxed = model.clock_model in treeflow_pipeline.model.RELAXED_CLOCK_MODELS
    trees_parsed = parse_beast_trees(
        tree_file, metadata_keys=(["rate"] if relaxed else []), taxon_order=taxon_names
    )

    branch_lengths = trees_parsed["branch_lengths"]
    heights = get_heights(
        branch_lengths, topology["parent_indices"], topology["preorder_indices"]
    )

    result = dict(
        {
            param: remove_burn_in(trace[param].values, burn_in)
            for param in model.free_params()
        },
        branch_lengths=remove_burn_in(branch_lengths, burn_in),
        heights=remove_burn_in(heights, burn_in),
    )

    if relaxed:
        absolute_rates = trees_parsed["metadata"]["rate"]
        if "clock_rate" in trace.columns:
            rates = absolute_rates / trace.clock_rate[:, np.newaxis]
        else:
            rates = absolute_rates
        result["absolute_rates"] = remove_burn_in(absolute_rates, burn_in)
        result["rates"] = remove_burn_in(rates, burn_in)

    return result


def get_variational_samples(
    variational_fit,
    topology_file,
    model,
    clock_approx,
    starting_values,
    trace_out,
    tree_out,
    seed,
    n=1000,
):
    if fit_successful(variational_fit):
        approx = treeflow_pipeline.model.reconstruct_approx(
            topology_file, variational_fit, model, clock_approx
        )

        samples = approx.sample(n, seed=seed)
    else:  # In case there is are numerical issues # TODO: Fix
        n = NUMERICAL_ISSUE_N
        approx = treeflow_pipeline.model.construct_approx(
            topology_file, model, clock_approx
        )
        samples = approx.sample(n, seed=seed)

    branch_lengths = treeflow.sequences.get_branch_lengths(samples["tree"]).numpy()

    result_dict = {key: samples[key].numpy() for key in model.free_params()}

    heights = samples["tree"]["heights"].numpy()
    trace_dict = {
        "tree.height": heights[:, -1],
        "tree.treeLength": np.sum(branch_lengths, axis=1),
        **result_dict,
    }

    result_dict["branch_lengths"] = branch_lengths
    result_dict["heights"] = heights
    if model.relaxed_clock():
        result_dict["rates"] = samples["rates"].numpy()
        if "clock_rate" in result_dict:
            absolute_rates = (
                result_dict["clock_rate"][:, np.newaxis] * result_dict["rates"]
            )
        else:
            absolute_rates = result_dict["rates"]
        result_dict["absolute_rates"] = absolute_rates

        trace_dict["rate_stats.mean"] = np.mean(absolute_rates, axis=1)
        trace_dict["rate_stats.variance"] = np.var(absolute_rates, axis=1)
        trace_dict["rate_stats.coefficientOfVariation"] = (
            np.sqrt(trace_dict["rate_stats.variance"]) / trace_dict["rate_stats.mean"]
        )

    trace = pd.DataFrame(trace_dict)
    trace.index.name = "Sample"
    trace.to_csv(trace_out, sep="\t")

    write_tensor_trees(
        topology_file,
        branch_lengths,
        tree_out,
        branch_metadata=(
            dict(rate=result_dict["absolute_rates"]) if model.relaxed_clock() else {}
        ),
    )

    return result_dict


def calculate_coverage(stat_file):
    df = pd.read_table(stat_file)
    return df.truth.between(
        df["95HPDlow"], df["95HPDup"]
    ).mean()  # TODO: Bias and error stats?


def parse_tree_coverage_file(tree_coverage_string):
    df = pd.read_table(io.StringIO(tree_coverage_string.split("\n\n")[1])).set_index(
        "metadata "
    )
    filtered = df[df.index != "posterior"]
    return (filtered["% matches"] / 100.0).to_dict()


def build_method_coverage_table(
    method, coverage_stat_dict, tree_coverage_string, output_file
):
    data_dict = {
        "method": method,
        **{
            stat: calculate_coverage(stat_file)
            for stat, stat_file in coverage_stat_dict.items()
        },
        **parse_tree_coverage_file(tree_coverage_string),
    }
    df = pd.DataFrame([data_dict]).set_index("method")
    df.to_csv(output_file)


def aggregate_coverage_tables(coverage_tables, output_file):
    df = pd.concat([pd.read_csv(file).set_index("method") for file in coverage_tables])
    df.to_csv(output_file)


def build_method_coverage_plot_table(method, coverage_stat_dict, output_file):
    df = pd.concat(
        [
            pd.read_table(stat_file).reset_index().assign(stat=stat)
            for stat, stat_file in coverage_stat_dict.items()
        ]
    ).assign(method=method)
    df.to_csv(output_file, index=False)


def aggregate_coverage_plot_tables(coverage_tables, output_file):
    df = pd.concat([pd.read_csv(file) for file in coverage_tables])
    df.to_csv(output_file, index=False)


def process_variational_key(key: str):
    splits = key.split("_")
    variational_param = splits[-1].split(":")[0]
    var = "_".join(splits[:-1])
    return var, variational_param


def get_formatted_variational_trace(variational_trace):
    parameters = variational_trace.parameters
    for tree_var in ["tree_loc:0", "tree_scale:0"]:
        if tree_var in parameters:
            value = parameters.pop(tree_var)
            parameters[tree_var] = value[:, -1]
    flat_variational_trace, variational_key_mapping = flatten_samples_to_dict(
        parameters,
    )
    variational_params_mapping = {
        key: process_variational_key(key) for key in parameters.keys()
    }
    flat_var_name_mapping = pd.DataFrame(
        [
            (
                i,
                flat_varname,
                varname,
            )
            for varname, flat_varnames in variational_key_mapping.items()
            for i, flat_varname in enumerate(flat_varnames)
        ],
        columns=["var_index", "flat_var_name", "var_name"],
    )
    var_type_mapping = pd.DataFrame(
        [
            (flat_var_name, var_name, var_type)
            for (
                flat_var_name,
                (var_name, var_type),
            ) in variational_params_mapping.items()
        ],
        columns=["var_name", "model_var_name", "var_type"],
    )
    var_metadata = flat_var_name_mapping.merge(var_type_mapping)
    variational_vars_df = pd.DataFrame(
        {key: tensor.numpy() for key, tensor in flat_variational_trace.items()}
    )
    variational_vars_melted = variational_vars_df.reset_index().melt(
        id_vars=("index",), var_name="flat_var_name"
    )
    variational_vars_with_metadata = variational_vars_melted.merge(var_metadata)
    elbo = -variational_trace.loss.numpy()
    elbo_df = (
        pd.DataFrame(dict(value=elbo))
        .reset_index()
        .assign(var_type="loss", model_var_name="elbo", var_index=0)
    )
    return pd.concat(
        [variational_vars_with_metadata[list(elbo_df.columns)], elbo_df]
    ).assign(method="VI")


def get_formatted_ml_trace(ml_trace):
    parameters = dict(ml_trace.parameters)
    tree_trace = parameters.pop("tree")
    ml_vars_dict, ml_params_mapping = flatten_samples_to_dict(parameters)
    raw_tree_stats = calculate_tree_stats("tree", tree_trace)
    tree_stats = dict(tree=raw_tree_stats["tree_height"])
    ml_vars_dict.update(tree_stats)
    flat_var_name_mapping = pd.DataFrame(
        [
            (
                i,
                flat_varname,
                varname,
            )
            for varname, flat_varnames in ml_params_mapping.items()
            for i, flat_varname in enumerate(flat_varnames)
        ]
        + [(0, key, key) for key in tree_stats],
        columns=["var_index", "flat_var_name", "model_var_name"],
    )
    ml_vars_df = pd.DataFrame(
        {key: tensor.numpy() for key, tensor in ml_vars_dict.items()}
    )
    ml_vars_melted = ml_vars_df.reset_index().melt(
        id_vars=("index",), var_name="flat_var_name"
    )
    ml_vars_with_metadata = ml_vars_melted.merge(flat_var_name_mapping).assign(
        var_type="variable"
    )
    ll = ml_trace.log_likelihood.numpy()
    ll_df = (
        pd.DataFrame(dict(value=ll))
        .reset_index()
        .assign(var_type="loss", model_var_name="log_likelihood", var_index=0)
    )
    return pd.concat([ml_vars_with_metadata[list(ll_df.columns)], ll_df]).assign(
        method="ML"
    )


def extract_trace_plot_data(
    variational_trace: VIResults, ml_trace: MLResults, output_file: str
):

    variational_df = get_formatted_variational_trace(variational_trace)
    ml_df = get_formatted_ml_trace(ml_trace)

    res = pd.concat([variational_df, ml_df])
    res.to_csv(output_file, index=False)


def compute_empirical_nucleotide_frequencies(fasta_file: str, output_file: str):
    from treeflow.evolution.substitution.nucleotide.alphabet import A, C, G, T

    frequencies_index_mapping = dict(A=A, C=C, G=G, T=T)

    sequence_dict = parse_fasta(fasta_file)
    sequence_df = pd.DataFrame(
        {taxon: list(sequence) for taxon, sequence in sequence_dict.items()}
    )
    sequences_melted = sequence_df.melt()
    nucleotide_counts = sequences_melted["value"][
        sequences_melted["value"].isin(list(frequencies_index_mapping.keys()))
    ].value_counts()
    empirical_frequencies = nucleotide_counts / nucleotide_counts.sum()
    frequencies_df = pd.DataFrame(
        {
            f"frequencies_{i}": [empirical_frequencies[base]]
            for base, i in frequencies_index_mapping.items()
        }
    )
    frequencies_df.to_csv(output_file, index=False)


def get_runtime_from_benchmark_file(benchmark_file):
    return pd.read_table(benchmark_file).s[0]


def compute_beast_ess(beast_trace, burn_in=0.1):
    burned_in = beast_trace.iloc[int(beast_trace.shape[0] * burn_in) :]
    ess_array = arviz.ess(burned_in.drop(columns="Sample").to_dict(orient="series"))
    return {key: array.values.item() for key, array in dict(ess_array).items()}


def compute_variational_convergence(variational_trace, rtol=0.1, window_size=100):
    loss_series = pd.Series(variational_trace.loss.numpy())
    stds = loss_series.rolling(window=window_size).std()
    final_std = stds.iloc[-1]
    return np.argmax(np.abs((stds - final_std) / final_std) < rtol)


def assemble_timing_data(
    vi_benchmark_file,
    vi_trace_file,
    beast_benchmark_file,
    beast_trace_file,
    output_file,
    beast_burn_in=0.1,
    min_ess=200,
):
    beast_runtime = get_runtime_from_benchmark_file(beast_benchmark_file)
    beast_trace = beast_log_input(beast_trace_file)
    beast_iter = len(beast_trace)
    beast_ess = compute_beast_ess(beast_trace, burn_in=beast_burn_in)
    min_beast_ess = min(beast_ess.values())

    vi_runtime = get_runtime_from_benchmark_file(vi_benchmark_file)
    vi_trace = pickle_input(vi_trace_file)
    vi_iter = len(vi_trace.loss)
    vi_converged_iter = compute_variational_convergence(vi_trace)

    base_vi_df = pd.DataFrame(
        dict(iteration=[0, vi_converged_iter, vi_iter], value=[0.0, 1.0, 1.0])
    )
    vi_df = base_vi_df.assign(
        time=vi_runtime * base_vi_df["iteration"] / vi_iter,
        variable="converged",
        method="vi",
    )

    base_beast_df = pd.DataFrame(
        dict(
            iteration=[0, int(beast_iter * min_ess / min_beast_ess), beast_iter],
            value=[0.0, min_ess, min_beast_ess],
        )
    )
    beast_df = base_beast_df.assign(
        time=beast_runtime * base_beast_df["iteration"] / beast_iter,
        variable="min_ess",
        method="beast",
    )
    res = pd.concat([vi_df, beast_df], axis="rows")
    res.to_csv(output_file)
