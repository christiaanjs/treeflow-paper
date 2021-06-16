import dendropy
import numpy as np
import treeflow.tree_processing
import treeflow.sequences
import treeflow_pipeline.model
import pandas as pd
import io


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


def tensor_to_dendro(
    topology, taxon_namespace, taxon_names, branch_lengths, branch_metadata={}
):
    taxon_count = len(taxon_names)
    leaves = [
        dendropy.Node(taxon=taxon_namespace.get_taxon(name)) for name in taxon_names
    ]
    nodes = leaves + [dendropy.Node() for _ in range(taxon_count - 1)]
    for i, node in enumerate(nodes[:-1]):
        node.edge_length = branch_lengths[i]
        for key, value in branch_metadata.items():
            node.annotations[key] = value[i]
        parent = nodes[topology["parent_indices"][i]]
        parent.add_child(node)
    return dendropy.Tree(
        taxon_namespace=taxon_namespace, seed_node=nodes[-1], is_rooted=True
    )


from dendropy.dataio import nexusprocessing


class CustomNewickWriter(dendropy.dataio.newickwriter.NewickWriter):
    def _write_node_body(self, node, out):
        out.write(self._render_node_tag(node))
        if not self.suppress_annotations:
            node_annotation_comments = (
                nexusprocessing.format_item_annotations_as_comments(
                    node,  # Place node annotations before colon
                    nhx=self.annotations_as_nhx,
                    real_value_format_specifier=self.real_value_format_specifier,
                )
            )
            out.write(node_annotation_comments)
        if node.edge and node.edge.length != None and not self.suppress_edge_lengths:
            out.write(":{}".format(self.edge_label_compose_fn(node.edge)))
        if not self.suppress_annotations:
            edge_annotation_comments = (
                nexusprocessing.format_item_annotations_as_comments(
                    node.edge,
                    nhx=self.annotations_as_nhx,
                    real_value_format_specifier=self.real_value_format_specifier,
                )
            )
            out.write(edge_annotation_comments)
        out.write(self._compose_comment_string(node))
        out.write(self._compose_comment_string(node.edge))


class CustomNexusWriter(dendropy.dataio.nexuswriter.NexusWriter):
    def __init__(self, **kwargs):
        super(CustomNexusWriter, self).__init__(**kwargs)

        kwargs_to_preserve = [
            "unquoted_underscores",
            "preserve_spaces",
            "annotations_as_nhx",
            "suppress_annotations",
            "suppress_item_comments",
        ]

        newick_kwargs = dict(
            unquoted_underscores=self.unquoted_underscores,
            preserve_spaces=self.preserve_spaces,
            annotations_as_nhx=self.annotations_as_nhx,
            suppress_annotations=self.suppress_annotations,
            suppress_item_comments=self.suppress_item_comments,
        )
        self._newick_writer = CustomNewickWriter(**newick_kwargs)


def fit_successful(variational_fit):
    return np.isfinite(variational_fit["loss"]).all()


NUMERICAL_ISSUE_N = 2


def write_tensor_trees(topology_file, branch_lengths, output_file, branch_metadata={}):
    taxon_namespace = dendropy.Tree.get(
        path=topology_file, schema="newick", preserve_underscores=True
    ).taxon_namespace
    tree, taxon_names = treeflow.tree_processing.parse_newick(topology_file)
    trees = dendropy.TreeList(
        [
            tensor_to_dendro(
                tree["topology"],
                taxon_namespace,
                taxon_names,
                branch_lengths[i],
                branch_metadata={
                    key: value[i] for key, value in branch_metadata.items()
                },
            )
            for i in range(branch_lengths.shape[0])
        ],
        taxon_namespace=taxon_namespace,
    )

    writer = CustomNexusWriter(unquoted_underscores=True)
    with open(output_file, "w") as f:
        writer.write_tree_list(trees, f)


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
