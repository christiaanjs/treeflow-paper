import Bio
import Bio.Phylo
import Bio.Phylo.TreeConstruction
import pathlib
import treeflow_pipeline.util as util
import treeflow_pipeline.model
import subprocess
import re
import os
import glob


def get_neighbor_joining_tree(msa):
    calculator = Bio.Phylo.TreeConstruction.DistanceCalculator("identity")
    constructor = Bio.Phylo.TreeConstruction.DistanceTreeConstructor(
        method="nj", distance_calculator=calculator
    )
    return constructor.build_tree(msa)


LSD_TREE_PATH = "distance-tree.newick"


def infer_topology_neighbor_joining(
    input_file, input_format, out_dir
):  # TODO: Refactor
    with open(input_file) as f:
        sequences = next(Bio.AlignIO.parse(f, format=input_format))
    tree = get_neighbor_joining_tree(sequences)

    for node in tree.get_nonterminals():
        node.name = None

    tree_path = pathlib.Path(out_dir) / LSD_TREE_PATH
    with open(tree_path, "w") as f:
        Bio.Phylo.write([tree], f, "newick")

    return tree_path


def build_raxml_command(
    alignment_file, working_directory, seed, site_model, subst_model
):
    args = []
    if subst_model == "hky":
        args.append("--HKY85")
    elif subst_model == "jc":
        args.append("--JC69")
    else:
        raise ValueError(
            "Unsupported substitution model for RAxML: {0}".format(subst_model)
        )
    if site_model == "none":
        args.append("-V")
    elif site_model == "discrete_gamma":
        pass
    else:
        raise ValueError("Unsupported site model for RAxML: {0}".format(site_model))
    kwargs = dict(p=seed, n="raxml", s=alignment_file, w=working_directory, m="GTRCAT" if subst_model == "jc" else "GTRCATX")
    return util.build_shell_command("raxmlHPC", args, kwargs)


RAXML_ALPHA_REGEX = r"alpha\[0\]: ([0-9.]+)"
RAXML_RATES_REGEX = r"rates\[0\] ([acgt ]+): ([0-9. ]+) "
RAXML_FREQS_REGEX = r"freqs\[0\]: ([0-9. ]+) "


def parse_raxml_info(filename, estimate_frequencies=True):
    with open(filename) as f:
        raxml_info = f.read()

    alpha = float(re.search(RAXML_ALPHA_REGEX, raxml_info).group(1))

    rates_res = re.search(RAXML_RATES_REGEX, raxml_info)
    rate_keys_string = rates_res.group(1)
    rate_keys = rate_keys_string.split(" ")
    rate_vals = [float(x) for x in rates_res.group(2).split(" ")]
    rate_dict = dict(zip(rate_keys, rate_vals))

    param_dict = {"alpha": alpha, "rates": rate_dict}

    if estimate_frequencies:
        freqs = [
            float(x)
            for x in re.search(RAXML_FREQS_REGEX, raxml_info).group(1).split(" ")
        ]
        sum_freqs = sum(freqs)
        param_dict["frequencies"] = [
            freq / sum_freqs for freq in freqs
        ]  # Normalise for rounding error

    return param_dict


def infer_topology_raxml(
    input_file,
    input_format,
    out_dir,
    subst_model=None,
    seed=123,
    run_raxml=True,
    estimate_frequencies=True,
):
    out_path = pathlib.Path(out_dir)

    if input_format not in ["phylip", "phylip-sequential", "phylip-relaxed", "fasta"]:
        raise ValueError("Format not supported by RaxML: " + input_format)

    raxml_args = [
        "-p",
        str(seed),
        "-m",
        "GTRGAMMAX" if estimate_frequencies else "GTRGAMMA",
        "-s",
        input_file,
        "-n",
        RAXML_ID,
        "-w",
        os.path.abspath(out_dir),
    ]
    if subst_model is not None:
        raxml_args += ["--" + subst_model]

    if run_raxml:
        for raxml_working_file in glob.glob(str(out_path / ("RAxML_*." + RAXML_ID))):
            os.remove(raxml_working_file)
        subprocess.run(["raxmlHPC"] + raxml_args)

    raxml_info = parse_raxml_info(out_path / ("RAxML_info." + RAXML_ID))

    return (out_path / ("RAxML_bestTree." + RAXML_ID)), parse_raxml_info(
        raxml_info, estimate_frequencies=estimate_frequencies
    )


def get_starting_values_raxml(raxml_info_file, subst_model):
    estimate_frequencies = subst_model != "jc"
    raxml_info = parse_raxml_info(raxml_info_file, estimate_frequencies=estimate_frequencies)
    res = dict()
    if estimate_frequencies:
        res["frequencies"] = raxml_info["frequencies"]
    if subst_model == "hky":
        res["kappa"] = raxml_info["rates"]["ag"]
    elif subst_model != "jc":
        raise ValueError(
            "Unsupported substitution model for RAxML: {0}".format(subst_model)
        )
    return res


LSD_DATE_PATH = "distance-tree.dates"
LSD_OUT_PATH = "lsd-tree"


def parse_dates(names):
    return {name: float(name.split("_")[-1]) for name in names}


def parse_leaf_heights(names):
    dates = parse_dates(names)
    max_date = max(dates.values())
    return {name: (max_date - date) for name, date in dates.items()}


def build_lsd_date_file(sequence_dict, output_file):
    date_trait_dict = parse_dates(sequence_dict)

    with open(output_file, "w") as f:
        f.write("{0}\n".format(len(date_trait_dict)))
        for taxon_name, date in date_trait_dict.items():
            f.write("{0} {1}\n".format(taxon_name, date))


def build_lsd_inputs(
    input_file, input_format, out_dir, tree_path
):  # TODO: Remove input format
    out_path = pathlib.Path(out_dir)

    date_path = out_path / LSD_DATE_PATH
    build_lsd_date_file(util.sequence_input(input_file, input_format), date_path)

    lsd_args = ["-c"] + util.cmd_kwargs(
        r="a", i=tree_path, d=date_path, o=(out_path / LSD_OUT_PATH)
    )

    return lsd_args


def estimate_rate(
    date_tree_file,
    distance_tree_file,
    date_tree_format="nexus",
    distance_tree_format="nexus",
):
    with open(date_tree_file) as f:
        date_tree = next(Bio.Phylo.parse(f, format=date_tree_format))

    with open(distance_tree_file) as f:
        distance_tree = next(Bio.Phylo.parse(f, format=distance_tree_format))

    return distance_tree.total_branch_length() / date_tree.total_branch_length()


def root_topology(
    input_file, input_format, out_dir, date_regex, tree_file
):  # TODO: Remove date_regex argument
    lsd_args = build_lsd_inputs(input_file, input_format, out_dir, tree_file)
    subprocess.run(["lsd"] + lsd_args)

    out_path = pathlib.Path(out_dir)
    lsd_out_path = out_path / LSD_OUT_PATH

    lsd_date_tree_file = str(lsd_out_path) + ".date.nexus"
    lsd_distance_tree_file = str(lsd_out_path) + ".nexus"

    estimated_rate = estimate_rate(lsd_date_tree_file, lsd_distance_tree_file)

    return lsd_date_tree_file, estimated_rate


def estimate_pop_size(tree_file, tree_format):
    with open(tree_file) as f:
        tree = next(Bio.Phylo.parse(f, tree_format))

    depths = tree.depths()
    leaf_depths = [depths[leaf] for leaf in tree.get_terminals()]
    return min(leaf_depths) * 2.0

def estimate_birth_rate(tree_file, tree_format):
    # Kendal-Moran estimator of the speciation rate 
    with open(tree_file) as f:
        tree = next(Bio.Phylo.parse(f, tree_format))

    n = tree.count_terminals()
    branch_sum = tree.total_branch_length()
    return (n - 2) / branch_sum

def get_taxon_count(tree_file, tree_format):
    with open(tree_file) as f:
        tree = next(Bio.Phylo.parse(f, tree_format))

    return tree.count_terminals()


def get_starting_values_lsd(date_tree_file, distance_tree_file, lsd_output_format, tree_model):
    res = dict(
        clock_rate=estimate_rate(
            date_tree_file,
            distance_tree_file,
            date_tree_format=lsd_output_format,
            distance_tree_format=lsd_output_format,
        ),
    )
    if tree_model == "constant_coalescent":
        res["pop_size"] = estimate_pop_size(date_tree_file, lsd_output_format)
    elif tree_model == "yule":
        res["birth_rate"] = estimate_birth_rate(date_tree_file, lsd_output_format)
    else:
        raise ValueError(f"Unknown tree model: {tree_model}")
    return res

def finalise_starting_values(
    tree_starting_values, rooting_starting_values, clock_model
):
    clock_model_defaults = treeflow_pipeline.model.get_non_rate_defaults(clock_model)
    return {**clock_model_defaults, **tree_starting_values, **rooting_starting_values}


EPSILON = 1e-6


def adjust_zero_branches(clade, epsilon=EPSILON):
    if not clade.is_terminal():
        for subclade in clade.clades:
            adjust_zero_branches(subclade, epsilon=epsilon)
        for subclade in clade.clades:
            if subclade.branch_length < epsilon:
                diff = epsilon - subclade.branch_length
                clade.branch_length -= diff
                for subclade_2 in clade.clades:
                    subclade_2.branch_length += diff


def preorder_traversal(tree):
    stack = [tree.clade]
    while len(stack):
        clade = stack.pop()
        yield clade
        for child in clade:
            stack.append(child)


def postorder_traversal(tree):
    stack = [tree.clade]
    output = []
    while len(stack):
        clade = stack.pop()
        output.append(clade)
        for child in clade:
            stack.append(child)
    while len(output):
        yield output.pop()


def get_node_heights(tree):
    distances = {}
    for node in preorder_traversal(tree):
        if not node.is_terminal():
            if node is tree.clade:
                distances[node] = 0.0
            for child in node:
                distances[child] = distances[node] + child.branch_length
    max_distance = max(distances.values())
    return {node: max_distance - distance for node, distance in distances.items()}


def fix_leaf_dates(tree):
    leaf_heights = parse_leaf_heights([node.name for node in tree.get_terminals()])
    heights = get_node_heights(tree)
    for node in postorder_traversal(tree):
        if node.is_terminal():
            heights[node] = leaf_heights[node.name]
        else:
            min_height = max([heights[child] for child in node])
            if heights[node] < min_height:
                heights[node] = min_height
            for child in node:
                child.branch_length = heights[node] - heights[child]


def convert_tree(
    input_file,
    input_format,
    output_file,
    output_format,
    strip_data=False,
    allow_zero_branches=True,
    fix_dates=False,
    epsilon=EPSILON,
):
    with open(input_file) as f:
        trees = list(Bio.Phylo.parse(input_file, input_format))

    if strip_data:
        for tree in trees:
            for clade in tree.find_clades():
                clade.comment = None

    if fix_dates:
        for tree in trees:
            fix_leaf_dates(tree)

    if not allow_zero_branches:
        for tree in trees:
            adjust_zero_branches(tree.clade, epsilon=epsilon)

    with open(output_file, "w") as f:
        Bio.Phylo.write(trees, f, output_format)
