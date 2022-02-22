import numpy as np
import treeflow_pipeline.model
import xml
import pandas as pd
import Bio.Seq
import Bio.SeqIO
from treeflow.tree.io import parse_newick


def rng(sim_config, hash=0):
    return np.random.default_rng(sim_config["seed"] + abs(hash))


def get_sampling_times(sim_config, taxon_count):
    times = np.linspace(0, sim_config["sampling_window"], taxon_count)
    return {"T{0}_{1}".format(i + 1, time): float(time) for i, time in enumerate(times)}


def sample_prior(sampling_times, model, seed):
    prior = treeflow_pipeline.model.get_phylo_prior(
        list(sampling_times.values()), model
    )
    sample = prior.sample(seed=seed)
    return {key: sample[key].numpy().item() for key in model.free_params()}


def build_sim_trace(tree_file, prior_sample, out_file, rate_trace=None):
    tree = parse_newick(tree_file)
    branch_lengths = tree.branch_lengths

    values = {
        "tree.height": tree["heights"][-1],
        "tree.treeLength": np.sum(branch_lengths),
        **prior_sample,
    }

    if rate_trace is not None:  # TODO: Safe to assume clock_rate?
        rates = pd.read_table(rate_trace, comment="#")
        if "clock_rate" in prior_sample:
            rate_values = (
                rates[rates.columns[1:]].values[0] * prior_sample["clock_rate"]
            )
        else:
            rate_values = rates[rates.columns[1:]].values[0]
        values["rate_stats.mean"] = np.mean(rate_values)
        values["rate_stats.variance"] = np.var(rate_values)
        values["rate_stats.coefficientOfVariation"] = (
            np.sqrt(values["rate_stats.variance"]) / values["rate_stats.mean"]
        )

    pd.DataFrame([values]).to_csv(out_file, index=False, sep="\t")


def aggregate_sim_traces(sim_trace_files, out_file):
    traces = [pd.read_table(x) for x in sim_trace_files]
    aggregate = pd.concat(traces, axis=0, ignore_index=True)
    aggregate.index.name = "Sample"
    aggregate.to_csv(out_file, sep="\t")


def parse_branch_rates(df):
    return df.iloc[0, 1:].tolist()


def convert_simulated_sequences(input_file, output_file, output_format):
    seq_xml_root = xml.etree.ElementTree.parse(input_file)
    records = [
        Bio.SeqIO.SeqRecord(
            Bio.Seq.Seq(tag.attrib["value"], Bio.Alphabet.generic_dna),
            tag.attrib["taxon"],
            description="",
        )
        for tag in seq_xml_root.findall("./sequence")
    ]
    with open(output_file, "w") as f:
        Bio.SeqIO.write(records, f, "fasta")


import dendropy


def aggregate_trees(input_files, input_format, output_file, output_format):
    trees = [
        dendropy.Tree.get(
            path=input_file,
            schema="nexus",
            rooting="default-rooted",
            preserve_underscores=True,
        )
        for input_file in input_files
    ]
    taxon_namespace = trees[0].taxon_namespace
    tree_list = dendropy.TreeList(trees, taxon_namespace=taxon_namespace)
    tree_list.write(path=output_file, schema=output_format)
