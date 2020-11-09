import numpy as np
import treeflow_pipeline.model
import xml
import Bio.Seq
import Bio.SeqIO

def rng(sim_config, hash=0):
    return  np.random.default_rng(sim_config['seed'] + abs(hash))

def get_sampling_times(sim_config, taxon_count):
    times = np.linspace(0, sim_config["sampling_window"], taxon_count)
    return { 'T{0}_{1}'.format(i + 1, time): float(time) for i, time in enumerate(times) }

def sample_prior(sampling_times, model, seed):
    prior = treeflow_pipeline.model.get_phylo_prior(list(sampling_times.values()), model)
    sample = prior.sample(seed=seed)
    return { key: sample[key].numpy().item() for key in model.free_params() }

def parse_branch_rates(df):
    return df.iloc[0, 1:].tolist()

def convert_simulated_sequences(input_file, output_file, output_format):
    seq_xml_root = xml.etree.ElementTree.parse(input_file)
    records = [Bio.SeqIO.SeqRecord(Bio.Seq.Seq(tag.attrib['value'], Bio.Alphabet.generic_dna), tag.attrib['taxon'], description="") for tag in seq_xml_root.findall('./sequence')]
    with open(output_file, 'w') as f:
        Bio.SeqIO.write(records, f, 'fasta')