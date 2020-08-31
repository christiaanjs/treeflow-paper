import Bio
import Bio.Phylo
import Bio.Phylo.TreeConstruction
import pathlib
import util
import subprocess
import re
import os
import glob

def get_neighbor_joining_tree(msa):
    calculator = Bio.Phylo.TreeConstruction.DistanceCalculator('identity')
    constructor = Bio.Phylo.TreeConstruction.DistanceTreeConstructor(method='nj', distance_calculator=calculator)
    return constructor.build_tree(msa)

LSD_TREE_PATH = 'distance-tree.newick'

def infer_topology_neighbor_joining(input_file, input_format, out_dir):
    with open(input_file) as f:
        sequences = next(Bio.AlignIO.parse(f, format=input_format))
    tree = get_neighbor_joining_tree(sequences)

    for node in tree.get_nonterminals():
        node.name = None

    tree_path = pathlib.Path(out_dir) / LSD_TREE_PATH
    with open(tree_path, 'w') as f:
        Bio.Phylo.write([tree], f, 'newick')
    
    return tree_path     

RAXML_ALPHA_REGEX = r'alpha\[0\]: ([0-9.]+)'
RAXML_RATES_REGEX = r'rates\[0\] ([acgt ]+): ([0-9. ]+) '
RAXML_FREQS_REGEX = r'freqs\[0\]: ([0-9. ]+) '

def parse_raxml_info(raxml_info, estimate_frequencies=True):
    alpha = float(re.search(RAXML_ALPHA_REGEX, raxml_info).group(1))

    rates_res = re.search(RAXML_RATES_REGEX, raxml_info)
    rate_keys_string = rates_res.group(1)
    rate_keys = rate_keys_string.split(' ')
    rate_vals = [float(x) for x in rates_res.group(2).split(' ')]
    rate_dict = dict(zip(rate_keys, rate_vals))

    param_dict = {
        'alpha': alpha,
        'rates': rate_dict
    }

    if estimate_frequencies:
        freqs = [float(x) for x in re.search(RAXML_FREQS_REGEX, raxml_info).group(1).split(' ')]
        param_dict['frequencies'] = freqs

    return param_dict

def infer_topology_raxml(input_file, input_format, out_dir, subst_model=None, seed=123, run_raxml=True, estimate_frequencies=True):
    out_path = pathlib.Path(out_dir)

    if input_format not in ['phylip', 'phylip-sequential', 'phylip-relaxed', 'fasta']:
        raise ValueError('Format not supported by RaxML: ' + input_format)

    raxml_args = ['-p', str(seed), '-m', 'GTRGAMMAX' if estimate_frequencies else 'GTRGAMMA', '-s', input_file, '-n', RAXML_ID, '-w', os.path.abspath(out_dir)]
    if subst_model is not None:
        raxml_args += ['--' + subst_model]

    if run_raxml:
        for raxml_working_file in glob.glob(str(out_path / ('RAxML_*.' + RAXML_ID))):
            os.remove(raxml_working_file)
        subprocess.run(['raxmlHPC'] + raxml_args)

    with open(out_path / ('RAxML_info.' + RAXML_ID)) as f:
        raxml_info = f.read()

    return (out_path / ('RAxML_bestTree.' + RAXML_ID)), parse_raxml_info(raxml_info, estimate_frequencies=estimate_frequencies)
    
LSD_DATE_PATH = 'distance-tree.dates'
LSD_OUT_PATH = 'lsd-tree'

def build_lsd_inputs(out_dir, tree_path, date_trait_dict):
    out_path = pathlib.Path(out_dir)

    date_path = out_path / LSD_DATE_PATH
    with open(date_path, 'w') as f:
        f.write('{0}\n'.format(len(date_trait_dict)))
        for taxon_name, date in date_trait_dict.items():
            f.write('{0} {1}\n'.format(taxon_name, date))

    lsd_args = ['-c'] + util.cmd_kwargs(
        r='a',
        i=tree_path,
        d=date_path,
        o=(out_path / LSD_OUT_PATH)
    )

    return lsd_args

def root_topology(input_file, input_format, out_dir, date_regex, tree_file):
    with open(input_file) as f:
        sequences = next(Bio.AlignIO.parse(f, format=input_format))
    
    date_trait_dict = { record.name: float(re.search(date_regex, record.name).group(0)) for record in sequences }
    lsd_args = build_lsd_inputs(out_dir, tree_file, date_trait_dict)
    subprocess.run(['lsd'] + lsd_args)

    out_path = pathlib.Path(out_dir)
    lsd_out_path = (out_path / LSD_OUT_PATH)

    lsd_date_tree_file = str(lsd_out_path) + '.date.nexus'
    lsd_distance_tree_file = str(lsd_out_path) + '.nexus'

    with open(lsd_date_tree_file) as f:
        date_tree = next(Bio.Phylo.parse(f, format='nexus'))

    with open(lsd_distance_tree_file) as f:
        distance_tree = next(Bio.Phylo.parse(f, format='nexus'))

    estimated_rate = distance_tree.total_branch_length() / date_tree.total_branch_length()

    return lsd_date_tree_file, estimated_rate 

def estimate_pop_size(tree_file, tree_format):
    with open(tree_file) as f:
        tree = next(Bio.Phylo.parse(f, tree_format))
    
    depths = tree.depths()
    leaf_depths = [depths[leaf] for leaf in tree.get_terminals()]
    return  min(leaf_depths) * 2.0

def get_taxon_count(tree_file, tree_format):
    with open(tree_file) as f:
        tree = next(Bio.Phylo.parse(f, tree_format))

    return tree.count_terminals()

EPSILON = 1e-4
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

def convert_tree(input_file, input_format, output_format, strip_data=False, allow_zero_branches=True, epsilon=EPSILON):
    with open(input_file) as f:
        trees = list(Bio.Phylo.parse(input_file, input_format))

    if strip_data:
        for tree in trees:
            for clade in tree.find_clades():
                clade.comment = None

    if not allow_zero_branches:
        adjust_zero_branches(tree.clade, epsilon=epsilon)

    output_path = pathlib.Path(input_file).with_suffix('.' + output_format)
    with open(output_path, 'w') as f:
        Bio.Phylo.write(trees, f, output_format)
    return str(output_path)