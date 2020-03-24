import Bio
import Bio.Phylo
import Bio.Phylo.TreeConstruction
import pathlib
import util
import subprocess
import re

def get_neighbor_joining_tree(msa):
    calculator = Bio.Phylo.TreeConstruction.DistanceCalculator('identity')
    constructor = Bio.Phylo.TreeConstruction.DistanceTreeConstructor(method='nj', distance_calculator=calculator)
    return constructor.build_tree(msa)

LSD_TREE_PATH = 'distance-tree.newick'
LSD_DATE_PATH = 'distance-tree.dates'
LSD_OUT_PATH = 'lsd-tree'

def build_lsd_inputs(out_dir, distance_tree, date_trait_dict):
    out_path = pathlib.Path(out_dir)

    tree_path = out_path / LSD_TREE_PATH
    with open(tree_path, 'w') as f:
        Bio.Phylo.NewickIO.write([distance_tree], f)

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
            
def infer_topology(input_file, input_format, out_dir, date_regex):
    with open(input_file) as f:
        sequences = next(Bio.AlignIO.parse(f, format=input_format))
    
    date_trait_dict = { record.name: float(re.search(date_regex, record.name).group(0)) for record in sequences }
    nj_tree = get_neighbor_joining_tree(sequences)
    lsd_args = build_lsd_inputs(out_dir, nj_tree, date_trait_dict)
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

    return estimated_rate, lsd_date_tree_file
