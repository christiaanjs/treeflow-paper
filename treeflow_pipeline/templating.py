import jinja2
import pathlib
import treeflow_pipeline.topology_inference

template_env = jinja2.Environment(
    loader=jinja2.PackageLoader('treeflow_pipeline', 'templates')
)

def build_date_string(date_dict):
    return ','.join(['{0}={1}'.format(name, date) for name, date in date_dict.items()])

def build_tree_sim(sim_config, sampling_times, prior_sample, out_file):
    out_path = pathlib.Path(out_file)
    template = template_env.get_template('tree-sim.j2.xml')
    date_trait_string = build_date_string(sampling_times)
    taxon_names = list(sampling_times.keys())
    return template.render(
        pop_size=prior_sample["pop_size"],
        date_trait_string=date_trait_string,
        taxon_names=taxon_names,
        out_file=out_path.parents[0] / "tree-sim.trees"
    )

def build_branch_rate_sim(newick_string, prior_sample, out_file):
    out_path = pathlib.Path(out_file)
    template = template_env.get_template('rate-sim.j2.xml')
    return template.render(
        newick_string=newick_string,
        trace_out_path=out_path.with_suffix('.log'),
        tree_out_path=out_path.with_suffix('.trees'),
        rate_sd=prior_sample["rate_sd"]
    )

def build_sequence_sim(sim_config, newick_string, prior_sample, branch_rates, sampling_times, sequence_length, out_file):
    out_path = pathlib.Path(out_file)
    template = template_env.get_template('sim-seq.j2.xml')
    return template.render(
        out_file=out_path.parents[0] / "sequences.xml",
        sequence_length=sequence_length,
        newick_string=newick_string,
        taxon_names=list(sampling_times.keys()),
        clock_rate=prior_sample["clock_rate"],
        kappa=sim_config["kappa"],
        frequencies=sim_config["frequencies"],
        relaxed_clock=True,
        branch_rates=branch_rates
    )
    

def build_beast_analysis(clock_model, tree_type, sequence_dict, newick_string, init_values, prior_params, beast_config, out_file):
    out_path = pathlib.Path(out_file)
    template = template_env.get_template('beast-analysis.j2.xml')
    date_dict = treeflow_pipeline.topology_inference.parse_dates(sequence_dict)
    date_trait_string = build_date_string(date_dict)
    if tree_type == 'fixed':
        estimate_topology = False
    else:
        estimate_topology = True
    if clock_model not in ['strict', 'relaxed']:
        raise ValueError('Clock model not known: ' + clock_model)
    return template.render(
            estimate_topology=estimate_topology,
            clock_model=clock_model,
            newick_string=newick_string,
            sequence_dict=sequence_dict,
            date_trait_string=date_trait_string,
            trace_out_path=out_path.with_suffix('.log'),
            tree_out_path=out_path.with_suffix('.trees'),
            prior_params=prior_params,
            init_values=init_values,
            **beast_config
        )