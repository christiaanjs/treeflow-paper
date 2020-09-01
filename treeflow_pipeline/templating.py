import jinja2
import pathlib
import treeflow_pipeline.topology_inference

template_env = jinja2.Environment(
    loader=jinja2.PackageLoader('treeflow_pipeline', 'templates')
)

def build_beast_analysis(clock_model, tree_type, sequence_dict, newick_string, init_values, prior_params, beast_config, out_file):
    out_path = pathlib.Path(out_file)
    template = template_env.get_template('beast-analysis.j2.xml')
    date_dict = treeflow_pipeline.topology_inference.parse_dates(sequence_dict)
    date_trait_string = ','.join(['{0}={1}'.format(name, date) for name, date in date_dict.items()])
    if tree_type == 'fixed':
        estimate_topology = False
    else:
        estimate_topology = True
    if clock_model != 'strict':
        raise ValueError('Clock model not known: ' + clock_model)
    return template.render(
            estimate_topology=estimate_topology,
            newick_string=newick_string,
            sequence_dict=sequence_dict,
            date_trait_string=date_trait_string,
            trace_out_path=out_path.with_suffix('.log'),
            tree_out_path=out_path.with_suffix('.trees'),
            prior_params=prior_params,
            init_values=init_values,
            **beast_config
        )