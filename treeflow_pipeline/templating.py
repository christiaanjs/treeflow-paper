import jinja2
import treeflow_pipeline.topology_inference

template_env = jinja2.Environment(
    loader=jinja2.PackageLoader('treeflow-paper', 'templates')
)

def build_beast_analysis(sequence_dict, newick_string, init_values, prior_params, out_file):
    out_path = pathlib.Path(out_file)
    template = template_env.get_template('beast-analysis.j2.xml')
    return template.render(
            newick_string=newick_string,
            sequence_dict=sequence_dict,
            date_trait_string=treeflow_pipeline.topology_inference.parse_dates(sequence_dict),
            trace_out_path=out_path.with_suffix('.log'),
            tree_out_path=out_path.with_suffix('.trees'),
            prior_params=prior_params,
            init_values=init_values
        )