import jinja2
import pathlib
import treeflow.model
import treeflow_pipeline.model
import treeflow_pipeline.topology_inference
import xml.etree.ElementTree as ET
import collections

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
        rate_sd=prior_sample["rate_sd"],
        clock_rate=prior_sample["clock_rate"]
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

def x2s(x):
    return ET.tostring(x, encoding="unicode")

def get_state_tag(name, dist_dict, init_value):
    dist_name, params = next(iter(dist_dict.items()))
    support = treeflow.model.distribution_class_supports[treeflow_pipeline.model.dists[dist_name]]
    element = ET.Element("parameter", attrib=dict(name="stateNode", id=name))
    element.text = str(init_value)
    if support == "nonnegative":
        element.set("lower", str(0.0))
    elif support != "real":
        raise ValueError(f"Support {support} not implemented for distribution {dist_name}")
    return x2s(element)

dist_functions = dict(
    lognormal = lambda loc, scale: ("LogNormal", dict(M=str(loc), S=str(scale)))
)

def wrap_in_prior(name, dist_tag):
    prior_tag = ET.Element("prior", name="distribution", x=f"@{name}")
    prior_tag.insert(0, dist_tag)
    return prior_tag

def get_prior_tag(name, dist_dict):
    dist_name, params = next(iter(dist_dict.items()))
    try:
        tag_name, params_attrib = dist_functions[dist_name](**params)
    except KeyError:
        raise ValueError(f"Unknown distribution {dist_name} for {name}")
    dist_tag = ET.Element(tag_name, attrib=params_attrib, name="distr")
    return x2s(wrap_in_prior(name, dist_tag))

def resolve_param_value(name, model, init_value):
    if model == "fixed":
        if isinstance(init_value, collections.Iterable): # TODO: Case of 0d Numpy as array
            return " ".join(map(str, init_value))
        else:
            return str(init_value)
    else:
        return f"@{name}"

def get_rate_prior_tag(clock_model, params):
    if clock_model == "strict":
        return None
    elif clock_model == "relaxed_lognormal":
        tag_name, params_attrib = dist_functions["lognormal"](1.0, resolve_param_value("rate_sd", params["rate_sd"], 1.0)) # TODO: Init values for rates
        dist_tag = ET.Element(tag_name, attrib=params_attrib, meanInRealSpace="true", name="distr", id="rate_prior")
        return x2s(wrap_in_prior("rates", dist_tag))
    else:
        raise ValueError(f"Unknown clock model: {clock_model}")

def get_tree_prior_tag(tree_model, params, init_values):
    if tree_model == "coalescent_constant":
        population_tag = ET.Element("populationModel", spec="ConstantPopulation", popSize=resolve_param_value("pop_size", params["pop_size"], init_values["pop_size"])) # TODO: Case of missing init value
        intervals_tag = ET.Element("treeIntervals", spec="TreeIntervals", tree="@tree")
        dist_tag = ET.Element("distribution", spec="Coalescent")
        dist_tag.insert(0, population_tag)
        dist_tag.insert(1, intervals_tag)
        return x2s(dist_tag)
    else:
        raise ValueError(f"Unknown tree model: {tree_model}")

subst_model_specs = dict(
    hky="HKY"
)

def get_subst_model_tag(subst_model, params, init_values):
    frequencies_tag = ET.Element(
        "frequencies",
        spec="Frequencies",
        frequencies=resolve_param_value("frequencies", params["frequencies"], init_values["frequencies"])
    )
    try:
        spec = subst_model_specs[subst_model]
    except KeyError:
        raise ValueError(f"Unknown substitution model {subst_model}")
    subst_tag = ET.Element(
        "substModel",
        attrib={ param: resolve_param_value(param, model, init_values[param]) for param, model in params.items() if param != "frequencies" },
        spec=spec
    )
    subst_tag.insert(0, frequencies_tag)
    return subst_tag

def get_site_model_tag(site_model, params, init_values, subst_model_tag):
    if site_model == "none":
        site_tag = ET.Element("siteModel", spec="SiteModel", mutationRate=str(1.0), shape=str(1.0), proportionInvariant=str(0.0))
        site_tag.insert(0, subst_model_tag)
        return x2s(site_tag)
    else:
        raise ValueError(f"Unknown site model: {site_model}")

def get_branch_rate_model_tag(clock_model, params, init_values):
    attrib = {
        "id": "branch_rate_model",
        "clock.rate": resolve_param_value("clock_rate", params["clock_rate"], init_values["clock_rate"])
    }
    if clock_model == "strict":
        attrib["spec"] = "StrictClockModel"
    elif clock_model == "relaxed_lognormal":
        attrib.update(dict(
            spec="UCRelaxedClockModel",
            rates="@rates",
            tree="@tree",
            distr="@rate_prior"
        ))
    else:
        raise ValueError(f"Unknown clock model: {clock_model}")
    return x2s(ET.Element("branchRateModel", attrib=attrib))

op_functions = dict(
    nonnegative=lambda scale_factor=0.75: dict(spec="ScaleOperator", scaleFactor=str(scale_factor))
)

def get_operator_tag(name, dist_dict, weight=3.0, **op_kwargs): # TODO: Deal with special 
    dist_name, params = next(iter(dist_dict.items()))
    support = treeflow.model.distribution_class_supports[treeflow_pipeline.model.dists[dist_name]]
    try:
        attrib = op_functions[support](**op_kwargs)
    except KeyError:
        raise ValueError(f"Unknown support {support} for distribution {dist_name} for {name}")
    return x2s(ET.Element("operator", attrib, parameter=f"@{name}", weight=str(weight)))
    
def get_clock_operator_tags(clock_model, params):
    clock_rate_dist = next(iter(params["clock_rate"]))
    ops = []

    if clock_rate_dist != "fixed":
        ops.append(x2s(ET.Element("operator",
            spec="UpDownOperator",
            scaleFactor=str(0.75),
            weight=str(3.0),
            up="@clock_rate",
            down="@tree"
        )))
    
    # TODO: Relaxed clock rate operators
    return ops


def get_log_tag(name):
    return x2s(ET.Element("log", idref=name))

def build_beast_analysis(sequence_dict, newick_string, init_values, model, beast_config, out_file, estimate_topology=False):
    out_path = pathlib.Path(out_file)
    template = template_env.get_template('beast-analysis.j2.xml')
    date_dict = treeflow_pipeline.topology_inference.parse_dates(sequence_dict)
    date_trait_string = build_date_string(date_dict)

    # TODO: Handle case when init_values missing
    # TODO: Starting value for rates?
    state_tags = [get_state_tag(name, dist_dict, init_values[name]) for name, dist_dict in model.free_params().items()]
    prior_tags = [get_prior_tag(name, dist_dict) for name, dist_dict in model.free_params().items()]
    tree_prior_tag = get_tree_prior_tag(model.tree_model, model.tree_params, init_values)
    rate_prior_tag = get_rate_prior_tag(model.clock_model, model.clock_params)
    subst_model_tag = get_subst_model_tag(model.subst_model, model.subst_params, init_values)
    site_model_tag = get_site_model_tag(model.site_model, model.site_params, init_values, subst_model_tag)
    branch_rate_model_tag = get_branch_rate_model_tag(model.clock_model, model.clock_params, init_values)

    operator_tags = ([get_operator_tag(name, dist_dict) for name, dist_dict in model.free_params().items()]
        + get_clock_operator_tags(model.clock_model, model.clock_params))
    log_tags = [get_log_tag(name) for name in model.free_params()]

    return template.render(
            estimate_topology=estimate_topology,
            relaxed=model.clock_model in treeflow_pipeline.model.RELAXED_CLOCK_MODELS,
            newick_string=newick_string,
            sequence_dict=sequence_dict,
            date_trait_string=date_trait_string,
            trace_out_path=out_path.with_suffix('.log'),
            tree_out_path=out_path.with_suffix('.trees'),
            state_tags=state_tags,
            prior_tags=prior_tags,
            tree_prior_tag=tree_prior_tag,
            rate_prior_tag=rate_prior_tag,
            site_model_tag=site_model_tag,
            branch_rate_model_tag=branch_rate_model_tag,
            operator_tags=operator_tags,
            log_tags=log_tags,
            **beast_config
        )