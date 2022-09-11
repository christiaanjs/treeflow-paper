import tensorflow as tf
import tensorflow_probability as tfp

tfd = tfp.distributions
import treeflow
from treeflow.model.phylo_model import (
    PhyloModel as Model,
    parse_model,
    RELAXED_CLOCK_MODELS,
    prior_distribution_classes as dists,
    phylo_model_to_joint_distribution,
)
from treeflow.vi.util import VIResults
from treeflow.evolution.seqio import Alignment
from treeflow.tree.rooted.tensorflow_rooted_tree import convert_tree_to_tensor
from treeflow.tree.io import parse_newick
from treeflow.model.approximation import get_fixed_topology_mean_field_approximation


APPROX_MODELS = [
    "mean_field",
    "scaled",
    "scaled_shrinkage",
    "scaled_conjugate",
    "scaled_shrinkage_conjugate",
    "scaled_shrinkage_local_conjugate",
]


def get_non_rate_defaults(clock_model):
    if clock_model == "relaxed_lognormal":
        return dict(rate_sd=1.0)
    elif clock_model == "strict":
        return {}
    else:
        raise ValueError(f"Clock model not known: {clock_model}")


def cast(x):
    return tf.convert_to_tensor(x, treeflow.DEFAULT_FLOAT_DTYPE_TF)


def get_dist(dist_dict):
    dist_name, params = parse_model(dist_dict)  # TODO: Handle case of fixed values
    try:
        return dists[dist_name](
            **{key: cast(val) for key, val in params.items()}
        )  # TODO: Non-float params?
    except KeyError:
        raise ValueError("Distribution not known: {0}".format(dist_name))


def get_phylo_prior(sampling_times, model):
    sampling_times = cast(sampling_times)
    taxon_count = sampling_times.shape[-1]
    model_dict = {}

    if model.tree_model == "coalescent_constant":
        model_dict["pop_size"] = get_dist(
            model.tree_params["pop_size"]
        )  # TODO: Make this more automatic
        model_dict["tree"] = lambda pop_size: treeflow.coalescent.ConstantCoalescent(
            taxon_count, pop_size, sampling_times
        )
    else:
        raise ValueError("Tree prior not known: {0}".format(model.tree_model))

    if model.clock_model == "relaxed_lognormal":
        model_dict["clock_rate"] = get_dist(model.clock_params["clock_rate"])
        model_dict["rate_sd"] = get_dist(model.clock_params["rate_sd"])
        model_dict["rates"] = lambda rate_sd: tfd.Sample(
            tfd.LogNormal(loc=-(rate_sd**2.0) / 2, scale=rate_sd),
            sample_shape=2 * taxon_count - 2,
        )
    elif model.clock_model == "relaxed_lognormal_conjugate":
        conjugate_prior_dict = treeflow.priors.get_normal_conjugate_prior_dict(
            **{
                key: cast(value)
                for key, value in model.clock_params["rate_loc_precision"][
                    "normalgamma"
                ].items()
            }
        )
        model_dict["rate_precision"] = conjugate_prior_dict["precision"]
        model_dict["rate_loc"] = lambda rate_precision: conjugate_prior_dict["loc"](
            rate_precision
        )
        model_dict["rates"] = lambda rate_loc, rate_precision: tfd.Sample(
            tfd.LogNormal(loc=rate_loc, scale=tf.math.sqrt(1.0 / rate_precision)),
            sample_shape=2 * taxon_count - 2,
        )
    elif model.clock_model == "strict":
        model_dict["clock_rate"] = get_dist(model.clock_params["clock_rate"])
    else:
        raise ValueError("Clock model not known: {0}".format(model.clock_model))

    # TODO: Substitution model
    # TODO: Site model

    return tfd.JointDistributionNamed(model_dict)


optimizers = dict(adam=tf.optimizers.Adam, sgd=tf.optimizers.SGD)


def fit_surrogate_posterior(
    log_p,
    q,
    vi_config,
    trace_fn=lambda x: (x.loss, x.parameters),
    vi_kwargs={},
    function_mode=True,
):
    optimizer = optimizers[vi_config["optimizer"]](**vi_config["optimizer_kwargs"])
    if function_mode:
        return tfp.vi.fit_surrogate_posterior(
            log_p,
            q,
            optimizer,
            vi_config["num_steps"],
            trace_fn=trace_fn,
            seed=vi_config["seed"],
            sample_size=vi_config["sample_size"],
            **vi_kwargs,
        )
    else:
        return treeflow.vi.fit_surrogate_posterior(
            log_p,
            q,
            optimizer,
            vi_config["num_steps"],
            trace_fn=trace_fn,
            seed=vi_config["seed"],
            sample_size=vi_config["sample_size"],
            **vi_kwargs,
        )


def get_likelihood(newick_file, fasta_file, starting_values, model, vi_config):
    if model.subst_model == "hky":
        subst_model = treeflow.substitution_model.HKY()
        if (
            model.subst_params["kappa"] == "fixed"
            and model.subst_params["frequencies"] == "fixed"
        ):
            return treeflow.beagle.log_prob_conditioned_branch_only(
                fasta_file,
                subst_model,
                cast(starting_values["frequencies"]),
                kappa=cast(starting_values["kappa"]),
                rescaling=vi_config["rescaling"],
                newick_file=newick_file,
                dated=True,
            )
        else:
            raise ValueError(
                "Only fixed substitution model parameters supported with Beagle likelihood"
            )
    else:
        raise ValueError("Unknown substitution model: {0}".format(model.subst_model))


def get_approx_dict(clock_approx, tree, model):
    if clock_approx in [
        "scaled_conjugate",
        "scaled_shrinkage_conjugate",
        "scaled_shrinkage_local_conjugate",
    ]:
        return dict(
            rate_loc=dict(
                approx="normal_conjugate",
                param="loc",
                prior_params=model.clock_params["rate_loc_precision"]["normalgamma"],
                input_callable=lambda f: (
                    lambda rates, rate_precision: f(tf.math.log(rates), rate_precision)
                ),
            ),
            rate_precision=dict(
                approx="normal_conjugate",
                param="precision",
                prior_params=model.clock_params["rate_loc_precision"]["normalgamma"],
                input_callable=lambda f: (lambda rates: f(tf.math.log(rates))),
            ),
        )
    else:
        return {}


def get_rate_approx_model(clock_approx):  # "mean_field",
    if clock_approx in ["scaled", "scaled_conjugate"]:
        return "scaled"
    elif clock_approx in ["scaled_shrinkage", "scaled_shrinkage_conjugate"]:
        return "scaled_shrinkage"
    elif clock_approx in ["scaled_shrinkage_local", "scaled_shrinkage_local_conjugate"]:
        return "scaled_shrinkage_local"
    else:
        return "mean_field"


def construct_approx(newick_file, model, clock_approx):
    tree, taxon_names = treeflow.tree_processing.parse_newick(newick_file)
    init_heights = tree["heights"]
    prior = get_phylo_prior(init_heights[: (init_heights.shape[0] + 1) // 2], model)
    prior_sample = prior.sample()
    q_dict, _ = treeflow.model.construct_prior_approximation(
        prior, prior_sample, approxs=get_approx_dict(clock_approx, tree, model)
    )
    q_tree, _ = treeflow.model.construct_tree_approximation(newick_file)
    q_dict["tree"] = q_tree

    if model.clock_model in RELAXED_CLOCK_MODELS:
        q_rates, _ = treeflow.model.construct_rate_approximation(
            treeflow.model.get_concrete_dist(prior.model["rates"], prior_sample),
            prior,
            approx_model=get_rate_approx_model(clock_approx),
        )
        q_dict["rates"] = q_rates
    return tfp.distributions.JointDistributionNamed(q_dict)


def get_variational_fit(
    newick_file,
    fasta_file,
    starting_values,
    model,
    vi_config,
    clock_approx,
    debug_log_dir=None,
    trace_fn=None,
    vi_kwargs={},
    function_mode=True,
):
    if debug_log_dir is not None:
        tf.debugging.experimental.enable_dump_debug_info(
            debug_log_dir,
            tensor_debug_mode="FULL_HEALTH",
            circular_buffer_size=-1,
        )

    likelihood, instance = get_likelihood(
        newick_file, fasta_file, starting_values, model, vi_config
    )
    tree_info = treeflow.libsbn.get_tree_info(instance)
    init_heights = tree_info.tree["heights"]
    prior = get_phylo_prior(init_heights[: (init_heights.shape[0] + 1) // 2], model)
    prior_sample = prior.sample()
    log_p = treeflow.model.get_log_posterior(prior, likelihood)

    q_dict, vars = treeflow.model.construct_prior_approximation(
        prior,
        prior_sample,
        init_mode={key: cast(value) for key, value in starting_values.items()},
        approxs=get_approx_dict(clock_approx, tree_info.tree, model),
    )

    q_tree, tree_vars = treeflow.model.construct_tree_approximation(
        newick_file, inst=instance
    )
    q_dict["tree"] = q_tree
    vars["tree"] = tree_vars

    if model.clock_model in RELAXED_CLOCK_MODELS:

        q_rates, rate_vars = treeflow.model.construct_rate_approximation(
            treeflow.model.get_concrete_dist(prior.model["rates"], prior_sample),
            prior,
            approx_model=get_rate_approx_model(clock_approx),
        )
        q_dict["rates"] = q_rates
        vars["rates"] = rate_vars

    q = tfp.distributions.JointDistributionNamed(q_dict)

    if trace_fn is None:  # TODO: Tidy this up
        loss, params = fit_surrogate_posterior(
            log_p, q, vi_config, vi_kwargs=vi_kwargs, function_mode=function_mode
        )
        return dict(loss=loss, vars=vars, params=params)
    else:
        return (
            fit_surrogate_posterior(
                log_p,
                q,
                vi_config,
                vi_kwargs=vi_kwargs,
                trace_fn=trace_fn,
                function_mode=function_mode,
            ),
            vars,
        )


def reconstruct_approx(
    newick_file, variational_fit, model, clock_approx
):  # TODO: Use libsbn parsing?
    vars = variational_fit["vars"]
    tree, taxon_names = treeflow.tree_processing.parse_newick(newick_file)
    init_heights = tree["heights"]
    prior = get_phylo_prior(init_heights[: (init_heights.shape[0] + 1) // 2], model)
    prior_sample = prior.sample()
    q_dict, _ = treeflow.model.construct_prior_approximation(
        prior,
        prior_sample,
        vars=vars,
        approxs=get_approx_dict(clock_approx, tree, model),
    )
    q_tree, _ = treeflow.model.construct_tree_approximation(
        newick_file, vars=vars["tree"]
    )
    q_dict["tree"] = q_tree

    if model.clock_model in RELAXED_CLOCK_MODELS:
        q_rates, _ = treeflow.model.construct_rate_approximation(
            treeflow.model.get_concrete_dist(prior.model["rates"], prior_sample),
            prior,
            approx_model=get_rate_approx_model(clock_approx),
            vars=vars["rates"],
        )
        q_dict["rates"] = q_rates

    return tfp.distributions.JointDistributionNamed(q_dict)


def build_init_value_string(init_value) -> str:
    if isinstance(init_value, list):
        return "|".join([str(x) for x in init_value])
    else:
        return str(init_value)


def build_init_values_string(init_values_dict):
    return ",".join(
        [
            f"{key}={build_init_value_string(value)}"
            for key, value in init_values_dict.items()
        ]
    )


def compute_trace_log_likelihoods(
    trace: VIResults,
    model: Model,
    tree_file,
    alignment_file,
):
    tree = convert_tree_to_tensor(parse_newick(tree_file))
    alignment = Alignment(alignment_file).get_compressed_alignment()
    encoded_sequences = alignment.get_encoded_sequence_tensor(tree.taxon_set)
    pattern_counts = alignment.get_weights_tensor()
    dist = phylo_model_to_joint_distribution(model, tree, alignment, pattern_counts)
    pinned_model = dist.experimental_pin(alignment=encoded_sequences)
    base_dict = trace.parameters
    approximation, variables_dict = get_fixed_topology_mean_field_approximation(
        model, topology_pins=dict(tree=tree.topology)
    )


def get_model_name(model):
    if isinstance(model, dict):
        return next(iter(model.keys()))
    else:
        return model


def generate_models_grid(models_grid_dict):

    datasets = models_grid_dict["datasets"]
    raw_tree_models = models_grid_dict["models"]["tree"]
    clock_models = models_grid_dict["models"]["clock"]
    subst_models = models_grid_dict["models"]["substitution"]
    site_models = models_grid_dict["models"]["site"]

    dated_tree_models = []
    undated_tree_models = []
    for tree_model_dict in raw_tree_models:
        tree_model_name, inner_tree_model_dict = next(iter(tree_model_dict.items()))
        model_dated = inner_tree_model_dict.pop("dates_possible")
        if model_dated:
            dated_tree_models.append(tree_model_dict)
        else:
            undated_tree_models.append(tree_model_dict)

    models_output = dict()
    for dataset_name, dataset_dict in datasets.items():
        dataset_dated = dataset_dict["dates"]
        if dataset_dated:
            tree_models = dated_tree_models
        else:
            tree_models = dated_tree_models + undated_tree_models
        for tree_model in tree_models:
            tree_model_name = get_model_name(tree_model)
            for clock_model in clock_models:
                clock_model_name = get_model_name(clock_model)
                if clock_model_name == "strictfixed":
                    final_clock_model = dict(strict=clock_model["strictfixed"])
                else:
                    final_clock_model = clock_model
                for subst_model in subst_models:
                    subst_model_name = get_model_name(subst_model)
                    for site_model in site_models:
                        site_model_name = get_model_name(site_model)
                        if site_model_name == "none":
                            site_model_name = "nosite"
                        grid_dataset_name = f"{dataset_name}-{tree_model_name}-{clock_model_name}-{subst_model_name}-{site_model_name}"
                        models_output[grid_dataset_name] = dict(
                            dataset_dict,
                            model=dict(
                                tree=tree_model,
                                clock=final_clock_model,
                                substitution=subst_model,
                                site=site_model,
                            ),
                        )
    return models_output
