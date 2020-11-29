import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions
import treeflow
import treeflow.coalescent
import treeflow.model
import treeflow.substitution_model
import treeflow.sequences
import treeflow.beagle
import treeflow.libsbn

def parse_model(model):
    if isinstance(model, dict):
        return next(iter(model.items()))
    else:
        return model, None

RELAXED_CLOCK_MODELS = ["relaxed_lognormal"]
APPROX_MODELS = ["mean_field", "scaled"]#, "scaled_all"]

class Model:
    def __init__(self, dict):
        self.tree_model, self.tree_params = parse_model(dict["tree"])
        self.clock_model, self.clock_params = parse_model(dict["clock"])
        self.subst_model, self.subst_params = parse_model(dict["substitution"])
        self.site_model, self.site_params = parse_model(dict["site"])

    def all_params(self):
        return { key: value for comp_params in [self.tree_params, self.clock_params, self.subst_params, self.site_params] if comp_params is not None for key, value in comp_params.items() }

    def free_params(self):
        return { key: value for key, value in self.all_params().items() if value != 'fixed' }

    def relaxed_clock(self):
        return self.clock_model in RELAXED_CLOCK_MODELS

def get_non_rate_defaults(clock_model):
    if clock_model == "relaxed_lognormal":
        return dict(rate_sd=1.0)
    elif clock_model == "strict":
        return {}
    else:
        raise ValueError(f"Clock model not known: {clock_model}")

def cast(x):
    return tf.convert_to_tensor(x, treeflow.DEFAULT_FLOAT_DTYPE_TF)

dists = dict(
    lognormal=tfd.LogNormal
)

def get_dist(dist_dict):
    dist_name, params = parse_model(dist_dict) # TODO: Handle case of fixed values
    try:
        return dists[dist_name](**{ key: cast(val) for key, val in params.items() }) # TODO: Non-float params?
    except KeyError:
        raise ValueError("Distribution not known: {0}".format(dist_name))

def get_phylo_prior(sampling_times, model):
    sampling_times = cast(sampling_times)
    taxon_count = sampling_times.shape[-1]
    model_dict = {}

    if model.tree_model == "coalescent_constant":
        model_dict["pop_size"] = get_dist(model.tree_params["pop_size"]) # TODO: Make this more automatic
        model_dict["tree"] = lambda pop_size: treeflow.coalescent.ConstantCoalescent(taxon_count, pop_size, sampling_times)
    else:
        raise ValueError("Tree prior not known: {0}".format(model.tree_model))

    if model.clock_model == "relaxed_lognormal":
        model_dict["clock_rate"] = get_dist(model.clock_params["clock_rate"])
        model_dict["rate_sd"] = get_dist(model.clock_params["rate_sd"])
        model_dict["rates"] = lambda rate_sd: tfd.Sample(tfd.LogNormal(loc=-rate_sd ** 2.0 / 2, scale=rate_sd), sample_shape=2*taxon_count - 2)
    elif model.clock_model == "strict":
        model_dict["clock_rate"] = get_dist(model.clock_params["clock_rate"])
    else:
        raise ValueError("Clock model not known: {0}".format(clock_model))

    # TODO: Substitution model
    # TODO: Site model
    
    return tfd.JointDistributionNamed(model_dict)

optimizers = dict(adam=tf.optimizers.Adam)

def fit_surrogate_posterior(log_p, q, vi_config):
    trace_fn = lambda x: (x.loss, x.parameters)
    return tfp.vi.fit_surrogate_posterior(
        log_p,
        q,
        optimizers[vi_config["optimizer"]](**vi_config["optimizer_kwargs"]),
        vi_config["num_steps"],
        trace_fn=trace_fn,
        seed=vi_config["seed"]
    ) # TODO: Convergence criterion

def get_likelihood(newick_file, fasta_file, starting_values, model, vi_config):
    if model.subst_model == "hky":
        subst_model = treeflow.substitution_model.HKY()
        if model.subst_params["kappa"] == "fixed" and model.subst_params["frequencies"] == "fixed":
            return treeflow.beagle.log_prob_conditioned_branch_only(
                fasta_file,
                subst_model,
                cast(starting_values["frequencies"]),
                kappa=cast(starting_values["kappa"]),
                rescaling=vi_config["rescaling"],
                newick_file=newick_file,
                dated=True
            )
        else:
            raise ValueError("Only fixed substitution model parameters supported with Beagle likelihood")
    else:
        raise ValueError("Unknown substitution model: {0}".format(model.subst_model))

def get_approx_dict(clock_approx, tree):
    if clock_approx == "scaled_all":
        return dict(clock_rate=dict(approx="scaled", tree_statistic="length", tree=tree))
    else:
        return {}

def get_rate_approx_model(clock_approx):
    if clock_approx in ["scaled", "scaled_all"]:
        return "scaled"
    else:
        return "mean_field"

def get_variational_fit(newick_file, fasta_file, starting_values, model, vi_config, clock_approx):
    likelihood, instance = get_likelihood(newick_file, fasta_file, starting_values, model, vi_config)
    tree_info = treeflow.libsbn.get_tree_info(instance)
    init_heights = tree_info.tree["heights"]
    prior = get_phylo_prior(init_heights[:(init_heights.shape[0] + 1)//2], model)
    log_p = treeflow.model.get_log_posterior(prior, likelihood, relaxed_clock=model.clock_model in RELAXED_CLOCK_MODELS)

    q_dict, vars = treeflow.model.construct_prior_approximation(
        prior,
        init_mode=dict(
            pop_size=cast(starting_values["pop_size"]),
            clock_rate=cast(starting_values["clock_rate"])
        ),
        approxs=get_approx_dict(clock_approx, tree_info.tree)
    )
    
    q_tree, tree_vars = treeflow.model.construct_tree_approximation(newick_file, inst=instance)
    q_dict["tree"] = q_tree
    vars["tree"] = tree_vars

    if model.clock_model in RELAXED_CLOCK_MODELS:

        q_rates, rate_vars = treeflow.model.construct_rate_approximation(prior.model["rates"](cast(1.0)), approx_model=get_rate_approx_model(clock_approx))
        q_dict["rates"] = q_rates
        vars["rates"] = rate_vars

    q = tfp.distributions.JointDistributionNamed(q_dict)
    loss, params = fit_surrogate_posterior(log_p, q, vi_config)
    return dict(loss=loss, vars=vars, params=params)

def reconstruct_approx(newick_file, variational_fit, model, clock_approx):
    vars = variational_fit["vars"]
    tree, taxon_names = treeflow.tree_processing.parse_newick(newick_file)
    init_heights = tree["heights"]
    prior = get_phylo_prior(init_heights[:(init_heights.shape[0] + 1)//2], model)
    q_dict, _ = treeflow.model.construct_prior_approximation(prior, vars=vars, approxs=get_approx_dict(clock_approx, tree))
    q_tree, _ = treeflow.model.construct_tree_approximation(newick_file, vars=vars["tree"])
    q_dict["tree"] = q_tree

    if model.clock_model in RELAXED_CLOCK_MODELS:
        q_rates, _ = treeflow.model.construct_rate_approximation(prior.model["rates"](cast(1.0)), approx_model=get_rate_approx_model(clock_approx), vars=vars["rates"])
        q_dict["rates"] = q_rates

    return tfp.distributions.JointDistributionNamed(q_dict)

    