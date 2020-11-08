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
        return next(iter(model.items())) # TODO: Parsing checks
    else:
        return model, None

class Model:
    def __init__(self, dict):
        self.tree_model, self.tree_params = parse_model(dict["tree"])
        self.clock_model, self.clock_params = parse_model(dict["clock"])
        self.subst_model, self.subst_params = parse_model(dict["substitution"])
        self.site_model, self.site_params = parse_model(dict["site"])

def cast(x):
    return tf.convert_to_tensor(x, treeflow.DEFAULT_FLOAT_DTYPE_TF)

def get_lognormal(loc, scale):
    return tfd.LogNormal(loc=cast(loc), scale=cast(scale))

def get_phylo_prior(sampling_times, prior_params, clock_model):
    sampling_times = cast(sampling_times)
    taxon_count = sampling_times.shape[-1]
    model_dict = {
        'pop_size': get_lognormal(prior_params['pop_size']['m'], prior_params['pop_size']['s']),
        'clock_rate': get_lognormal(prior_params['clock_rate']['m'], prior_params['clock_rate']['s']),
        'tree': lambda pop_size: treeflow.coalescent.ConstantCoalescent(taxon_count, pop_size, sampling_times)
    }
    if clock_model == 'relaxed':
        model_dict['rate_sd'] = get_lognormal(prior_params['rate_sd']['m'], prior_params['rate_sd']['s'])
        model_dict['rates'] = lambda rate_sd: tfd.Sample(tfd.LogNormal(loc=-rate_sd ** 2.0 / 2, scale=rate_sd), sample_shape=2*taxon_count - 2)
    elif clock_model != 'strict':
        raise ValueError('Clock model not yet implemented: {0}'.format(clock_model))
    return tfd.JointDistributionNamed(model_dict)

optimizers = dict(adam=tf.optimizers.Adam)

def fit_surrogate_posterior(log_p, q, vi_config):
    trace_fn = lambda x: (x.loss, x.parameters)
    return tfp.vi.fit_surrogate_posterior(
        log_p,
        q,
        optimizers[vi_config['optimizer']](**vi_config['optimizer_kwargs']),
        vi_config['num_steps'],
        trace_fn=trace_fn
    ) # TODO: Convergence criterion

def get_variational_fit(newick_file, fasta_file, starting_values, prior_params, vi_config, clock_model, clock_approx):
    subst_model = treeflow.substitution_model.HKY()
    likelihood, instance = treeflow.beagle.log_prob_conditioned_branch_only(
        fasta_file,
        subst_model,
        cast(starting_values['frequencies']),
        kappa=cast(starting_values['kappa']),
        rescaling=vi_config['rescaling'],
        newick_file=newick_file)
    tree_info = treeflow.libsbn.get_tree_info(instance)

    init_heights = tree_info.tree['heights']
    prior = get_phylo_prior(init_heights[:(init_heights.shape[0] + 1)//2], prior_params, clock_model)
    log_p = treeflow.model.get_log_posterior(prior, likelihood, relaxed_clock=clock_model == 'relaxed')

    q_dict, vars = treeflow.model.construct_prior_approximation(
        prior,
        init_mode=dict(
            pop_size=cast(starting_values['pop_size']),
            clock_rate=cast(starting_values['clock_rate'])
        )
    )
    
    q_tree, tree_vars = treeflow.model.construct_tree_approximation(newick_file, inst=instance)
    q_dict['tree'] = q_tree
    vars['tree'] = tree_vars

    if clock_model == 'relaxed':
        q_rates, rate_vars = treeflow.model.construct_rate_approximation(prior.model['rates'](cast(1.0)), approx_model=clock_approx)
        q_dict['rates'] = q_rates
        vars['rates'] = rate_vars

    q = tfp.distributions.JointDistributionNamed(q_dict)
    loss, params = fit_surrogate_posterior(log_p, q, vi_config)
    return dict(loss=loss, vars=vars, params=params)

def reconstruct_approx(newick_file, variational_fit, prior_params, clock_model, clock_approx):
    vars = variational_fit['vars']
    tree, taxon_names = treeflow.tree_processing.parse_newick(newick_file)
    init_heights = tree['heights']
    prior = get_phylo_prior(init_heights[:(init_heights.shape[0] + 1)//2], prior_params, clock_model)
    q_dict, _ = treeflow.model.construct_prior_approximation(prior, vars=vars)
    q_tree, _ = treeflow.model.construct_tree_approximation(newick_file, vars=vars['tree'])
    q_dict['tree'] = q_tree

    if clock_model == 'relaxed':
        q_rates, _ = treeflow.model.construct_rate_approximation(prior.model['rates'](cast(1.0)), approx_model=clock_approx, vars=vars['rates'])
        q_dict['rates'] = q_rates

    return tfp.distributions.JointDistributionNamed(q_dict)

    