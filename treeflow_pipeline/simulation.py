import numpy as np
import treeflow_pipeline.simulation

def rng(sim_config):
    return  np.random.default_rng(sim_config['seed'])

def get_sampling_times(sim_config):
    times = rng(sim_config).uniform(0.0, sim_config['sampling_window'], size=sim_config['n_taxa'])
    return { 'T{0}_{1}'.format(i + 1, time): float(time) for i, time in enumerate(times) }

def sample_lognormal(sim_config, prior_params, param):
    return rng(sim_config).lognormal(prior_params[param]['m'], prior_params[param]['s'])

get_pop_size, get_rate_sd, get_clock_rate = [
    lambda x, y: sample_lognormal(x, y, p) for p in ['pop_size', 'rate_sd', 'clock_rate']
]