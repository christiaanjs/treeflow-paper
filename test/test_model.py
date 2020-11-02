import pytest
import treeflow
import tensorflow as tf
import treeflow_pipeline.model as mod
import tensorflow_probability as tfp

@pytest.mark.parametrize("clock_model", ["relaxed", "strict"])
def test_get_phylo_prior(clock_model):
    default_params = dict(m=0.0, s=1.0)
    prior_params = dict(
        pop_size=default_params,
        clock_rate=default_params
    )
    if clock_model == "relaxed":
        prior_params["rate_sd"] = default_params
    sampling_times = tf.zeros(3, dtype=treeflow.DEFAULT_FLOAT_DTYPE_TF)
    dist = mod.get_phylo_prior(sampling_times, prior_params, clock_model)
    dist_sample = dist.sample()
    for key in prior_params:
        assert key in dist_sample
    if clock_model == "relaxed":
        assert "rates" in dist_sample
