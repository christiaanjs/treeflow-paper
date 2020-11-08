import pytest
import treeflow
import tensorflow as tf
import treeflow_pipeline.model as mod
import tensorflow_probability as tfp

@pytest.mark.parametrize("clock_model", ["relaxed_lognormal", "strict"])
def test_get_phylo_prior(clock_model):
    default_prior = dict(lognormal=dict(loc=0.0, scale=1.0))
    model_dict = dict(
        tree=dict(coalescent_constant=dict(pop_size=default_prior)),
        clock={ clock_model: dict(clock_rate=default_prior) },
        substitution=dict(hky=dict(kappa="fixed", frequencies="fixed")),
        site="none"
    )
    if clock_model == "relaxed_lognormal":
        model_dict["clock"][clock_model]["rate_sd"] = default_prior
    sampling_times = tf.zeros(3, dtype=treeflow.DEFAULT_FLOAT_DTYPE_TF)
    model = mod.Model(model_dict)
    dist = mod.get_phylo_prior(sampling_times, model)
    dist_sample = dist.sample()
    for param in ["pop_size", "clock_rate"]:
        assert param in dist_sample
    if clock_model == "relaxed_lognormal":
        assert "rate_sd" in dist_sample
        assert "rates" in dist_sample
