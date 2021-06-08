import pytest
import treeflow
import tensorflow as tf
import treeflow_pipeline.model as mod
import tensorflow_probability as tfp


default_prior = dict(lognormal=dict(loc=0.0, scale=1.0))
clock_model_dicts_and_keys = [
    (dict(strict=dict(clock_rate=default_prior)), ["clock_rate"]),
    (
        dict(relaxed_lognormal=dict(clock_rate=default_prior, rate_sd=default_prior)),
        ["clock_rate", "rate_sd", "rates"],
    ),
    (
        dict(
            relaxed_lognormal_conjugate=dict(
                loc_precision=dict(
                    normalgamma=dict(
                        loc=-8.117321296021004,
                        precision_scale=0.6145984625809157,
                        concentration=2.0264094157104164,
                        rate=0.055894814133092476,
                    )
                )
            )
        ),
        ["rate_loc", "rate_precision", "rates"],
    ),
]


@pytest.mark.parametrize(
    "clock_model_dict,clock_model_keys", clock_model_dicts_and_keys
)
def test_get_phylo_prior(clock_model_dict, clock_model_keys):
    model_dict = dict(
        tree=dict(coalescent_constant=dict(pop_size=default_prior)),
        clock=clock_model_dict,
        substitution=dict(hky=dict(kappa="fixed", frequencies="fixed")),
        site="none",
    )
    sampling_times = tf.zeros(3, dtype=treeflow.DEFAULT_FLOAT_DTYPE_TF)
    model = mod.Model(model_dict)
    dist = mod.get_phylo_prior(sampling_times, model)
    dist_sample = dist.sample()
    for param in ["pop_size"] + clock_model_keys:
        assert param in dist_sample
