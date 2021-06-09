import pytest
import treeflow
import tensorflow as tf
import treeflow_pipeline.model as mod
import tensorflow_probability as tfp


def c(x):
    return tf.convert_to_tensor(x, dtype=treeflow.DEFAULT_FLOAT_DTYPE_TF)


default_prior = dict(lognormal=dict(loc=c(0.0), scale=c(1.0)))
strict_clock_dict = dict(strict=dict(clock_rate=default_prior))
strict_clock_keys = ["clock_rate"]
relaxed_clock_dict = dict(
    relaxed_lognormal=dict(clock_rate=default_prior, rate_sd=default_prior)
)
relaxed_clock_keys = ["clock_rate", "rate_sd", "rates"]
relaxed_conjugate_dict = dict(
    relaxed_lognormal_conjugate=dict(
        rate_loc_precision=dict(
            normalgamma=dict(
                loc=c(-8.117321296021004),
                precision_scale=c(0.6145984625809157),
                concentration=c(2.0264094157104164),
                rate=c(0.055894814133092476),
            )
        )
    )
)
relaxed_conjugate_keys = ["rate_loc", "rate_precision", "rates"]
clock_model_dicts_and_keys = [
    (strict_clock_dict, strict_clock_keys),
    (relaxed_clock_dict, relaxed_clock_keys),
    (relaxed_conjugate_dict, relaxed_conjugate_keys),
]


def clock_model_dict_to_model(clock_model_dict):
    model_dict = dict(
        tree=dict(coalescent_constant=dict(pop_size=default_prior)),
        clock=clock_model_dict,
        substitution=dict(hky=dict(kappa="fixed", frequencies="fixed")),
        site="none",
    )

    model = mod.Model(model_dict)
    return model


@pytest.fixture(params=clock_model_dicts_and_keys)
def model_and_keys(request):
    clock_model_dict, clock_model_keys = request.param

    return clock_model_dict_to_model(clock_model_dict), clock_model_keys + ["pop_size"]


def test_get_phylo_prior(model_and_keys):
    model, param_keys = model_and_keys
    sampling_times = tf.zeros(3, dtype=treeflow.DEFAULT_FLOAT_DTYPE_TF)
    dist = mod.get_phylo_prior(sampling_times, model)
    dist_sample = dist.sample()
    for param in param_keys:
        assert param in dist_sample
    assert "tree" in dist_sample


def test_construct_approx_mean_field(model_and_keys, test_newick_file):
    model, param_keys = model_and_keys
    dist = mod.construct_approx(str(test_newick_file), model, "mean_field")
    dist_sample = dist.sample()
    dist.log_prob(dist_sample)
    assert set(dist_sample.keys()) == set(param_keys + ["tree"])


relaxed_clock_model_dicts_and_keys = [
    (relaxed_clock_dict, relaxed_clock_keys),
    (relaxed_conjugate_dict, relaxed_conjugate_keys),
]


@pytest.fixture(params=relaxed_clock_model_dicts_and_keys)
def relaxed_clock_model_and_keys(request):
    clock_model_dict, clock_model_keys = request.param

    return clock_model_dict_to_model(clock_model_dict), clock_model_keys + ["pop_size"]


def test_construct_approx_scaled(relaxed_clock_model_and_keys, test_newick_file):
    model, param_keys = relaxed_clock_model_and_keys
    dist = mod.construct_approx(str(test_newick_file), model, "scaled")
    dist_sample = dist.sample()
    dist.log_prob(dist_sample)
    assert set(dist_sample.keys()) == set(param_keys + ["tree"])


@pytest.fixture
def conjugate_model_and_keys():
    return clock_model_dict_to_model(relaxed_conjugate_dict), relaxed_conjugate_keys + [
        "pop_size"
    ]


def test_construct_approx_conjugate(conjugate_model_and_keys, test_newick_file):
    model, param_keys = conjugate_model_and_keys

    dist = mod.construct_approx(str(test_newick_file), model, "scaled_conjugate")
    dist_sample = dist.sample()
    dist.log_prob(dist_sample)
    assert set(dist_sample.keys()) == set(param_keys + ["tree"])


@pytest.fixture()
def conjugate_model(conjugate_model_and_keys):
    return conjugate_model_and_keys[0]


@pytest.mark.parametrize("approx", ["mean_field", "scaled", "scaled_conjugate"])
def test_variational_fit_conjugate(
    conjugate_model, test_newick_file, test_fasta_file, approx
):

    sampling_times = tf.zeros(3, dtype=treeflow.DEFAULT_FLOAT_DTYPE_TF)
    prior = mod.get_phylo_prior(sampling_times, conjugate_model)
    prior_sample = prior.sample()
    prior_sample.pop("rates")
    prior_sample.pop("tree")
    starting_values = {
        "frequencies": [0.26, 0.24, 0.27, 0.23],
        "kappa": 2.0,
        **prior_sample,
    }
    vi_config = dict(
        optimizer="adam",
        optimizer_kwargs=dict(learning_rate=0.005),
        num_steps=10,
        rescaling=False,
        seed=3,
    )

    res = mod.get_variational_fit(
        str(test_newick_file),
        str(test_fasta_file),
        starting_values,
        conjugate_model,
        vi_config,
        approx,
    )


def test_reconstruct_approx_scaled_conjugate(
    conjugate_model_and_keys, test_newick_file, test_data_dir
):
    import pickle

    with open(test_data_dir / "variational-fit-scaled_conjugate.pickle", "rb") as f:
        variational_fit = pickle.load(f)
    model, param_keys = conjugate_model_and_keys
    approx = mod.reconstruct_approx(
        str(test_newick_file), variational_fit, model, "scaled_conjugate"
    )
    approx_sample = approx.sample(10, seed=4)
    print(approx_sample)
    assert set(approx_sample.keys()) == set(param_keys + ["tree"])
