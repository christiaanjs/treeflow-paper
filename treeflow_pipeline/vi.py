import typing as tp
import pickle
from functools import partial
import tensorflow as tf
from tensorflow_probability.python.util import DeferredTensor
from tensorflow_probability.python.bijectors.softplus import _stable_grad_softplus
from tensorflow_probability.python.distributions import (
    Distribution,
    Independent,
    Normal,
    JointDistributionNamed,
    TransformedDistribution,
)
from tensorflow_probability.python.vi import fit_surrogate_posterior
from tensorflow_probability.python.math.minimize import minimize
from treeflow import DEFAULT_FLOAT_DTYPE_TF
from treeflow.tree.rooted.tensorflow_rooted_tree import convert_tree_to_tensor
from treeflow.tree.io import parse_newick
from treeflow.evolution.seqio import Alignment
from treeflow.model.io import write_samples_to_file
from treeflow.cli.inference_common import (
    optimizer_builders,
    ROBUST_ADAM_KEY,
    parse_init_values,
    EXAMPLE_PHYLO_MODEL_DICT,
    get_tree_vars,
    write_trees,
)
from treeflow.cli.vi import convergence_criterion_classes
import yaml
from treeflow.model.phylo_model import (
    phylo_model_to_joint_distribution,
    PhyloModel,
    DEFAULT_TREE_VAR_NAME,
)
from treeflow.model.event_shape_bijector import (
    get_fixed_topology_event_shape_and_space_bijector,
    get_fixed_topology_event_shape,
    get_unconstrained_init_values,
)
from treeflow.vi.fixed_topology_advi import default_vi_trace_fn


def build_base_approx_for_var(
    name, event_shape, init_loc: tp.Optional[tf.Tensor]
) -> tp.Tuple[Distribution, tp.Dict[str, tf.Variable]]:
    loc_name = f"{name}_loc"
    loc = tf.Variable(
        tf.zeros(event_shape, dtype=DEFAULT_FLOAT_DTYPE_TF)
        if (init_loc is None)
        else init_loc,
        name=loc_name,
    )
    scale_name = f"{name}_scale"
    scale = tf.Variable(
        tf.zeros(event_shape, dtype=DEFAULT_FLOAT_DTYPE_TF), name=scale_name
    )
    return Independent(
        Normal(loc, DeferredTensor(scale, _stable_grad_softplus)),
        reinterpreted_batch_ndims=tf.shape(event_shape)[0],
    ), {loc.name: loc, scale.name: scale}


def fit_vi_alternate(
    input,
    topology,
    num_steps,
    model_file,
    seed=None,
    init_values=None,
    trace_output=None,
    samples_output=None,
    tree_samples_output=None,
    learning_rate=1e-3,
    n_output_samples=200,
    optimizer="robust_adam",
    convergence_criterion="nonfinite",
):
    optimizer = optimizer_builders[optimizer](learning_rate=learning_rate)

    print(f"Parsing topology {topology}")
    tree = convert_tree_to_tensor(parse_newick(topology))

    print(f"Parsing alignment {input}")
    alignment = Alignment(input).get_compressed_alignment()
    encoded_sequences = alignment.get_encoded_sequence_tensor(tree.taxon_set)
    pattern_counts = alignment.get_weights_tensor()
    if model_file is None:
        model_dict = EXAMPLE_PHYLO_MODEL_DICT
    else:
        with open(model_file) as f:
            model_dict = yaml.safe_load(f)
    phylo_model = PhyloModel(model_dict)
    model = phylo_model_to_joint_distribution(
        phylo_model, tree, alignment, pattern_counts=pattern_counts
    )
    pinned_model = model.experimental_pin(alignment=encoded_sequences)

    topologies = {DEFAULT_TREE_VAR_NAME: tree.topology}
    bijector, base_event_shape = get_fixed_topology_event_shape_and_space_bijector(
        pinned_model, topology_pins=topologies
    )

    init_values_dict = (
        None
        if init_values is None
        else {
            key: tf.constant(value, dtype=DEFAULT_FLOAT_DTYPE_TF)
            for key, value in parse_init_values(init_values).items()
        }
    )
    model_names = set(pinned_model._flat_resolve_names())
    if init_values_dict is None:
        init_loc = None
    else:
        init_loc = {
            key: value for key, value in init_values_dict.items() if key in model_names
        }
        init_loc[DEFAULT_TREE_VAR_NAME] = tree

    init_unconstrained = get_unconstrained_init_values(
        model,
        bijector,
        event_shape_fn=partial(
            get_fixed_topology_event_shape, topology_pins=topologies
        ),
        init=init_loc,
    )
    dists_and_variables = {
        name: build_base_approx_for_var(name, var_event_shape, init_unconstrained[name])
        for name, var_event_shape in base_event_shape.items()
    }
    base_dist = JointDistributionNamed(
        {name: dist for name, (dist, _) in dists_and_variables.items()}
    )
    variables_dict = {
        var_name: var
        for name, (_, var_dict) in dists_and_variables.items()
        for var_name, var in var_dict.items()
    }
    trace_fn = partial(default_vi_trace_fn, variables_dict=variables_dict)
    if convergence_criterion is not None:
        convergence_criterion_instance = convergence_criterion_classes[
            convergence_criterion
        ]()
    else:
        convergence_criterion_instance = None

    approx = TransformedDistribution(base_dist, bijector)

    def loss():
        approx_sample = approx.sample()
        return approx.log_prob(approx_sample) - pinned_model.unnormalized_log_prob(
            approx_sample
        )

    print(f"Running VI for {num_steps} iterations...")
    # trace = fit_surrogate_posterior(
    #     pinned_model.unnormalized_log_prob,
    #     approx,
    #     optimizer,
    #     num_steps,
    #     trace_fn=trace_fn,
    # )
    trace = minimize(loss, num_steps, optimizer, trace_fn=trace_fn)
    print("Inference complete")

    inference_steps = trace.loss.shape[0]
    print(f"Ran inference for {inference_steps} iterations")

    if trace_output is not None:
        print(f"Saving trace to {trace_output}...")
        with open(trace_output, "wb") as f:
            pickle.dump(trace, f)

    if samples_output is not None or tree_samples_output is not None:
        print("Sampling fitted approximation...")
        samples = approx.sample(n_output_samples)
        samples_dict = samples._asdict()
        tree_vars = get_tree_vars(phylo_model)

        tree_samples = dict()
        for var in tree_vars:
            tree_samples[var] = samples_dict.pop(var)

        if samples_output is not None:
            print(f"Saving samples to {samples_output}...")
            write_samples_to_file(
                samples,
                pinned_model,
                samples_output,
                vars=samples_dict.keys(),
                tree_vars={DEFAULT_TREE_VAR_NAME: tree_samples[DEFAULT_TREE_VAR_NAME]},
            )

        if tree_samples_output is not None:
            print(f"Saving tree samples to {tree_samples_output}...")
            write_trees(tree_samples, topology, tree_samples_output)

    print("Exiting...")
