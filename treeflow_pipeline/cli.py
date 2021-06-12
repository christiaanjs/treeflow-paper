import click
import click_config_file
import pathlib
import importlib.resources
import snakemake
import yaml
import treeflow_pipeline.model
import pickle

# TODO: Error handling - click.Abort on Exception

DEFAULT_SEED = 123


class Context:
    def __init__(self, alignment_path, model_file, output_path, seed):
        self.alignment_path = alignment_path
        self.model = treeflow_pipeline.model.Model(yaml.safe_load(model_file))
        self.output_path = pathlib.Path(output_path)
        if seed is None:
            self.seed = DEFAULT_SEED  # TODO: What's the best approach when we don't have a seed?
        else:
            self.seed = seed


@click.group()
@click.argument("alignment", type=click.Path(exists=True), required=False)
@click.argument("model", type=click.File(), required=False)
@click.argument("output", type=click.Path(), required=False)
@click.option("-s", "--seed", type=int)
@click.pass_context
def cli(ctx, alignment, model, output, seed):
    ctx.obj = Context(alignment, model, output, seed)


DEFAULT_TREE_METHOD = "raxml"
DEFAULT_ROOTING_METHOD = "lsd"


@cli.command()
@click.option(
    "--tree-method",
    type=click.Choice(["raxml"], case_sensitive=False),
    default=DEFAULT_TREE_METHOD,
)
@click.option(
    "--rooting-method",
    type=click.Choice(["lsd"], case_sensitive=False),
    default=DEFAULT_ROOTING_METHOD,
)
@click.option("-w", "--working-directory", type=click.Path())
@click.pass_obj
def infer_topology(obj, working_directory, tree_method, rooting_method):
    if working_directory is None:
        working_directory = obj.output_path.parents[0]

    with importlib.resources.path("treeflow_pipeline", "topology.smk") as snakefile:
        snakemake.snakemake(
            snakefile,
            config=dict(
                alignment=obj.alignment_path,
                output=obj.output_path,
                working_directory=working_directory,
                tree_method=tree_method,
                rooting_method=rooting_method,
                subst_model=obj.model.subst_model,
                site_model=obj.model.site_model,
                clock_model=obj.model.clock_model,
                seed=obj.seed,
            ),
            targets=["tree", "starting_values"],
            lock=False,
        )


def infer_topology_not_provided(ctx):
    fit_output = ctx.obj.output_path
    topology = fit_output.parents[0] / "topology.nwk"
    ctx.obj.output_path = topology
    click.echo(
        "No topology supplied. Inferring with {0} and rooting with {1}".format(
            DEFAULT_TREE_METHOD, DEFAULT_ROOTING_METHOD
        )
    )
    ctx.invoke(infer_topology)
    ctx.obj.output_path = fit_output

    with open(fit_output.parents[0] / "starting-values.yaml") as f:
        ml_start_values = yaml.safe_load(f)

    return topology, ml_start_values


def get_vi_config(ctx, optimizer, learning_rate, num_steps, rescaling):
    return dict(
        seed=ctx.obj.seed,
        optimizer=optimizer,
        num_steps=num_steps,
        rescaling=rescaling,
        optimizer_kwargs=dict(learning_rate=learning_rate),
    )


def yaml_config_provider(file_path, cmd_name):
    with open(file_path) as f:
        return yaml.safe_load(f)[cmd_name]


# TODO: What's the best way to specify inference configuration?
# TODO: Options for likelihood computation
@cli.command()
@click.option("-t", "--topology", type=click.Path(exists=True))
@click.option("-s", "--starting-values", type=click.File())
@click.option("-o", "--optimizer", type=click.Choice(["adam", "sgd"]), default="adam")
@click.option("-l", "--learning-rate", type=float, default=1e-2)
@click.option("-n", "--num-steps", type=int, default=40000)
@click.option(
    "-a",
    "--clock-approx",
    type=click.Choice(treeflow_pipeline.model.APPROX_MODELS),
    default="scaled",
)
@click.option("-r", "--rescaling", type=bool, default=True)
@click_config_file.configuration_option(provider=yaml_config_provider)
@click.pass_context
def variational_fit(
    ctx,
    topology,
    starting_values,
    optimizer,
    learning_rate,
    num_steps,
    clock_approx,
    rescaling,
):
    if starting_values is None:  # TODO: Infer start values
        start_values_dict = {}
    else:
        start_values_dict = yaml.safe_load(starting_values)
    if topology is None:
        topology, ml_start_values = infer_topology_not_provided(ctx)
        ml_start_values.update(start_values_dict)
        start_values_dict = ml_start_values

    res = treeflow_pipeline.model.get_variational_fit(
        str(topology),
        str(ctx.obj.alignment_path),
        start_values_dict,
        ctx.obj.model,
        get_vi_config(ctx, optimizer, learning_rate, num_steps, rescaling),
        clock_approx,
    )

    with open(ctx.obj.output_path, "wb") as f:
        pickle.dump(res, f)
