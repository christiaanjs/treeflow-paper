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
DEFAULT_ROOTING_METHOD = "lsd-dates"
DEFAULT_LSD_OUTPUT_FORMAT = "newick"


@cli.command()
@click.option(
    "--tree-method",
    type=click.Choice(["raxml"], case_sensitive=False),
    default=DEFAULT_TREE_METHOD,
)
@click.option(
    "--rooting-method",
    type=click.Choice(["lsd", "lsd-dates"], case_sensitive=False),
    default=DEFAULT_ROOTING_METHOD,
)
@click.option("-w", "--working-directory", type=click.Path())
@click.option(
    "--lsd-output-format",
    type=click.Choice(["newick", "nexus"]),
    default=DEFAULT_LSD_OUTPUT_FORMAT,
)
@click.pass_obj
def infer_topology(
    obj, working_directory, tree_method, rooting_method, lsd_output_format
):
    if working_directory is None:
        working_directory = obj.output_path.parents[0]

    with importlib.resources.path("treeflow_pipeline", "topology.smk") as snakefile:
        success = snakemake.snakemake(
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
                tree_model=obj.model.tree_model,
                lsd_output_format=lsd_output_format,
                seed=obj.seed,
            ),
            targets=["tree", "starting_values"],
            lock=False,
            forceall=True,
        )
    if not success:
        raise click.UsageError(
            "Topology inference pipeline was unsuccessful, check inputs"
        )
