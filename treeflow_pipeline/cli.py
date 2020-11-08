import click
import pathlib
import importlib.resources
import snakemake
import yaml
import treeflow_pipeline.model

DEFAULT_SEED = 123

class Context:
    def __init__(self, alignment_path, model_file, output_path, seed):
        self.alignment_path = alignment_path
        self.model = treeflow_pipeline.model.Model(yaml.load(model_file))
        self.output_path = pathlib.Path(output_path)
        if seed is None:
            self.seed = DEFAULT_SEED # TODO: What's the best approach when we don't have a seed?

@click.group()
@click.argument("alignment", type=click.Path(exists=True), required=True)
@click.argument("model", type=click.File(), required=True)
@click.argument("output", type=click.Path(), required=True)
@click.option("-s", "--seed", type=int)
@click.pass_context
def cli(ctx, alignment, model, output, seed):
    ctx.obj = Context(alignment, model, output, seed)

@cli.command()
@click.option("--tree-method", type=click.Choice(["raxml"], case_sensitive=False), default="raxml")
@click.option("--rooting-method", type=click.Choice(["lsd"], case_sensitive=False), default="lsd")
@click.option("-w", "--working-directory", type=click.Path())
@click.pass_obj
def infer_topology(ctx, working_directory, tree_method, rooting_method):
    if working_directory is None:
        working_directory = ctx.output_path.parents[0]
    
    with importlib.resources.path("treeflow_pipeline", "topology.smk") as snakefile:
        snakemake.snakemake(
            snakefile,
            config=dict(
                alignment=ctx.alignment_path,
                output=ctx.output_path,
                working_directory=working_directory,
                tree_method=tree_method,
                rooting_method=rooting_method,
                subst_model=ctx.model.subst_model,
                site_model=ctx.model.site_model,
                seed=ctx.seed
            ),
            targets=["tree", "starting_values"]
        )

@cli.command() # TODO: How can we automatically infer topology if needed?
@click.option("-t", "--topology", type=click.Path(exists=True))
@click.option("-s", "--starting-values", type=click.File())
@click.pass_obj
def variational_fit(ctx, topology, alignment, output):
    if topology is None:
        print("There is no topology")
    else:
        print("The topology is " + topology)
    print("The alignment is " + alignment)