import click
import pathlib
import importlib.resources
import snakemake

@click.group() # TODO: What should we share with context?
def cli():
    pass

@cli.command()
@click.option("--tree-method", type=click.Choice(["raxml"], case_sensitive=False), default="raxml")
@click.option("--rooting-method", type=click.Choice(["lsd"], case_sensitive=False), default="lsd")
@click.option("-w", "--working-directory", type=click.Path())
@click.argument("alignment", type=click.Path(exists=True), required=True)
@click.argument("output", type=click.Path(), required=True)
def infer_topology(alignment, output, working_directory, tree_method, rooting_method):
    if working_directory is None:
        working_directory = pathlib.Path(output).parents[0]
    
    with importlib.resources.path("treeflow_pipeline", "topology.smk") as snakefile:
        snakemake.snakemake(
            snakefile,
            config=dict(
                alignment=alignment,
                output=output,
                working_directory=working_directory,
                tree_method=tree_method,
                rooting_method=rooting_method
            ),
            targets=["tree", "starting_values"]
        )

@cli.command()
@click.option("-t", "--topology", type=click.Path(exists=True))
@click.argument("alignment", type=click.Path(exists=True), required=True)
@click.argument("output", type=click.Path(), required=True)
def variational_fit(topology, alignment, output):
    if topology is None:
        print("There is no topology")
    else:
        print("The topology is " + topology)
    print("The alignment is " + alignment)