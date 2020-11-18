from treeflow_pipeline.util import yaml_input, yaml_output
import treeflow_pipeline.model
import pathlib

configfile: "config/data-config.yaml"

wd = pathlib.Path(config["working_directory"])
all_models = yaml_input(config["model_file"])
datasets = list(all_models.keys())

models = { dataset: treeflow_pipeline.model.Model(all_models[dataset]["model"]) for dataset in datasets }

dataset_dir = "{dataset}"

rule data:
    input:
        expand(wd / dataset_dir / "variational-fit-scaled.pickle", dataset=datasets)

rule model_files:
    output:
        wd / dataset_dir / "model.yaml"
    run:
        yaml_output(all_models[wildcards.dataset]["model"], output[0])

rule topology:
    input:
        fasta = lambda wildcards: all_models[wildcards.dataset]["alignment"],
        model_file = wd / dataset_dir / "model.yaml"
    params:
        wd = str(wd / dataset_dir)
    output:
        topology = wd / dataset_dir / "topology.nwk",
        starting_values = wd / dataset_dir / "starting-values.yaml"
    shell:
        """
        treeflow_pipeline -s {config[seed]} \
            {input.fasta} {input.model_file} {output.topology} \
            infer-topology -w {params.wd}
        """ 

rule variational_fit:
    input:
        fasta =lambda wildcards: all_models[wildcards.dataset]["alignment"],
        topology = wd / dataset_dir / "topology.nwk",
        starting_values = wd / dataset_dir / "starting-values.yaml",
        model_file = wd / dataset_dir / "model.yaml"
    output:
        wd / dataset_dir / "variational-fit-{clock_approx}.pickle"
    shell:
        """
        treeflow_pipeline -s {config[seed]} \
            {input.fasta} {input.model_file} {output} \
            variational-fit \
            -t {input.topology} \
            -s {input.starting_values} \
            -c {wildcards.clock_approx} \
            --config {config[vi_config]}
        """