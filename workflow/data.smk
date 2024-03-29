from treeflow_pipeline.util import yaml_input, yaml_output, text_input, text_output, sequence_input, pickle_input
from treeflow.model.phylo_model import PhyloModel
import treeflow_pipeline.templating as tem
from treeflow_pipeline.model import build_init_values_string
from treeflow_pipeline.results import (
    extract_trace_plot_data,
    compute_empirical_nucleotide_frequencies,
    assemble_timing_data
)
from treeflow_pipeline.data import convert_dates_to_numeric, extract_xml_sequences, remove_identical_sequences
import pathlib

configfile: "config/data-config.yaml"

wd = pathlib.Path(config["working_directory"])
all_models = yaml_input(config["model_file"])
datasets = list(all_models.keys())

models = { dataset: PhyloModel(all_models[dataset]["model"]) for dataset in datasets }

dataset_dir = "{dataset}"
data_dir = pathlib.Path("data")
default_out_dir = pathlib.Path("out")

rule data:
    input:
        expand(wd / dataset_dir / "marginals.png", dataset=datasets),
        expand(wd / dataset_dir / "traces.png", dataset=datasets),
        expand(wd / dataset_dir / "timing-data.csv", dataset=datasets),

rule carnivores_data_xml:
    output:
        wd / "carnivores-data.xml"
    shell:
        "curl {config[carnivores_xml_url]} -o {output}"

rule xml_to_fasta:
    input:
        wd / "{dataset}-data.xml"
    output:
        data_dir / "{dataset}.fasta"
    run:
        extract_xml_sequences(input[0], output[0], "fasta", reformat_taxon_name=True)

rule flu_fasta:
    input:
        data_dir / "h3n2.xml"
    output:
        data_dir / "h3n2.fasta"
    run:
        extract_xml_sequences(input[0], output[0], "fasta", date_from_trait_string=True)

rule big_flu_download:
    output:
        data_dir / "H3N2_HA_2011_2013.fasta"
    run:
        "curl {config[h3n2_fasta_url]} -o {output}"

rule big_flu_fasta:
    input:
        rules.big_flu_download.output[0]
    output:
        data_dir / "h3n2_big.fasta"
    run:
        convert_dates_to_numeric(
            input[0],
            "fasta",
            output[0],
            "fasta",
            date_format="%m/%d/%Y",
            split_char="|",
            split_index=-2,
            replace_old_date=False,
            exclude_bad_dates=True
        )

rule big_flu_reduction:
    input:
        data_dir / "h3n2_big.fasta"
    output:
        data_dir / "h3n2_big_reduced.fasta"
    run:
        remove_identical_sequences(
            input[0],
            "fasta",
            output[0],
            "fasta"
        )
        

rule carnivores_beast_run:
    input:
        xml = data_dir / "carnivores-beast2.xml"
    output:
        trees = default_out_dir / "carnivores-beast2.trees",
        trace = default_out_dir / "carnivores-beast2.log"
    shell:
        "beast -overwrite {input.xml}"

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
        wd = str(wd / dataset_dir),
        rooting_method = lambda wildcards: "lsd-dates" if all_models[wildcards.dataset]["dates"] else "lsd"
    output:
        topology = wd / dataset_dir / "topology.nwk",
        starting_values = wd / dataset_dir / "starting-values.yaml"
    shell:
        """
        treeflow_pipeline -s {config[seed]} \
            {input.fasta} {input.model_file} {output.topology} \
            infer-topology -w {params.wd} \
            --rooting-method {params.rooting_method} \
            --lsd-output-format {config[lsd_output_format]}
        """

rule beast_xml:
    input:
        fasta = lambda wildcards: all_models[wildcards.dataset]["alignment"],
        topology = wd / dataset_dir / "topology.nwk",
        starting_values = wd / dataset_dir / "starting-values.yaml",
        beast_config = "config/beast-config.yaml"
    output:
        wd / dataset_dir / "beast.xml"
    run:
        text_output(tem.build_beast_analysis(
            sequence_input(input.fasta),
            text_input(input.topology),
            yaml_input(input.starting_values),
            models[wildcards.dataset],
            yaml_input(input.beast_config),
            output[0],
            dated=all_models[wildcards.dataset]["dates"]
        ), output[0])

rule beast_run:
    input:
        wd / dataset_dir / "beast.xml"
    output:
        trace = wd / dataset_dir / "beast.log",
        trees = wd / dataset_dir / "beast.trees"
    benchmark:
        wd / dataset_dir / "beast-benchmark.txt"
    log:
        wd / dataset_dir / "beast-log.txt"
    shell:
        "beast -seed {config[seed]} {input} 2>&1 | tee {log}"

rule variational_fit:
    input:
        fasta = lambda wildcards: all_models[wildcards.dataset]["alignment"],
        topology = wd / dataset_dir / "topology.nwk",
        starting_values = wd / dataset_dir / "starting-values.yaml",
        model_file = wd / dataset_dir / "model.yaml"
    params:
        starting_values_string = lambda wildcards, input: build_init_values_string(yaml_input(input.starting_values))
    output:
        trace = wd / dataset_dir / "variational-trace.pickle",
        samples = wd / dataset_dir / "variational-samples.csv",
        tree_samples = wd / dataset_dir / "variational-tree-samples.nexus"
    benchmark:
        wd / dataset_dir / "variational-benchmark.txt"
    log:
        wd / dataset_dir / "variational-log.txt"
    shell:
        '''
        treeflow_vi -s {config[seed]} \
            -i {input.fasta} \
            -m {input.model_file} \
            -t {input.topology} \
            -n 30000 \
            --learning-rate 0.001 \
            --init-values "{params.starting_values_string}" \
            --trace-output {output.trace} \
            --samples-output {output.samples} \
            --tree-samples-output {output.tree_samples} \
            --n-output-samples {config[n_variational_samples]} \
            2>&1 | tee {log}
        '''


rule ml_fit:
    input:
        fasta = lambda wildcards: all_models[wildcards.dataset]["alignment"],
        topology = wd / dataset_dir / "topology.nwk",
        starting_values = wd / dataset_dir / "starting-values.yaml",
        model_file = wd / dataset_dir / "model.yaml"
    params:
        starting_values_string = lambda wildcards, input: build_init_values_string(yaml_input(input.starting_values))
    output:
        trace = wd / dataset_dir / "ml-trace.pickle",
        variables = wd / dataset_dir / "ml-variables.csv",
        tree = wd / dataset_dir / "ml-tree.nexus"
    shell:
        '''
        treeflow_ml \
            -i {input.fasta} \
            -m {input.model_file} \
            -t {input.topology} \
            -n 100 \
            --learning-rate 0.01 \
            --init-values "{params.starting_values_string}" \
            --trace-output {output.trace} \
            --variables-output {output.variables} \
            --tree-output {output.tree}
        '''

rule trace_plot_data:
    input:
        variational_trace = rules.variational_fit.output.trace,
        ml_trace = rules.ml_fit.output.trace
    output:
        wd / dataset_dir / "trace-plot-data.csv"
    run:
        extract_trace_plot_data(
            pickle_input(input.variational_trace),
            pickle_input(input.ml_trace),
            output[0]
        )

rule trace_plot:
    input:
        rules.trace_plot_data.output[0]
    output:
        wd / dataset_dir / "traces.png"
    script:
        "../scripts/trace-plot.R"

rule empirical_frequencies:
    input:
        lambda wildcards: all_models[wildcards.dataset]["alignment"]
    output:
        csv = wd / dataset_dir / "empirical-frequencies.csv"
    run:
        compute_empirical_nucleotide_frequencies(input[0], output.csv)


rule marginals_plot:
    input:
        vi_samples = rules.variational_fit.output.samples,
        beast_samples = rules.beast_run.output.trace,
        ml_variables = rules.ml_fit.output.variables,
        empirical_frequencies = rules.empirical_frequencies.output.csv
    output:
        wd / dataset_dir / "marginals.png"
    script:
        "../scripts/data-marginals-plot.R"

rule timing_data:
    input:
        vi_benchmark = wd / dataset_dir / "variational-benchmark.txt",
        vi_trace = wd / dataset_dir / "variational-trace.pickle",
        beast_benchmark = wd / dataset_dir / "beast-benchmark.txt",
        beast_trace = wd / dataset_dir / "beast.log"
    output:
        wd / dataset_dir / "timing-data.csv"
    run:
        assemble_timing_data(
            input.vi_benchmark,
            input.vi_trace,
            input.beast_benchmark,
            input.beast_trace,
            output[0]
        )
        