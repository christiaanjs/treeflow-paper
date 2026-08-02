from treeflow_pipeline.util import yaml_input, yaml_output, text_input, text_output, sequence_input, pickle_input
from treeflow.model.phylo_model import PhyloModel
import treeflow_pipeline.templating as tem
from treeflow_pipeline.model import build_init_values_string
from treeflow_pipeline.results import (
    extract_trace_plot_data,
    compute_empirical_nucleotide_frequencies,
    assemble_timing_data,
    get_runtime_from_benchmark_file,
)
from treeflow_pipeline.data import convert_dates_to_numeric, extract_xml_sequences, remove_identical_sequences
import pathlib
import sys
import pandas as pd
import dendropy

import treeflow
treeflow_dir = pathlib.Path(treeflow.__file__).parents[1]

configfile: "config/data-config.yaml"

wd = pathlib.Path(config["working_directory"])
all_models = yaml_input(config["model_file"])
datasets = list(all_models.keys())

models = { dataset: PhyloModel(all_models[dataset]["model"]) for dataset in datasets }

dataset_dir = "{dataset}"
data_dir = pathlib.Path("data")
default_out_dir = pathlib.Path("out")

# Both the flu (H3N2) and carnivores base-model VI fits are run multiple times
# with different seeds -- like examples/h3n2-vi-multi-run.sh in the treeflow
# repo, and matching the settings the carnivores example notebook uses for its
# own multi-run posterior (num_steps=60000, 4 runs, learning rate 0.001) -- so
# the manuscript figures can show a pooled posterior estimate with an inter-run
# Monte Carlo error band. root_full_rank is the approximation used for the
# main-text figures; full_rank and mean_field are fit as well so the
# supplementary material can compare the three families (see
# approximation_comparison_plot in workflow/ms.smk, which uses run 1 of each).
MULTI_RUN_DATASETS = ["h3n2", "carnivores"]
MULTI_RUN_APPROXES = ["full_rank", "mean_field", "root_full_rank"]
MULTI_RUN_SEEDS = [1, 2, 3, 4]
MULTI_RUN_NUM_STEPS = 60000
# full_rank's trace records its whole D x D scale matrix (D ~= 990 free dimensions
# for the H3N2 model) at every step by default, which is infeasible at 60,000
# steps (empirically ~1.5TB). --max-trace-coords bounds the trace to a fixed
# number of sampled coordinates per variable instead, which keeps memory to a
# few GB (empirically verified: ~5GB extrapolated at 60,000 steps for H3N2)
# without changing the optimisation itself. Applied to mean_field too for
# consistency, though its trace is much smaller regardless.
MULTI_RUN_MAX_TRACE_COORDS = 50

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
        # -beagle_CPU pins BEAGLE to the CPU resource. BEAGLE also exposes the
        # GPU via OpenCL on this machine, and letting it choose would make the
        # runtime depend on what hardware happened to be visible; the CPU
        # resource is what the reported timing comparison is based on, and
        # matches the single-threaded MCMC the manuscript describes.
        "beast -seed {config[seed]} -beagle_CPU {input} 2>&1 | tee {log}"

rule variational_fit:
    input:
        fasta = lambda wildcards: all_models[wildcards.dataset]["alignment"],
        topology = wd / dataset_dir / "topology.nwk",
        starting_values = wd / dataset_dir / "starting-values.yaml",
        model_file = wd / dataset_dir / "model.yaml"
    params:
        starting_values_string = lambda wildcards, input: build_init_values_string(
            {k: v for k, v in yaml_input(input.starting_values).items() if k in models[wildcards.dataset].free_params()}
        )
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
        treeflow_vi run -s {config[seed]} \
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

wildcard_constraints:
    approx = "|".join(MULTI_RUN_APPROXES)

rule multi_run_variational_fit:
    input:
        fasta = lambda wildcards: all_models[wildcards.dataset]["alignment"],
        topology = wd / dataset_dir / "topology.nwk",
        starting_values = wd / dataset_dir / "starting-values.yaml",
        model_file = wd / dataset_dir / "model.yaml"
    params:
        starting_values_string = lambda wildcards, input: build_init_values_string(
            {k: v for k, v in yaml_input(input.starting_values).items() if k in models[wildcards.dataset].free_params()}
        ),
        num_steps = MULTI_RUN_NUM_STEPS,
        max_trace_coords = MULTI_RUN_MAX_TRACE_COORDS
    output:
        trace = wd / dataset_dir / "variational-multi-run" / "{approx}" / "trace-run{run}.pickle",
        samples = wd / dataset_dir / "variational-multi-run" / "{approx}" / "samples-run{run}.csv",
        tree_samples = wd / dataset_dir / "variational-multi-run" / "{approx}" / "tree-samples-run{run}.nexus"
    benchmark:
        wd / dataset_dir / "variational-multi-run" / "{approx}" / "benchmark-run{run}.txt"
    log:
        wd / dataset_dir / "variational-multi-run" / "{approx}" / "log-run{run}.txt"
    shell:
        '''
        treeflow_vi run -s {wildcards.run} \
            -va {wildcards.approx} \
            -i {input.fasta} \
            -m {input.model_file} \
            -t {input.topology} \
            -n {params.num_steps} \
            --learning-rate 0.001 \
            --init-values "{params.starting_values_string}" \
            --max-trace-coords {params.max_trace_coords} \
            --trace-output {output.trace} \
            --samples-output {output.samples} \
            --tree-samples-output {output.tree_samples} \
            --n-output-samples {config[n_variational_samples]} \
            2>&1 | tee {log}
        '''

rule multi_run_variational_samples:
    # Pools the per-run parameter samples into one table annotated with a `run`
    # column, matching the format the carnivores example notebook originally
    # used for its own multi-run posterior (examples/carnivores.ipynb), so the
    # same plotting code can show a pooled density plus an inter-run error band.
    input:
        samples = expand(
            wd / dataset_dir / "variational-multi-run" / "{approx}" / "samples-run{run}.csv",
            run=MULTI_RUN_SEEDS, allow_missing=True
        )
    output:
        wd / dataset_dir / "variational-multi-run" / "{approx}" / "samples.csv"
    run:
        frames = []
        for run, path in zip(MULTI_RUN_SEEDS, input.samples):
            frame = pd.read_csv(path)
            frame.insert(0, "run", run)
            frames.append(frame)
        pd.concat(frames, ignore_index=True).to_csv(output[0], index=False)

rule multi_run_variational_tree_samples:
    # Pools the per-run tree samples into one Nexus file so downstream summary
    # statistics (e.g. per-node height mean/SD) reflect both within- and
    # between-run variability.
    input:
        tree_samples = expand(
            wd / dataset_dir / "variational-multi-run" / "{approx}" / "tree-samples-run{run}.nexus",
            run=MULTI_RUN_SEEDS, allow_missing=True
        )
    output:
        wd / dataset_dir / "variational-multi-run" / "{approx}" / "tree-samples.nexus"
    run:
        combined = dendropy.TreeList()
        for path in input.tree_samples:
            trees = dendropy.TreeList.get(
                path=path, schema="nexus", taxon_namespace=combined.taxon_namespace
            )
            combined.extend(trees)
        combined.write(path=output[0], schema="nexus")

rule multi_run_variational_timing:
    input:
        benchmarks = expand(
            wd / dataset_dir / "variational-multi-run" / "{approx}" / "benchmark-run{run}.txt",
            run=MULTI_RUN_SEEDS, allow_missing=True
        )
    output:
        wd / dataset_dir / "variational-multi-run" / "{approx}" / "timing.csv"
    params:
        seeds = MULTI_RUN_SEEDS,
        num_steps = MULTI_RUN_NUM_STEPS
    run:
        pd.DataFrame([
            dict(
                seed=seed,
                num_steps=params.num_steps,
                elapsed_seconds=get_runtime_from_benchmark_file(benchmark_file)
            )
            for seed, benchmark_file in zip(params.seeds, input.benchmarks)
        ]).to_csv(output[0], index=False)

rule multi_run_variational:
    input:
        samples = rules.multi_run_variational_samples.output[0],
        tree_samples = rules.multi_run_variational_tree_samples.output[0],
        timing = rules.multi_run_variational_timing.output[0]

rule multi_run_variational_all:
    # Convenience target: every approximation's full 4-seed campaign on both
    # datasets. The main text only needs `main_approximation_runs` below; this
    # target additionally runs all 4 seeds of the two comparison families,
    # which the supplementary figures do not use.
    input:
        expand(
            wd / dataset_dir / "variational-multi-run" / "{approx}" / "{output}",
            dataset=MULTI_RUN_DATASETS,
            approx=MULTI_RUN_APPROXES,
            output=["samples.csv", "tree-samples.nexus", "timing.csv"]
        )

MAIN_APPROX = "root_full_rank"

rule main_approximation_runs:
    # Convenience target: the 4-seed root_full_rank campaign behind the
    # main-text marginals and tree figures, for both datasets.
    input:
        expand(
            wd / dataset_dir / "variational-multi-run" / MAIN_APPROX / "{output}",
            dataset=MULTI_RUN_DATASETS,
            output=["samples.csv", "tree-samples.nexus", "timing.csv"]
        )

rule approximation_comparison_runs:
    # Convenience target: the single (seed 1) run of each approximation family
    # behind the supplementary comparison figures and tables.
    input:
        expand(
            wd / dataset_dir / "variational-multi-run" / "{approx}" / "samples-run1.csv",
            dataset=MULTI_RUN_DATASETS,
            approx=MULTI_RUN_APPROXES
        )

# The carnivores example notebook is the source of the per-lineage kappa model
# results: the manuscript's kappa-vs-branch-age figure, the base/alt tree
# comparison figure, and the marginal likelihoods quoted in the text. It is a
# notebook rather than a CLI analysis because the per-lineage kappa model is
# specified through TreeFlow's Python API (it cannot be expressed in the YAML
# model format), which is itself one of the points the figure makes. This rule
# executes it non-interactively so those outputs are reproducible from the
# workflow.
rule carnivores_example_notebook:
    input:
        notebook = treeflow_dir / "examples" / "carnivores.ipynb",
        alignment = treeflow_dir / "examples" / "demo-data" / "carnivores.fasta",
        newick = treeflow_dir / "examples" / "demo-data" / "carnivores.newick"
    output:
        base_samples = treeflow_dir / "examples" / "demo-out" / "carnivores-base-samples.csv",
        base_trees = treeflow_dir / "examples" / "demo-out" / "carnivores-base-trees.nexus",
        alt_samples = treeflow_dir / "examples" / "demo-out" / "carnivores-alt-samples.csv",
        alt_trees = treeflow_dir / "examples" / "demo-out" / "carnivores-alt-trees.nexus",
        marginal_likelihoods = treeflow_dir / "examples" / "demo-out" / "carnivores-marginal-log-likelihoods.yaml"
    benchmark:
        wd / "carnivores" / "example-notebook-benchmark.txt"
    log:
        wd / "carnivores" / "example-notebook-log.txt"
    params:
        runner = treeflow_dir / "examples" / "run_example.py",
        python_executable = sys.executable
    shell:
        "{params.python_executable} {params.runner} {input.notebook} --inplace 2>&1 | tee {log}"


rule ml_fit:
    input:
        fasta = lambda wildcards: all_models[wildcards.dataset]["alignment"],
        topology = wd / dataset_dir / "topology.nwk",
        starting_values = wd / dataset_dir / "starting-values.yaml",
        model_file = wd / dataset_dir / "model.yaml"
    params:
        starting_values_string = lambda wildcards, input: build_init_values_string(
            {k: v for k, v in yaml_input(input.starting_values).items() if k in models[wildcards.dataset].free_params()}
        )
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

# The H3N2 timing comparison quoted in the manuscript uses the multi-run
# root_full_rank campaign -- the same runs the main-text figures come from --
# rather than the single 30,000-iteration full_rank fit the generic timing_data
# rule uses. The BEAST 2 side is read from recorded values (see
# config/h3n2-beast-timing.yaml) because re-running it is not part of this
# workflow.
ruleorder: flu_timing_data > timing_data

rule flu_timing_data:
    input:
        vi_benchmarks = expand(
            wd / "h3n2" / "variational-multi-run" / "root_full_rank" / "benchmark-run{run}.txt",
            run=MULTI_RUN_SEEDS
        ),
        vi_traces = expand(
            wd / "h3n2" / "variational-multi-run" / "root_full_rank" / "trace-run{run}.pickle",
            run=MULTI_RUN_SEEDS
        ),
        beast_timing = "config/h3n2-beast-timing.yaml"
    output:
        wd / "h3n2" / "timing-data.csv"
    run:
        from treeflow_pipeline.results import compute_variational_convergence

        runtimes = [get_runtime_from_benchmark_file(f) for f in input.vi_benchmarks]
        converged_iters = []
        vi_iters = []
        for trace_file in input.vi_traces:
            trace = pickle_input(trace_file)
            converged_iters.append(compute_variational_convergence(trace))
            vi_iters.append(len(trace.loss))
        if len(set(vi_iters)) != 1:
            raise ValueError(f"Runs have differing iteration counts: {vi_iters}")
        vi_iter = vi_iters[0]
        # Mean wall clock over the independent runs, and the iteration by which
        # the slowest of them had converged -- so "converged in N iterations
        # taking T" holds for every run rather than only the luckiest one.
        vi_runtime = sum(runtimes) / len(runtimes)
        vi_converged_iter = max(converged_iters)
        print(
            f"VI runtimes (s): {runtimes}; mean {vi_runtime:.1f}\n"
            f"Convergence iterations: {converged_iters}; using {vi_converged_iter}"
        )

        beast_timing = yaml_input(input.beast_timing)
        vi_df = pd.DataFrame(
            dict(iteration=[0, vi_converged_iter, vi_iter], value=[0.0, 1.0, 1.0])
        )
        vi_df = vi_df.assign(
            time=vi_runtime * vi_df["iteration"] / vi_iter,
            variable="converged",
            method="vi",
        )
        beast_df = pd.DataFrame(
            dict(
                iteration=[0, beast_timing["iterations"]],
                value=[0.0, beast_timing["min_ess"]],
            )
        )
        beast_df = beast_df.assign(
            time=beast_timing["runtime_seconds"]
            * beast_df["iteration"]
            / beast_timing["iterations"],
            variable="min_ess",
            method="beast",
        )
        pd.concat([vi_df, beast_df], ignore_index=True).to_csv(output[0], index=False)

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
        