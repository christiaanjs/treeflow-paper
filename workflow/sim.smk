from treeflow_pipeline.util import yaml_input, yaml_output, text_input, text_output, beast_log_input, sequence_input, pickle_input, pickle_output
import treeflow_pipeline.simulation as sim
import treeflow_pipeline.templating as tem
import treeflow_pipeline.model as mod
import treeflow_pipeline.topology_inference as top
import treeflow_pipeline.results as res
import treeflow_pipeline.priors as priors
import pathlib

configfile: "config/sim-config.yaml"

wd = pathlib.Path(config["working_directory"])
#model = mod.Model(yaml_input(config["model_file"]))
beast_config = yaml_input(config["beast_config"])

TAXON_COUNTS = [20]
SEQUENCE_LENGTHS = [10000]
APPROXES = ["mean_field", "scaled"]
SEEDS = list(range(1, config["replicates"]+1))
DEMO_SEED = 4

taxon_dir = "{taxon_count}taxa"
seed_dir = "{seed}seed"
sequence_dir = "{sequence_length}sites"
aggregate_dir = "aggregate"

model_file = config["model_file"]

rule test:
    input:
        model_file

rule well_calibrated_study:
    input:
        expand(str(wd / aggregate_dir / taxon_dir / sequence_dir / "coverage.csv"), sequence_length=SEQUENCE_LENGTHS, taxon_count=TAXON_COUNTS),
        expand(str(wd / aggregate_dir / taxon_dir / sequence_dir / "coverage-plot.csv"), sequence_length=SEQUENCE_LENGTHS, taxon_count=TAXON_COUNTS)


rule demo:
    input:
        expand(str(wd / "10taxa" / seed_dir / sequence_dir / "plot-posterior-relaxed-{clock_approx}.html"), sequence_length=SEQUENCE_LENGTHS, clock_approx=APPROXES, seed=[DEMO_SEED])

rule sampling_times:
    output:
        wd / taxon_dir / "sampling-times.yaml"
    run:
        yaml_output(sim.get_sampling_times(config, int(wildcards.taxon_count)), output[0])

rule model_from_spec:
    input:
        model_spec = config["model_spec_file"]
    output:
        model = model_file
    run:
        yaml_output(
            priors.get_priors_from_spec(yaml_input(input.model_spec)),
            output.model
        )

rule sample_prior:
    input:
        sampling_times = wd / taxon_dir / "sampling-times.yaml",
        model = model_file
    output:
        wd / taxon_dir / seed_dir / "prior-sample.yaml"
    run:
        yaml_output(sim.sample_prior(yaml_input(input.sampling_times), yaml_input(input.model), int(wildcards.seed)), output[0])
    
rule tree_sim_xml:
    input:
        sampling_times = wd / taxon_dir / "sampling-times.yaml",
        prior_sample = wd / taxon_dir / seed_dir / "prior-sample.yaml"
    output:
        wd / taxon_dir / seed_dir / "tree-sim.xml"
    run:
        text_output(
            tem.build_tree_sim(
                config,
                yaml_input(input.sampling_times),
                yaml_input(input.prior_sample),
                output[0]
            ),
            output[0]
        )

rule tree_sim:
    input:
        wd / taxon_dir / seed_dir / "tree-sim.xml"
    output:
        wd / taxon_dir / seed_dir / "tree-sim.trees"
    shell:
        "beast -seed {wildcards.seed} {input}"

rule tree_sim_newick:
    input:
        "{dir}/tree-sim.trees"
    output:
        "{dir}/tree-sim.newick"
    run:
        top.convert_tree(input[0], 'nexus', output[0], 'newick')

rule branch_rate_sim_xml:
    input:
        prior_sample = wd / taxon_dir / seed_dir / "prior-sample.yaml",
        tree = wd / taxon_dir / seed_dir / "tree-sim.newick"
    output:
        wd / taxon_dir / seed_dir / "branch-rate-sim.xml"
    run:
        text_output(
            tem.build_branch_rate_sim(
                text_input(input.tree),
                yaml_input(input.prior_sample),
                output[0]
            ),
            output[0]
        )

rule branch_rate_sim:
    input:
        wd / taxon_dir / seed_dir / "branch-rate-sim.xml"
    output:
        wd / taxon_dir / seed_dir / "branch-rate-sim.log",
        wd / taxon_dir / seed_dir / "branch-rate-sim.trees"
    shell:
        "beast -seed {wildcards.seed} {input}"

rule convert_branch_rates:
    input:
        wd / taxon_dir / seed_dir / "branch-rate-sim.log"
    output:
        wd / taxon_dir / seed_dir / "branch-rates.yaml"
    run:
        yaml_output(
            sim.parse_branch_rates(beast_log_input(input[0])),
            output[0]
        )

rule sim_trace:
    input:
        prior_sample = wd / taxon_dir / seed_dir / "prior-sample.yaml",
        tree = wd / taxon_dir / seed_dir / "tree-sim.newick",
        rates = wd / taxon_dir / seed_dir / "branch-rate-sim.log"
    output:
        wd / taxon_dir / seed_dir / "sim.log"
    run:
        sim.build_sim_trace(input.tree, yaml_input(input.prior_sample), output[0], rate_trace=input.rates)

rule sequence_sim_xml:
    input:
        tree = wd / taxon_dir / seed_dir / "tree-sim.newick",
        prior_sample = wd / taxon_dir / seed_dir / "prior-sample.yaml",
        branch_rates = wd / taxon_dir / seed_dir / "branch-rates.yaml",
        sampling_times = wd / taxon_dir / "sampling-times.yaml"
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "sim-seq.xml"
    run:
        text_output(
            tem.build_sequence_sim(
                config,
                text_input(input.tree),
                yaml_input(input.prior_sample),
                yaml_input(input.branch_rates),
                yaml_input(input.sampling_times),
                int(wildcards.sequence_length),
                output[0]
            ),
            output[0]
        )

rule sequence_sim:
    input:
        wd / taxon_dir / seed_dir / sequence_dir / "sim-seq.xml"
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "sequences.xml"
    shell:
        "beast -seed {wildcards.seed} {input}"

rule fasta_sim:
    input:
        "{wd}/sequences.xml"
    output:
        "{wd}/sequences.fasta"
    run:
        sim.convert_simulated_sequences(input[0], output[0], 'fasta')

rule sim_starting_values:
    input:
        prior_sample = wd / taxon_dir / seed_dir / "prior-sample.yaml"
    output:
        wd / taxon_dir / seed_dir / "starting-values.yaml"
    run:
        yaml_output(dict(
            frequencies=config["frequencies"],
            kappa=config["kappa"],
            **yaml_input(input.prior_sample)
        ), output[0])

rule beast_xml:
    input:
        fasta = wd / taxon_dir / seed_dir / sequence_dir / "sequences.fasta",
        tree = wd / taxon_dir / seed_dir / "tree-sim.newick",
        starting_values = wd / taxon_dir / seed_dir / "starting-values.yaml",
        beast_config = config["beast_config"]
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "beast.xml"
    run:
        text_output(tem.build_beast_analysis(
            sequence_input(input.fasta),
            text_input(input.tree),
            yaml_input(input.starting_values),
            yaml_input,
            beast_config,
            output[0]
        ), output[0])

rule beast_run:
    input:
        wd / taxon_dir / seed_dir / sequence_dir / "beast.xml"
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "beast.log",
        wd / taxon_dir / seed_dir / sequence_dir / "beast.trees"
    shell:
        "beast -seed {wildcards.seed} {input}"

rule beast_results:
    input:
        topology = wd / taxon_dir / seed_dir / "tree-sim.newick",
        trees = wd / taxon_dir / seed_dir / sequence_dir / "beast.trees",
        trace = wd / taxon_dir / seed_dir / sequence_dir / "beast.log",
        model = model_file
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "beast.pickle"
    run:
        pickle_output(
            res.process_beast_results(
                input.trees,
                input.trace,
                input.topology,
                beast_config,
                yaml_input(input.model)
            ),
            output[0]
        )



rule variational_fit:
    input:
        fasta = wd / taxon_dir / seed_dir / sequence_dir / "sequences.fasta",
        tree = wd / taxon_dir / seed_dir / "tree-sim.newick",
        starting_values = wd / taxon_dir / seed_dir / "starting-values.yaml",
        model = model_file
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "variational-fit-{clock_approx}.pickle"
    shell:
        """
        treeflow_pipeline -s {wildcards.seed} \
            {input.fasta} {input.model} {output} \
            variational-fit \
            -t {input.tree} \
            -s {input.starting_values} \
            -a {wildcards.clock_approx} \
            --config {config[vi_config]}
        """

rule variational_samples: # TODO: Include this in CLI
    input:
        fit = wd / taxon_dir / seed_dir / sequence_dir / "variational-fit-{clock_approx}.pickle",
        topology = wd / taxon_dir / seed_dir / "tree-sim.newick",
        starting_values = wd / taxon_dir / seed_dir / "starting-values.yaml",
        model = model_file
    output:
        trace = wd / taxon_dir / seed_dir / sequence_dir / "variational-samples-{clock_approx}.log",
        trees = wd / taxon_dir / seed_dir / sequence_dir / "variational-samples-{clock_approx}.trees",
        samples = wd / taxon_dir / seed_dir / sequence_dir / "variational-samples-{clock_approx}.pickle"
    run:
        pickle_output(
            res.get_variational_samples(
                pickle_input(input.fit),
                input.topology,
                yaml_input(input.model),
                wildcards.clock_approx,
                yaml_input(input.starting_values),
                output.trace,
                output.trees,
                wildcards.seed
            ),
            output.samples
        )

rule relaxed_plot:
    input:
        topology = wd / taxon_dir / seed_dir / "tree-sim.newick",
        beast_result = wd / taxon_dir / seed_dir / sequence_dir / "beast.pickle",
        variational_fit = wd / taxon_dir / seed_dir / sequence_dir / "variational-fit-{clock_approx}.pickle",
        notebook = "notebook/plot-posterior-relaxed.ipynb",
        model = model_file
    output:
        notebook = wd / taxon_dir / seed_dir / sequence_dir / "plot-posterior-relaxed-{clock_approx}.ipynb",
        correlation_plot = wd / taxon_dir / seed_dir / sequence_dir / "rate-correlations-{clock_approx}.png",
        marginal_plot = wd / taxon_dir / seed_dir / sequence_dir / "marginals-{clock_approx}.png",
        rate_marginal_plot = wd / taxon_dir / seed_dir / sequence_dir / "rate-marginals-{clock_approx}.png"
    shell:
        """
        papermill {input.notebook} {output.notebook} \
            -p beast_result_file {input.beast_result} \
            -p topology_file {input.topology} \
            -p model_file {input.model} \
            -p variational_fit_file {input.variational_fit} \
            -p clock_approx {wildcards.clock_approx} \
            -p correlation_plot_out_file {output.correlation_plot} \
            -p marginal_plot_out_file {output.marginal_plot} \
            -p rate_marginal_plot_out_file {output.rate_marginal_plot}
        """

rule relaxed_report:
    input:
        wd / taxon_dir / seed_dir / sequence_dir / "plot-posterior-relaxed-{clock_approx}.ipynb"
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "plot-posterior-relaxed-{clock_approx}.html"
    shell:
        "jupyter nbconvert --to html {input}"

rule aggregate_sim_trees:
    input:
        expand(wd / taxon_dir / seed_dir / "tree-sim.trees", seed=SEEDS, allow_missing=True)
    output:
        wd / aggregate_dir / taxon_dir / "tree-sim.trees"
    run:
        sim.aggregate_trees(input, "nexus", output[0], "nexus")

rule log_analyser:
    input:
        expand(wd / taxon_dir / seed_dir / sequence_dir / "{result}.log", seed=SEEDS, allow_missing=True)
    output:
        wd / aggregate_dir / taxon_dir / sequence_dir / "{result}.log"
    params:
        burn_in = int(beast_config["burn_in"] * 100)
    shell:
        "loganalyser -b {params.burn_in} -oneline {input} > {output}"

rule aggregate_sim_trace:
    input:
        expand(wd / taxon_dir / seed_dir / "sim.log", seed=SEEDS, allow_missing=True)
    output:
        wd / aggregate_dir / taxon_dir / "sim.log"
    run:
        sim.aggregate_sim_traces(input, output[0])

free_params = priors.get_free_params_from_model_spec(yaml_input(config["model_spec_file"]))
stats = (
    free_params + 
    [f"rate_stats.{stat}" for stat in ["mean", "coefficientOfVariation"]] +
    [f"tree.{stat}" for stat in ["height", "treeLength"]]    
)

rule coverage:
    input:
        sim_trace = wd / aggregate_dir / taxon_dir / "sim.log",
        analyser_trace = wd / aggregate_dir / taxon_dir / sequence_dir / "{result}.log"
    output:
        report = wd / aggregate_dir / taxon_dir / sequence_dir / "{result}" / "coverage.html",
        stats = [wd / aggregate_dir / taxon_dir / sequence_dir / "{result}" / f"{stat}.tsv" for stat in stats]
    params:
        output_dir =  lambda wildcards, output: pathlib.Path(output.report).parents[0]
    shell:
        """
        applauncher CoverageCalculator \
            -log {input.sim_trace} \
            -logAnalyser {input.analyser_trace} \
            -out {params.output_dir} \
            -skip 0
        """

rule tree_annotator: # TODO: Methods in their own directories?
    input:
        trees = wd / taxon_dir / seed_dir / sequence_dir / "{result}.trees",
        tree_sim = wd / taxon_dir / seed_dir / "branch-rate-sim.trees"
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "mcc-{result}.trees"
    params:
        burn_in = int(beast_config["burn_in"] * 100)
    shell:
        "treeannotator -b {params.burn_in} -target {input.tree_sim} {input.trees} {output}"

rule tree_coverage:
    input:
        mcc_trees = expand(wd / taxon_dir / seed_dir / sequence_dir / "mcc-{result}.trees", seed=SEEDS, allow_missing=True),
        rate_sims = expand(wd / taxon_dir / seed_dir / "branch-rate-sim.trees", seed=SEEDS, allow_missing=True)
    params:
        mcc_file_template = str(wd / taxon_dir / "$(n)seed" / sequence_dir / "mcc-{result}.trees"),
        tree_file_template = str(wd / taxon_dir / "$(n)seed" / "branch-rate-sim.trees"),
        from_index = min(SEEDS),
        to_index = max(SEEDS)
    output:
        wd / aggregate_dir / taxon_dir / sequence_dir / "{result}" / "tree-coverage.txt"
    shell:
        """
        applauncher MCCTreeComparator \
            -mcc '{params.mcc_file_template}' \
            -tree '{params.tree_file_template}' \
            -from {params.from_index} \
            -to {params.to_index} \
            -out {output}
        """

rule method_coverage_table:
    input:
        coverage_stats = rules.coverage.output.stats,
        tree_coverage = rules.tree_coverage.output[0]
    output:
        wd / aggregate_dir / taxon_dir / sequence_dir / "{result}" / "coverage.csv"
    run:
        res.build_method_coverage_table(wildcards.result, dict(zip(stats, input.coverage_stats)), text_input(input.tree_coverage), output[0])

methods = expand("variational-samples-{approx}", approx=APPROXES) + ["beast"]
rule coverage_table:
    input:
        expand(rules.method_coverage_table.output[0], result=methods, allow_missing=True)
    output:
        wd / aggregate_dir / taxon_dir / sequence_dir / "coverage.csv"
    run:
        res.aggregate_coverage_tables(input, output[0])
    
rule method_coverage_plot_table:
    input:
        coverage_stats = rules.coverage.output.stats
    output:
        wd / aggregate_dir / taxon_dir / sequence_dir / "{result}" / "coverage-plot.csv"
    run:
        res.build_method_coverage_plot_table(wildcards.result, dict(zip(stats, input.coverage_stats)), output[0])

rule coverage_plot_table:
    input:
        expand(rules.method_coverage_plot_table.output[0], result=methods, allow_missing=True)
    output:
        wd / aggregate_dir / taxon_dir / sequence_dir / "coverage-plot.csv"
    run:
        res.aggregate_coverage_tables(input, output[0])