from treeflow_pipeline.util import yaml_input, yaml_output, text_input, text_output, beast_log_input
import treeflow_pipeline.simulation as sim
import treeflow_pipeline.templating as tem
import treeflow_pipeline.model as mod
import treeflow_pipeline.topology_inference as top
import pathlib

configfile: "config/sim-config.yaml"

wd = pathlib.Path(config["working_directory"])
model = mod.Model(yaml_input(config["model_file"]))

taxon_dir = "{taxon_count}taxa"

rule test_sim:
    input:
        "out/sim/10taxa/2seed/1000sites/sequences.fasta"

rule sampling_times:
    output:
        wd / taxon_dir / "sampling-times.yaml"
    run:
        yaml_output(sim.get_sampling_times(config, int(wildcards.taxon_count)), output[0])

seed_dir = "{seed}seed"

rule sample_prior:
    input:
        sampling_times = wd / taxon_dir / "sampling-times.yaml"
    output:
        wd / taxon_dir / seed_dir / "prior-sample.yaml"
    run:
        yaml_output(sim.sample_prior(yaml_input(input.sampling_times), model, int(wildcards.seed)), output[0])
    
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
        "beast -seed {config[seed]} {input}"

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
        "beast -seed {config[seed]} {input}"

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

sequence_dir = "{sequence_length}sites"

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
        "beast -seed {config[seed]} {input}"

rule fasta_sim:
    input:
        "{wd}/sequences.xml"
    output:
        "{wd}/sequences.fasta"
    run:
        sim.convert_simulated_sequences(input[0], output[0], 'fasta')
