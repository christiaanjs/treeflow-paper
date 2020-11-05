from treeflow_pipeline.util import yaml_input, yaml_output, text_input, text_output, beast_log_input
import treeflow_pipeline.simulation as sim
import treeflow_pipeline.templating as tem

rule sim_config:
    input:
        "config/sim-config.yaml"
    output:
        "out/sim/sequence_length{sequence_length}/sim-config.yaml"
    run:
        yaml_output(sim.prepare_config(yaml_input(input[0]), wildcards), output[0])


rule sampling_times:
    input:
        "out/sim/{dataset}/sim-config.yaml"
    output:
        "out/sim/{dataset}/sampling-times.yaml"
    run:
        yaml_output(sim.get_sampling_times(yaml_input(input[0])), output[0])

rule pop_size:
    input:
        prior_config = "config/prior-params.yaml",
        sim_config = "out/sim/{dataset}/sim-config.yaml"
    output:
        "out/sim/{dataset}/pop-size.yaml"
    run:
        yaml_output(sim.get_pop_size(yaml_input(input.sim_config), yaml_input(input.prior_config)), output[0])
    
rule tree_sim_xml:
    input:
        config = "out/sim/{dataset}/sim-config.yaml",
        sampling_times = "out/sim/{dataset}/sampling-times.yaml",
        pop_size = "out/sim/{dataset}/pop-size.yaml"
    output:
        "out/sim/{dataset}/tree-sim.xml"
    run:
        text_output(
            tem.build_tree_sim(
                yaml_input(input.config),
                yaml_input(input.sampling_times),
                yaml_input(input.pop_size),
                output[0]
            ),
            output[0]
        )

rule tree_sim: # TODO: Seed
    input:
        "out/sim/{dataset}/tree-sim.xml"
    output:
        "out/sim/{dataset}/tree-sim.trees"
    shell:
        "beast {input}"

rule tree_sim_newick:
    input:
        "out/sim/{dataset}/tree-sim.trees"
    output:
        "out/sim/{dataset}/tree-sim.newick"
    run:
        top.convert_tree(input[0], 'nexus', 'newick')


rule rate_sd:
    input:
        prior_config = "config/prior-params.yaml",
        sim_config = "out/sim/{dataset}/sim-config.yaml"
    output:
        "out/sim/{dataset}/rate-sd.yaml"
    run:
        yaml_output(sim.get_rate_sd(yaml_input(input.sim_config), yaml_input(input.prior_config)), output[0])

rule branch_rate_sim_xml:
    input:
        rate_sd = "out/sim/{dataset}/rate-sd.yaml",
        tree = "out/sim/{dataset}/tree-sim.newick"
    output:
        "out/sim/{dataset}/branch-rate-sim.xml"
    run:
        text_output(
            tem.build_branch_rate_sim(
                text_input(input.tree),
                yaml_input(input.rate_sd),
                output[0]
            ),
            output[0]
        )

rule branch_rate_sim:
    input:
        "out/sim/{dataset}/branch-rate-sim.xml"
    output:
        "out/sim/{dataset}/branch-rate-sim.log",
        "out/sim/{dataset}/branch-rate-sim.trees"
    shell:
        "beast {input}"

rule convert_branch_rates:
    input:
        "out/sim/{dataset}/branch-rate-sim.log"
    output:
        "out/sim/{dataset}/branch-rates.yaml"
    run:
        yaml_output(
            sim.parse_branch_rates(beast_log_input(input[0])),
            output[0]
        )

rule clock_rate:
    input:
        prior_config = "config/prior-params.yaml",
        sim_config = "out/sim/{dataset}/sim-config.yaml"
    output:
        "out/sim/{dataset}/clock-rate.yaml"
    run:
        yaml_output(sim.get_clock_rate(yaml_input(input.sim_config), yaml_input(input.prior_config)), output[0])


rule sequence_sim_xml:
    input:
        tree = "out/sim/{dataset}/tree-sim.newick",
        sim_config = "out/sim/{dataset}/sim-config.yaml",
        clock_rate = "out/sim/{dataset}/clock-rate.yaml",
        branch_rates = "out/sim/{dataset}/branch-rates.yaml",
        sampling_times = "out/sim/{dataset}/sampling-times.yaml",
    output:
        "out/sim/{dataset}/sim-seq.xml"
    run:
        text_output(
            tem.build_sequence_sim(
                yaml_input(input.sim_config),
                text_input(input.tree),
                yaml_input(input.clock_rate),
                yaml_input(input.branch_rates),
                yaml_input(input.sampling_times),
                output[0]
            ),
            output[0]
        )

rule sequence_sim: # TODO: Seed
    input:
        "out/sim/{dataset}/sim-seq.xml"
    output:
        "out/sim/{dataset}/sequences.xml"
    shell:
        "beast {input}"

rule fasta_sim:
    input:
        "out/sim/{dataset}/sequences.xml"
    output:
        "out/sim/{dataset}/sequences.fasta"
    run:
        sim.convert_simulated_sequences(input[0], output[0], 'fasta')
        
# rule sim_data_dir:
#     output:
#         directory("data/sim/{dataset}")
#     shell:
#         "mkdir -p {output}"

rule fasta_copy:
    input:
        "out/sim/{dataset}/sequences.fasta"
    output:
        "data/sim/{dataset}.fasta"
    shell:
        "cp {input} {output}"