from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
import pathlib
import pandas as pd
import treeflow_pipeline.model
from treeflow_pipeline.util import yaml_input, text_output
import treeflow_pipeline.manuscript

FTP = FTPRemoteProvider()

configfile: "config/ms-config.yaml"
model = treeflow_pipeline.model.Model(yaml_input(config["model_file"]))

def result(path):
    if config["remote_results"]:
        project_dir = pathlib.Path(config["project_dir"])
        return FTP.remote(str(project_dir / path))
    else:
        return pathlib.Path(path)

aggregate_result_dir = pathlib.Path(config["aggregate_result_dir"])
taxa_dir = config["taxa_dir"]
sequence_dir = config["sequence_dir"]
manuscript_dir = pathlib.Path("manuscript")

rule test:
    input:
        manuscript_dir / "figures" / "coverage.png"

rule ms:
    input:
        manuscript_dir / "out" / "main.pdf"

rule coverage_table:
    input: result(aggregate_result_dir / taxa_dir / sequence_dir / "coverage.csv")
    output: manuscript_dir / "tables" / "coverage-table.tex"
    run:
        pd.read_csv(input[0]).to_latex(output[0])

APPROXES = ["mean_field", "scaled"] # TODO: Where to store these in common?
methods = ["beast"] + expand("variational-samples-{approx}", approx=APPROXES)
stats = (
    [f"rate_stats.{stat}" for stat in ["mean", "coefficientOfVariation"]] +
    [f"tree.{stat}" for stat in ["height", "treeLength"]]    
)

rule coverage_plot:
    input:
        logs = expand(aggregate_result_dir / taxa_dir / sequence_dir / "{result}.log", result=methods, allow_missing=True),
        sim_trace = aggregate_result_dir / taxa_dir / "sim.log"
    output:
        manuscript_dir / "figures" / "coverage.png"
    run:
        treeflow_pipeline.manuscript.coverage_plot(
            dict(zip(methods, input.logs)),
            input.sim_trace,
            model,
            stats,
            output[0]
        )
            
rule compile_ms:
    input:
        rules.coverage_table.output[0],
        rules.coverage_plot.output[0],
        main = manuscript_dir / "main.tex"
    output:
        manuscript_dir / "out" / "main.pdf"
    params:
        output_dir =  lambda wildcards, output: output[0].parents[0]
    shell:
        "pdflatex -output-directory={params.output_dir} {input.main}"


