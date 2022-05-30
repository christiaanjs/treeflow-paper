import pathlib
import pandas as pd
import treeflow_pipeline.model
from treeflow_pipeline.util import yaml_input, text_input, text_output
import treeflow_pipeline.manuscript

import treeflow_benchmarks
treeflow_benchmarks_dir = pathlib.Path(treeflow_benchmarks.__file__).parents[1]

configfile: "config/ms-config.yaml"
model = treeflow_pipeline.model.Model(yaml_input(config["model_file"]))

def result(path):
    if config["remote_results"]:
        raise NotImplemented("Remote files note implemented")
    else:
        return pathlib.Path(path)

aggregate_result_dir = pathlib.Path(config["aggregate_result_dir"])
taxa_dir = config["taxa_dir"]
sequence_dir = config["sequence_dir"]
manuscript_dir = pathlib.Path("manuscript")

rule ms:
    input:
        manuscript_dir / "out" / "treeflow.pdf"

APPROXES = ["mean_field", "scaled"] # TODO: Where to store these in common?
methods = ["beast"] + expand("variational-samples-{approx}", approx=APPROXES)
stats = (
    [f"rate_stats.{stat}" for stat in ["mean", "coefficientOfVariation"]] +
    [f"tree.{stat}" for stat in ["height", "treeLength"]]    
)

rule coverage_table:
    input: result(aggregate_result_dir / taxa_dir / sequence_dir / "coverage.csv")
    output: manuscript_dir / "tables" / "coverage-table.tex"
    run:
        treeflow_pipeline.manuscript.coverage_table(input[0], list(model.free_params().keys()) + ["height", "rate"] + stats, output[0])

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

tex_template = manuscript_dir / "tex" / ("plos-template.j2.tex" if config["submission"] else "plain-template.j2.tex")

rule template_relaxed_clock_ms:
    input:
        coverage_table = rules.coverage_table.output[0],
        coverage_plot = rules.coverage_plot.output[0],
        template = tex_template,
        body_template = manuscript_dir / "tex" / "main.j2.tex"
    output:
        manuscript_dir / "out" / "main.tex"
    run:
        text_output(
            treeflow_pipeline.manuscript.build_manuscript(
                input.template,
                input.body_template,
                dict(coverage=input.coverage_plot),
                dict(coverage=input.coverage_table),
                dict(),
                submission=config["submission"]
            ),
            output[0]
        )

rule template_treeflow_ms:
    input:
        coverage_table = rules.coverage_table.output[0],
        coverage_plot = rules.coverage_plot.output[0],
        template = tex_template,
        body_template = manuscript_dir / "tex" / "treeflow.j2.tex",
        treeflow_benchmarks_config = treeflow_benchmarks_dir / "config.yaml",
    output:
        manuscript_dir / "out" / "treeflow.tex"
    run:
        text_output(
            treeflow_pipeline.manuscript.build_manuscript(
                input.template,
                input.body_template,
                dict(),
                dict(),
                treeflow_pipeline.manuscript.get_treeflow_manuscript_vars(yaml_input(input.treeflow_benchmarks_config)),
                submission=config["submission"]
            ),
            output[0]
        )


            
rule compile_ms:
    input:
        main = manuscript_dir / "out" / "{manuscript}.tex",
        bib = manuscript_dir / "tex" / "main.bib",
        bst = manuscript_dir / "tex" / "plos2015.bst"
    output:
        manuscript_dir / "out" / "{manuscript}.pdf"
    params:
        output_dir =  lambda _, output: pathlib.Path(output[0]).parents[0],
        aux_file = lambda _, output: pathlib.Path(output[0]).with_suffix(".aux"),
        bib_dir = lambda _, input: pathlib.Path(input.bib).parents[0],
        bst_dir = lambda _, input: pathlib.Path(input.bst).parents[0]
    shell:
        """
        pdflatex -output-directory={params.output_dir} {input.main}
        export BIBINPUTS={params.bib_dir}
        export BSTINPUTS={params.bst_dir}
        bibtex {params.aux_file}
        pdflatex -output-directory={params.output_dir} {input.main}
        pdflatex -output-directory={params.output_dir} {input.main}
        """


