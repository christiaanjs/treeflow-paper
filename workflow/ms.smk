import pathlib
import sys
import pandas as pd
import treeflow_pipeline.model
from treeflow_pipeline.util import yaml_input, text_input, text_output
import treeflow_pipeline.manuscript

import treeflow
import treeflow_benchmarks
treeflow_benchmarks_dir = pathlib.Path(treeflow_benchmarks.__file__).parents[1]
treeflow_dir = pathlib.Path(treeflow.__file__).parents[1]

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
out_dir = pathlib.Path("out")
submission_dir = manuscript_dir / "submission"
minted_cache_dir = "minted-cache"

rule ms:
    input:
        manuscript_dir / "out" / "submission.zip",
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

tex_templates = dict(
    plain="plain-template.j2.tex",
    plos="plos-template.j2.tex",
    systematic_biology="samplebibtex_sse.j2.tex"
)
tex_template = manuscript_dir / "tex" / tex_templates[config["tex_template"]]

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

rule benchmark_summary_table:
    input:
        plot_data = treeflow_benchmarks_dir / "out" / "plot-data.csv",
        fit_table = treeflow_benchmarks_dir / "out" / "fit-table.csv"
    output:
        tex = manuscript_dir / "tables" / "benchmark-table.tex"
    run:
        treeflow_pipeline.manuscript.benchmark_summary_table(input.plot_data, input.fit_table, output.tex)

rule benchmark_plot:
    input:
        plot_data = treeflow_benchmarks_dir / "out" / "plot-data.csv",
    output:
        plot = manuscript_dir / "figures" / "benchmark-log-scale-plot.png"
    params:
        python_executable = sys.executable
    script:
        "../scripts/improved-benchmark-plot.R"

rule data_marginals_plot:
    input:
        vi_samples = out_dir / "{dataset}" / "variational-samples.csv",
        beast_samples = out_dir / "{dataset}" / "beast.log"
    output:
        plot = manuscript_dir / "figures" / "{dataset}-marginals.png"
    params:
        python_executable = sys.executable
    script:
        "../scripts/improved-marginals-plot.R"

rule carnivores_kappa_plot:
    input:
        tree_samples = treeflow_dir / "examples" / "demo-out" / "carnivores-alt-trees.nexus"
    output:
        manuscript_dir / "figures" / "carnivores-kappa.png"
    script:
        "../scripts/carnivores-kappa-plot.R"

rule carnivores_tree_plot:
    input:
        alt_tree_samples = treeflow_dir / "examples" / "demo-out" / "carnivores-alt-trees.nexus",
        base_tree_samples = treeflow_dir / "examples" / "demo-out" / "carnivores-base-trees.nexus"
    output:
        manuscript_dir / "figures" / "carnivores-model-trees.png"
    script:
        "../scripts/carnivores-tree-plot.R"

rule data_tree_plot:
    input:
        vi_tree_samples = out_dir / "{dataset}" / "variational-tree-samples.nexus",
        beast_tree_samples = out_dir / "{dataset}" / "beast.trees"
    output:
        plot = manuscript_dir / "figures" / "{dataset}-trees.png"
    script:
        "../scripts/data-tree-plot.R"

rule nf_data_tree_plot:
    input:
        vi_tree_samples = out_dir / "{dataset}" / "nf-tree-samples.nexus",
        beast_tree_samples = out_dir / "{dataset}" / "beast.trees"
    output:
        plot = manuscript_dir / "figures" / "{dataset}-nf-trees.png"
    script:
        "../scripts/data-tree-plot.R"

rule template_treeflow_ms:
    input:
        template = tex_template,
        body_template = manuscript_dir / "tex" / "treeflow.j2.tex",
        treeflow_benchmarks_config = treeflow_benchmarks_dir / "config.yaml",
        benchmark_plot = rules.benchmark_plot.output.plot,
        benchmark_summary_table = rules.benchmark_summary_table.output[0],
        carnivores_marginals_plot = manuscript_dir / "figures" / "carnivores-marginals.png",
        carnivores_kappa_plot = rules.carnivores_kappa_plot.output[0],
        carnivores_tree_plot = rules.carnivores_tree_plot.output[0],
        flu_marginals_plot = manuscript_dir / "figures" / f"{config['flu_dataset']}-marginals.png",
        flu_tree_plot = manuscript_dir / "figures" / f"{config['flu_dataset']}-trees.png",
        flu_timing_csv = out_dir / config["flu_dataset"] / "timing-data.csv",
        flu_model_file = out_dir / config["flu_dataset"] / "model.yaml",
        flu_tree_file = out_dir / config["flu_dataset"] / "topology.nwk",
        bib = manuscript_dir / "tex" / "main.bib",
    output:
        manuscript_dir / "out" / "treeflow.tex"
    params:
        output_dir =  lambda _, output: pathlib.Path(output[0]).parents[0],
    run:
        text_output(
            treeflow_pipeline.manuscript.build_manuscript(
                input.template,
                input.body_template,
                figures_dict=dict(
                    benchmark=input.benchmark_plot,
                    carnivores_marginals=input.carnivores_marginals_plot,
                    carnivores_kappa=input.carnivores_kappa_plot,
                    carnivores_tree=input.carnivores_tree_plot,
                    flu_marginals=input.flu_marginals_plot,
                    flu_tree=input.flu_tree_plot
                ),
                tables_dict=dict(benchmark_summary=input.benchmark_summary_table),
                vars=dict(
                    treeflow_pipeline.manuscript.get_treeflow_manuscript_vars(
                        yaml_input(input.treeflow_benchmarks_config),
                        timing_csv_file=input.flu_timing_csv,
                        flu_model_file=input.flu_model_file,
                        flu_tree_file=input.flu_tree_file,
                        minted_cache_dir=params.output_dir / "minted-cache",
                        bibliography_file=input.bib
                    ),
                    output_dir = params.output_dir,
                ),
                submission=config["submission"]
            ),
            output[0]
        )

rule treeflow_submission_dir:
    input:
        benchmark_plot = rules.benchmark_plot.output.plot,
        carnivores_marginals_plot = manuscript_dir / "figures" / "carnivores-marginals.png",
        carnivores_kappa_plot = rules.carnivores_kappa_plot.output[0],
        carnivores_tree_plot = rules.carnivores_tree_plot.output[0],
        flu_marginals_plot = manuscript_dir / "figures" / f"{config['flu_dataset']}-marginals.png",
        flu_tree_plot = manuscript_dir / "figures" / f"{config['flu_dataset']}-trees.png",
        flu_model_file = out_dir / config["flu_dataset"] / "model.yaml",
        bib = manuscript_dir / "tex" / "main.bib",
    output:
        benchmark_plot = submission_dir / "benchmark-plot.png",
        carnivores_marginals_plot = submission_dir / "carnivores-marginals.png",
        carnivores_kappa_plot = submission_dir / "carnivores-kappa.png",
        carnivores_tree_plot = submission_dir / "carnivores-trees.png",
        flu_marginals_plot = submission_dir / "flu-marginals.png",
        flu_tree_plot = submission_dir / "flu-trees.png",
        flu_model_file = submission_dir / "flu-model.yaml",
        bib = submission_dir / "treeflow.bib",
    run:
        for key in input.keys():
            shell(f"cp {input[key]} {output[key]}")

rule template_treeflow_submission_ms:
    input:
        template = tex_template,
        body_template = manuscript_dir / "tex" / "treeflow.j2.tex",
        treeflow_benchmarks_config = treeflow_benchmarks_dir / "config.yaml",
        benchmark_plot = rules.treeflow_submission_dir.output.benchmark_plot,
        benchmark_summary_table = rules.benchmark_summary_table.output[0],
        carnivores_marginals_plot = rules.treeflow_submission_dir.output.carnivores_marginals_plot,
        carnivores_kappa_plot = rules.treeflow_submission_dir.output.carnivores_kappa_plot,
        carnivores_tree_plot = rules.treeflow_submission_dir.output.carnivores_tree_plot,
        flu_marginals_plot = rules.treeflow_submission_dir.output.flu_marginals_plot,
        flu_tree_plot = rules.treeflow_submission_dir.output.flu_tree_plot,
        flu_timing_csv = out_dir / config["flu_dataset"] / "timing-data.csv",
        flu_model_file = rules.treeflow_submission_dir.output.flu_model_file,
        flu_tree_file = out_dir / config["flu_dataset"] / "topology.nwk",
        bib = rules.treeflow_submission_dir.output.bib
    output:
        submission_dir / "treeflow{minted_state,.*}.tex"
    params:
        output_dir =  lambda _, output: pathlib.Path(output[0]).parents[0],
    run:
        text_output(
            treeflow_pipeline.manuscript.build_manuscript(
                input.template,
                input.body_template,
                figures_dict={ key: pathlib.Path(path).relative_to(params.output_dir)
                    for key, path in dict(
                        benchmark=input.benchmark_plot,
                        carnivores_marginals=input.carnivores_marginals_plot,
                        carnivores_kappa=input.carnivores_kappa_plot,
                        carnivores_tree=input.carnivores_tree_plot,
                        flu_marginals=input.flu_marginals_plot,
                        flu_tree=input.flu_tree_plot
                    ).items() },
                tables_dict=dict(benchmark_summary=input.benchmark_summary_table),
                vars=dict(
                    treeflow_pipeline.manuscript.get_treeflow_manuscript_vars(
                        yaml_input(input.treeflow_benchmarks_config),
                        timing_csv_file=input.flu_timing_csv,
                        flu_model_file=pathlib.Path(input.flu_model_file).relative_to(params.output_dir),
                        flu_tree_file=input.flu_tree_file,
                        minted_cache_dir=minted_cache_dir,
                        bibliography_file=input.bib,
                        frozen_minted_cache=(wildcards.minted_state != "-mintedcache")
                    ),
                    output_dir = ".",
                ),
                submission=config["submission"]
            ),
            output[0]
        )

rule treeflow_submission_minted_cache:
    input:
        tex = submission_dir / "treeflow-mintedcache.tex"
    output:
        directory(submission_dir / minted_cache_dir)
    params:
        submission_dir = str(submission_dir),
        tex = lambda _, input: pathlib.Path(input.tex).relative_to(submission_dir)
    shell:
        """
        cd {params.submission_dir}
        pdflatex --shell-escape {params.tex}
        """

# rule treeflow_submission_minted_sty:
#     output:
#         submission_dir / "minted.sty"
#     shell:
#         "curl {config[minted_sty_url]} -o {output}"

rule treeflow_submission_bbl:
    input:
        tex = submission_dir / "treeflow.tex",
        minted_cache = rules.treeflow_submission_minted_cache.output,
    output:
        aux = submission_dir / "treeflow.aux",
        bbl = submission_dir / "treeflow.bbl"
    params:
        submission_dir = str(submission_dir),
        tex = lambda _, input: pathlib.Path(input.tex).relative_to(submission_dir),
        aux = lambda _, output: pathlib.Path(output.aux).relative_to(submission_dir)
    shell:
        """
        cd {params.submission_dir}
        pdflatex {params.tex}
        bibtex {params.aux}
        """

rule treeflow_submission_zip:
    input:
        tex = submission_dir / "treeflow.tex",
        minted_cache = rules.treeflow_submission_minted_cache.output,
        # minted_sty = rules.treeflow_submission_minted_sty.output,
        bbl = rules.treeflow_submission_bbl.output.bbl,
    output:
        zip = manuscript_dir / "out" / "submission.zip"
    params:
        submission_dir = str(submission_dir),
        minted_cache_prefix = str(pathlib.Path(rules.treeflow_submission_minted_cache.input.tex).stem),
        output = lambda _, output: str(pathlib.Path(output.zip).resolve())
    shell:
        """
        cd {params.submission_dir}
        zip -r {params.output} . -x {params.minted_cache_prefix}*
        """

rule compile_ms:
    input:
        main = manuscript_dir / "out" / "{manuscript}.tex",
        bib = manuscript_dir / "tex" / "main.bib",
        #bst = manuscript_dir / "tex" / "plos2015.bst",
    output:
        manuscript_dir / "out" / "{manuscript}.pdf"
    params:
        output_dir =  lambda _, output: pathlib.Path(output[0]).parents[0],
        aux_file = lambda _, output: pathlib.Path(output[0]).with_suffix(".aux"),
        bib_dir = lambda _, input: pathlib.Path(input.bib).parents[0],
        #bst_dir = lambda _, input: pathlib.Path(input.bst).parents[0],
        tex_inputs = manuscript_dir / "tex"
    shell:
        """
        export TEXINPUTS=.:{params.tex_inputs}:
        export BIBINPUTS={params.bib_dir}
        pdflatex --shell-escape -output-directory={params.output_dir} {input.main}
        bibtex {params.aux_file}
        pdflatex --shell-escape -output-directory={params.output_dir} {input.main}
        pdflatex --shell-escape -output-directory={params.output_dir} {input.main}
        """

# export BSTINPUTS=.:{params.bst_dir}: