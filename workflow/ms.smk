import pathlib
import shutil
import sys
import pandas as pd
import treeflow_pipeline.model
from treeflow_pipeline.util import yaml_input, text_input, text_output
import treeflow_pipeline.manuscript
import treeflow_pipeline.diff
import treeflow_pipeline.submission

import treeflow
treeflow_dir = pathlib.Path(treeflow.__file__).parents[1]
# The benchmark was inlined into the treeflow repo (experiments/benchmarks); its
# notebook writes the manuscript-schema CSVs and config consumed below. This
# replaces the previously separate ``treeflow_benchmarks`` package.
treeflow_benchmark_data_dir = treeflow_dir / "experiments" / "benchmarks" / "data"

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
# Everything the journal asks to be uploaded for the final submission, gathered
# into one directory (see the final_submission_* rules at the end of this file).
final_submission_dir = manuscript_dir / "out" / "final-submission"
# The .bib file is renamed in the bundle, so the submitted main text refers to
# it by this stem rather than the working copy's "main".
submission_bibliography = "treeflow"
minted_cache_dir = "minted-cache"
dataset_dir = "{dataset}"
supplementary_data_dir = pathlib.Path("supplementary-data")
diff_base_commit = "07911c7c8a71461e9a912dc4baad6066fa27e8d1"

rule ms:
    input:
        #manuscript_dir / "out" / "submission.zip",
        manuscript_dir / "out" / "treeflow.pdf",
        manuscript_dir / "out" / "supplementary.pdf",
        manuscript_dir / "out" / "response-letter.pdf",
        #supplementary_data_dir / config["flu_dataset"] / "beast.xml",
        #supplementary_data_dir / "carnivores" / "beast.xml",
        manuscript_dir / "out" / "treeflow-diff.pdf",
        manuscript_dir / "out" / "final-submission.zip"

rule ms_diff:
    input:
        manuscript_dir / "out" / "treeflow-diff.pdf"

rule extract_old_treeflow_tex:
    output: temp(manuscript_dir / "out" / "treeflow-old.tex")
    params:
        commit = diff_base_commit,
        tex_path = "manuscript/out/treeflow.tex"
    shell:
        "git show {params.commit}:{params.tex_path} > {output}"

rule latexdiff_treeflow:
    input:
        old = manuscript_dir / "out" / "treeflow-old.tex",
        new = manuscript_dir / "out" / "treeflow.tex"
    output:
        manuscript_dir / "out" / "treeflow-diff.tex"
    run:
        text_output(treeflow_pipeline.diff.run_latexdiff(input.old, input.new), output[0])

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
        plot_data = treeflow_benchmark_data_dir / "manuscript-plot-data.csv",
        fit_table = treeflow_benchmark_data_dir / "manuscript-fit-table.csv"
    output:
        tex = manuscript_dir / "tables" / "benchmark-table.tex"
    run:
        treeflow_pipeline.manuscript.benchmark_summary_table(input.plot_data, input.fit_table, output.tex)

rule benchmark_plot:
    input:
        plot_data = treeflow_benchmark_data_dir / "manuscript-plot-data.csv",
    output:
        plot = manuscript_dir / "figures" / "benchmark-log-scale-plot.png"
    params:
        python_executable = sys.executable
    script:
        "../scripts/improved-benchmark-plot.R"

# The same plot rendered to PDF for the separate figure file upload. Systematic
# Biology requires figures as tif/eps/pdf/gif/jpg at a minimum of 600 dpi; the
# PNG above is written at ggplot2's default 300 dpi, so the benchmark figure is
# submitted as vector PDF, which has no resolution limit. Kept as a separate
# output rather than changing benchmark_plot's extension so the manuscript build
# is unaffected.
rule benchmark_plot_pdf:
    input:
        plot_data = treeflow_benchmark_data_dir / "manuscript-plot-data.csv",
    output:
        plot = manuscript_dir / "out" / "benchmark-log-scale-plot.pdf"
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

# The main-text marginals figures compare BEAST 2 against TreeFlow VI fit with
# the root_full_rank approximation, run 4 times independently (rule
# multi_run_variational_fit in workflow/data.smk, matching the settings the
# carnivores example notebook uses: num_steps 60,000, learning rate 0.001) so a
# Monte Carlo error band can be shown. The comparison against the mean-field and
# full-rank approximations is a supplementary figure
# (approximation_comparison_plot below). This takes precedence over the generic
# data_marginals_plot rule.
MAIN_APPROX = "root_full_rank"

ruleorder: carnivores_marginals_plot > data_marginals_plot

rule carnivores_marginals_plot:
    input:
        vi_samples_root_full_rank = out_dir / "carnivores" / "variational-multi-run" / MAIN_APPROX / "samples.csv",
        beast_samples = out_dir / "carnivores" / "beast.log"
    output:
        manuscript_dir / "figures" / "carnivores-marginals.png"
    params:
        python_executable = sys.executable
    script:
        "../scripts/multi-run-marginals-plot.R"

# The flu marginals figure gets the same multi-run treatment as carnivores, so
# it takes precedence over the generic data_marginals_plot rule for the flu
# dataset too.
ruleorder: flu_marginals_plot > data_marginals_plot

rule flu_marginals_plot:
    input:
        vi_samples_root_full_rank = out_dir / config["flu_dataset"] / "variational-multi-run" / MAIN_APPROX / "samples.csv",
        beast_samples = out_dir / config["flu_dataset"] / "beast.log"
    output:
        manuscript_dir / "figures" / f"{config['flu_dataset']}-marginals.png"
    params:
        python_executable = sys.executable
    script:
        "../scripts/multi-run-marginals-plot.R"

# Supplementary: how the choice of variational approximation family affects the
# fitted posterior. One run (seed 1) of each of mean_field, full_rank and
# root_full_rank, against the same BEAST 2 reference. Single runs, so no
# inter-run band -- this figure is about the systematic differences between the
# families, not their Monte Carlo error.
COMPARISON_APPROXES = ["mean_field", "full_rank", "root_full_rank"]

rule approximation_comparison_plot:
    input:
        vi_samples = expand(
            out_dir / "{{dataset}}" / "variational-multi-run" / "{approx}" / "samples-run1.csv",
            approx=COMPARISON_APPROXES
        ),
        beast_samples = out_dir / "{dataset}" / "beast.log"
    output:
        manuscript_dir / "figures" / "{dataset}-approximation-comparison.png"
    params:
        python_executable = sys.executable,
        approxes = COMPARISON_APPROXES
    script:
        "../scripts/approximation-comparison-plot.R"

# Companion to approximation_comparison_plot: the fitted ELBO and wall-clock
# runtime of the same runs, so the visual comparison can be read alongside the
# quantity being optimised. ELBOs are comparable within a dataset but not
# across datasets -- the log likelihood sums over alignment patterns, of which
# carnivores has far more than H3N2 despite having far fewer taxa.
rule approximation_comparison_table:
    input:
        logs = expand(
            out_dir / "{{dataset}}" / "variational-multi-run" / "{approx}" / "log-run1.txt",
            approx=COMPARISON_APPROXES
        ),
        benchmarks = expand(
            out_dir / "{{dataset}}" / "variational-multi-run" / "{approx}" / "benchmark-run1.txt",
            approx=COMPARISON_APPROXES
        )
    output:
        manuscript_dir / "tables" / "{dataset}-approximation-comparison.tex"
    params:
        approxes = COMPARISON_APPROXES,
        # `treeflow_vi run` prints the *sum* of the last --elbo-samples per-step
        # loss values, so divide by that to recover a per-step ELBO estimate.
        elbo_samples = 100
    run:
        import re
        labels = {
            "mean_field": "Mean field",
            "full_rank": "Full rank",
            "root_full_rank": "Root full rank",
        }
        rows = []
        for approx, log_path, benchmark_path in zip(
            params.approxes, input.logs, input.benchmarks
        ):
            elbo_matches = re.findall(
                r"ELBO estimate: (\S+)", pathlib.Path(log_path).read_text()
            )
            if not elbo_matches:
                raise ValueError(f"No ELBO estimate found in {log_path}")
            elbo = float(elbo_matches[-1]) / params.elbo_samples
            seconds = float(
                pathlib.Path(benchmark_path).read_text().splitlines()[1].split("\t")[0]
            )
            rows.append(
                f"{labels.get(approx, approx)} & {elbo:,.0f} & {seconds / 60:.1f} \\\\"
            )
        text_output(
            "\n".join(
                [
                    r"\begin{tabular}{lrr}",
                    r"\hline",
                    r"Approximation & ELBO & Runtime (minutes) \\",
                    r"\hline",
                ]
                + rows
                + [r"\hline", r"\end{tabular}"]
            ),
            output[0],
        )

# Accessory convergence diagnostic, not a manuscript figure: the variational
# parameter traces for each seed of the main-text H3N2 campaign, so the
# "parameters have stopped drifting" check described in the Scalable inference
# section can be reproduced.
MULTI_RUN_SEEDS = [1, 2, 3, 4]

rule flu_multi_run_trace_plot:
    input:
        traces = expand(
            out_dir / config["flu_dataset"] / "variational-multi-run" / MAIN_APPROX / "trace-run{run}.pickle",
            run=MULTI_RUN_SEEDS
        )
    output:
        manuscript_dir / "figures" / f"{config['flu_dataset']}-multi-run-traces.png"
    params:
        seeds = MULTI_RUN_SEEDS,
        approx = MAIN_APPROX
    script:
        "../scripts/h3n2-multi-run-trace-plot.py"

rule data_tree_plot:
    input:
        vi_tree_samples = out_dir / "{dataset}" / "variational-tree-samples.nexus",
        beast_tree_samples = out_dir / "{dataset}" / "beast.trees"
    output:
        plot = manuscript_dir / "figures" / "{dataset}-trees.png"
    script:
        "../scripts/data-tree-plot.R"

# Analogous to flu_marginals_plot above: use the pooled multi-run tree samples
# (all 4 runs' trees combined into one file), so the per-node height mean/SD
# comparison against BEAST reflects between-run as well as within-run
# variability. Takes precedence over the generic data_tree_plot rule for flu.
ruleorder: flu_tree_plot > data_tree_plot

rule flu_tree_plot:
    input:
        vi_tree_samples_root_full_rank = out_dir / config["flu_dataset"] / "variational-multi-run" / MAIN_APPROX / "tree-samples.nexus",
        beast_tree_samples = out_dir / config["flu_dataset"] / "beast.trees"
    output:
        plot = manuscript_dir / "figures" / f"{config['flu_dataset']}-trees.png"
    script:
        "../scripts/data-tree-plot.R"

# The per-node agreement statistics quoted in the text alongside flu_tree_plot,
# computed from the same inputs so the two cannot disagree.
rule flu_tree_stats:
    input:
        vi_tree_samples = out_dir / config["flu_dataset"] / "variational-multi-run" / MAIN_APPROX / "tree-samples.nexus",
        beast_tree_samples = out_dir / config["flu_dataset"] / "beast.trees"
    output:
        out_dir / config["flu_dataset"] / "tree-comparison-stats.yaml"
    script:
        "../scripts/tree-comparison-stats.R"

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
        treeflow_benchmarks_config = treeflow_benchmark_data_dir / "benchmark-config.yaml",
        benchmark_plot = rules.benchmark_plot.output.plot,
        benchmark_summary_table = rules.benchmark_summary_table.output[0],
        carnivores_marginals_plot = manuscript_dir / "figures" / "carnivores-marginals.png",
        carnivores_kappa_plot = rules.carnivores_kappa_plot.output[0],
        carnivores_tree_plot = rules.carnivores_tree_plot.output[0],
        carnivores_marginal_likelihoods = treeflow_dir / "examples" / "demo-out" / "carnivores-marginal-log-likelihoods.yaml",
        flu_marginals_plot = manuscript_dir / "figures" / f"{config['flu_dataset']}-marginals.png",
        flu_tree_plot = manuscript_dir / "figures" / f"{config['flu_dataset']}-trees.png",
        flu_timing_csv = out_dir / config["flu_dataset"] / "timing-data.csv",
        flu_tree_stats = rules.flu_tree_stats.output[0],
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
                        carnivores_marginal_likelihoods=yaml_input(input.carnivores_marginal_likelihoods),
                        minted_cache_dir=params.output_dir / "minted-cache",
                        bibliography_file=input.bib
                    ),
                    output_dir = params.output_dir,
                    **yaml_input(input.flu_tree_stats),
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
        style = manuscript_dir / "tex" / "sysbio_sse.cls",
    output:
        benchmark_plot = submission_dir / "benchmark-plot.png",
        carnivores_marginals_plot = submission_dir / "carnivores-marginals.png",
        carnivores_kappa_plot = submission_dir / "carnivores-kappa.png",
        carnivores_tree_plot = submission_dir / "carnivores-trees.png",
        flu_marginals_plot = submission_dir / "flu-marginals.png",
        flu_tree_plot = submission_dir / "flu-trees.png",
        flu_model_file = submission_dir / "flu-model.yaml",
        bib = submission_dir / "treeflow.bib",
        style = submission_dir / "sysbio_sse.cls",
    run:
        for key in input.keys():
            shell(f"cp {input[key]} {output[key]}")

rule template_treeflow_submission_ms:
    input:
        template = tex_template,
        body_template = manuscript_dir / "tex" / "treeflow.j2.tex",
        treeflow_benchmarks_config = treeflow_benchmark_data_dir / "benchmark-config.yaml",
        benchmark_plot = rules.treeflow_submission_dir.output.benchmark_plot,
        benchmark_summary_table = rules.benchmark_summary_table.output[0],
        carnivores_marginals_plot = rules.treeflow_submission_dir.output.carnivores_marginals_plot,
        carnivores_kappa_plot = rules.treeflow_submission_dir.output.carnivores_kappa_plot,
        carnivores_tree_plot = rules.treeflow_submission_dir.output.carnivores_tree_plot,
        flu_marginals_plot = rules.treeflow_submission_dir.output.flu_marginals_plot,
        flu_tree_plot = rules.treeflow_submission_dir.output.flu_tree_plot,
        flu_timing_csv = out_dir / config["flu_dataset"] / "timing-data.csv",
        flu_tree_stats = rules.flu_tree_stats.output[0],
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
        zip -r {params.output} . -x {params.minted_cache_prefix}* *.pdf
        """

submission_figures_dir = manuscript_dir / "out" / "figures"

# Ordered mapping of submission figure names to source files.
# Subfigures within the same figure environment use letter suffixes (e.g. 3a, 3b).
# Figure 5 is taken from the vector PDF rendering (see benchmark_plot_pdf) to
# meet the journal's 600 dpi minimum for separately uploaded figure files.
submission_figures = {
    "figure-1": manuscript_dir / "out" / "architecture.pdf",
    "figure-2": manuscript_dir / "figures" / "carnivores-marginals.png",
    "figure-3": manuscript_dir / "out" / "figure-3.pdf",
    "figure-4": manuscript_dir / "out" / "figure-4.pdf",
    "figure-5": rules.benchmark_plot_pdf.output.plot,
}

submission_figure_descriptions = {
    "figure-1": "Figure 1: TreeFlow package architecture",
    "figure-2": "Figure 2: carnivores marginal posteriors",
    "figure-3": "Figure 3: carnivores per-branch kappa (a) and node ages (b)",
    "figure-4": "Figure 4: influenza marginal posteriors (a) and node heights (b)",
    "figure-5": "Figure 5: phylogenetic likelihood benchmark (vector)",
}

rule compile_architecture_figure:
    input: manuscript_dir / "tex" / "architecture.tex"
    output: manuscript_dir / "out" / "architecture.pdf"
    params:
        output_dir = str(manuscript_dir / "out"),
        tex_inputs = manuscript_dir / "tex"
    shell:
        """
        export TEXINPUTS=.:{params.tex_inputs}:
        pdflatex -output-directory={params.output_dir} {input}
        """

rule extract_compound_figure_3:
    input: manuscript_dir / "out" / "treeflow.tex"
    output: manuscript_dir / "out" / "figure-3.tex"
    run:
        text_output(treeflow_pipeline.diff.extract_compound_figure_tex(input[0], 0), output[0])

rule extract_compound_figure_4:
    input: manuscript_dir / "out" / "treeflow.tex"
    output: manuscript_dir / "out" / "figure-4.tex"
    run:
        text_output(treeflow_pipeline.diff.extract_compound_figure_tex(input[0], 1), output[0])

rule compile_compound_figure:
    input: manuscript_dir / "out" / "figure-{n}.tex"
    output: manuscript_dir / "out" / "figure-{n}.pdf"
    params:
        output_dir = str(manuscript_dir / "out"),
        tex_inputs = manuscript_dir / "tex"
    shell:
        """
        export TEXINPUTS=.:{params.tex_inputs}:
        pdflatex -output-directory={params.output_dir} {input}
        """

ruleorder: compile_compound_figure > compile_ms

rule copy_submission_figures:
    input:
        **{name: src for name, src in submission_figures.items()}
    output:
        **{name: submission_figures_dir / (name + pathlib.Path(str(src)).suffix)
           for name, src in submission_figures.items()}
    run:
        # Copy in Python rather than via shell(): a shell() format string is
        # expanded by snakemake, where "{input[key]}" indexes with the literal
        # string "key" rather than with this loop variable's value.
        pathlib.Path(str(submission_figures_dir)).mkdir(parents=True, exist_ok=True)
        for key in input.keys():
            shutil.copy(input[key], output[key])

rule ms_figures:
    input: list(rules.copy_submission_figures.output)

rule extract_benchmark_table:
    input: manuscript_dir / "out" / "treeflow.tex"
    output: manuscript_dir / "out" / "table-1.tex"
    run:
        text_output(treeflow_pipeline.diff.extract_table_tex(input[0], remove_caption=True), output[0])

rule ms_tables:
    input: manuscript_dir / "out" / "table-1.pdf"

rule compile_table_tex:
    input: manuscript_dir / "out" / "table-1.tex"
    output: manuscript_dir / "out" / "table-1.pdf"
    params:
        output_dir = str(manuscript_dir / "out"),
        tex_inputs = manuscript_dir / "tex"
    shell:
        """
        export TEXINPUTS=.:{params.tex_inputs}:
        pdflatex -output-directory={params.output_dir} {input}
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
        export BSTINPUTS={params.bib_dir}
        export BIBINPUTS={params.bib_dir}
        pdflatex --shell-escape -output-directory={params.output_dir} {input.main}
        bibtex {params.aux_file}
        pdflatex --shell-escape -output-directory={params.output_dir} {input.main}
        pdflatex --shell-escape -output-directory={params.output_dir} {input.main}
        """

rule compile_supplementary:
    input:
        tex = manuscript_dir / "tex" / "supplementary.tex",
        approximation_comparison_plots = expand(
            manuscript_dir / "figures" / "{dataset}-approximation-comparison.png",
            dataset=["carnivores", config["flu_dataset"]]
        ),
        approximation_comparison_tables = expand(
            manuscript_dir / "tables" / "{dataset}-approximation-comparison.tex",
            dataset=["carnivores", config["flu_dataset"]]
        )
    output:
        manuscript_dir / "out" / "supplementary.pdf"
    params:
        output_dir = lambda _, output: pathlib.Path(output[0]).parents[0],
        tex_inputs = manuscript_dir / "tex"
    shell:
        """
        export TEXINPUTS=.:{params.tex_inputs}:
        pdflatex --shell-escape -output-directory={params.output_dir} {input.tex}
        pdflatex --shell-escape -output-directory={params.output_dir} {input.tex}
        """

rule compile_response_letter:
    input:
        tex = manuscript_dir / "response-letter.tex"
    output:
        manuscript_dir / "out" / "response-letter.pdf"
    params:
        output_dir = lambda _, output: pathlib.Path(output[0]).parents[0],
        tex_inputs = manuscript_dir / "tex"
    shell:
        """
        export TEXINPUTS=.:{params.tex_inputs}:
        pdflatex --shell-escape -output-directory={params.output_dir} {input.tex}
        """

rule supplementary_data:
    input:
        beast_xml = out_dir / dataset_dir / "beast.xml",
        model_file = out_dir / dataset_dir / "model.yaml",
        topology = out_dir / dataset_dir / "topology.nwk"
    output:
        beast_xml = supplementary_data_dir / dataset_dir / "beast.xml",
        model_file = supplementary_data_dir / dataset_dir / "model.yaml",
        topology = supplementary_data_dir / dataset_dir / "topology.nwk"
    run:
        for key in input.keys():
            shell(f"cp {input[key]} {output[key]}")

# --- Final submission bundle -------------------------------------------------
#
# Systematic Biology asks for, at the final submission stage:
#   * a clean copy of the main text in an editable format, i.e. the LaTeX
#     source together with the compiled PDF and the accompanying files (.bib,
#     the journal class, the bibliography style);
#   * the figures removed from the main text and uploaded as separate files
#     (tables may stay embedded as long as they remain editable, which they are:
#     the table source is inlined into the manuscript, not included as an image);
#   * the figure legends placed after the reference list;
#   * every figure in one of tif/eps/pdf/gif/jpg at no less than 600 dpi.
#
# The rules below assemble exactly those files in one directory, plus a manifest
# describing what each file is, and zip it up for upload.

final_submission_main_text_stem = "treeflow-main-text"

# The main text, taken from the built manuscript (so the submitted source is the
# source the manuscript PDF was compiled from) with the figures lifted out into
# a figure legends section after the reference list.
rule final_submission_main_text:
    input:
        tex = manuscript_dir / "out" / "treeflow.tex"
    output:
        tex = final_submission_dir / f"{final_submission_main_text_stem}.tex"
    params:
        bibliography = submission_bibliography
    run:
        text_output(
            treeflow_pipeline.submission.build_main_text(
                input.tex, bibliography=params.bibliography
            ),
            output.tex,
        )

# The files the main text source needs to compile on its own: the bibliography
# database, the journal class, and the bibliography style.
rule final_submission_support_files:
    input:
        bib = manuscript_dir / "tex" / "main.bib",
        cls = manuscript_dir / "tex" / "sysbio_sse.cls",
        bst = manuscript_dir / "tex" / "CSE.bst",
    output:
        bib = final_submission_dir / f"{submission_bibliography}.bib",
        cls = final_submission_dir / "sysbio_sse.cls",
        bst = final_submission_dir / "CSE.bst",
    run:
        pathlib.Path(str(final_submission_dir)).mkdir(parents=True, exist_ok=True)
        for key in input.keys():
            shutil.copy(input[key], output[key])

# Compiled inside the submission directory, with neither TEXINPUTS/BIBINPUTS nor
# --shell-escape set: this checks that the bundle is self-contained and that it
# builds the way a production system would build it.
rule compile_final_submission_main_text:
    input:
        tex = rules.final_submission_main_text.output.tex,
        bib = rules.final_submission_support_files.output.bib,
        cls = rules.final_submission_support_files.output.cls,
        bst = rules.final_submission_support_files.output.bst,
    output:
        pdf = final_submission_dir / f"{final_submission_main_text_stem}.pdf",
        bbl = final_submission_dir / f"{final_submission_main_text_stem}.bbl",
    params:
        submission_dir = str(final_submission_dir),
        tex = lambda _, input: pathlib.Path(input.tex).name,
        aux = lambda _, input: pathlib.Path(input.tex).with_suffix(".aux").name,
    shell:
        """
        cd {params.submission_dir}
        pdflatex {params.tex}
        bibtex {params.aux}
        pdflatex {params.tex}
        pdflatex {params.tex}
        """

ruleorder: compile_final_submission_main_text > compile_ms

# One file per figure, as uploaded.
rule final_submission_figures:
    input:
        **dict(rules.copy_submission_figures.output.items())
    output:
        **{name: final_submission_dir / pathlib.Path(str(path)).name
           for name, path in rules.copy_submission_figures.output.items()}
    run:
        pathlib.Path(str(final_submission_dir)).mkdir(parents=True, exist_ok=True)
        for key in input.keys():
            shutil.copy(input[key], output[key])

rule final_submission_supplementary:
    input:
        supplementary = manuscript_dir / "out" / "supplementary.pdf"
    output:
        supplementary = final_submission_dir / "supplementary-appendix.pdf"
    run:
        shutil.copy(input.supplementary, output.supplementary)

rule final_submission_manifest:
    input:
        tex = rules.final_submission_main_text.output.tex,
        pdf = rules.compile_final_submission_main_text.output.pdf,
        bbl = rules.compile_final_submission_main_text.output.bbl,
        bib = rules.final_submission_support_files.output.bib,
        cls = rules.final_submission_support_files.output.cls,
        bst = rules.final_submission_support_files.output.bst,
        supplementary = rules.final_submission_supplementary.output.supplementary,
        figures = list(rules.final_submission_figures.output),
    output:
        final_submission_dir / "MANIFEST.txt"
    run:
        text_output(
            treeflow_pipeline.submission.manifest(
                main_text={
                    "main text, LaTeX source (figures removed, legends after the "
                    "reference list, tables embedded and editable)": input.tex,
                    "main text, compiled PDF": input.pdf,
                },
                figures={
                    submission_figure_descriptions[name]: path
                    for name, path in rules.final_submission_figures.output.items()
                },
                supporting={
                    "bibliography database (BibTeX)": input.bib,
                    "compiled bibliography": input.bbl,
                    "journal document class": input.cls,
                    "bibliography style": input.bst,
                    "supplementary appendix": input.supplementary,
                },
            ),
            output[0],
        )

rule final_submission_zip:
    input:
        manifest = rules.final_submission_manifest.output,
        tex = rules.final_submission_main_text.output.tex,
        pdf = rules.compile_final_submission_main_text.output.pdf,
        bbl = rules.compile_final_submission_main_text.output.bbl,
        support = list(rules.final_submission_support_files.output),
        figures = list(rules.final_submission_figures.output),
        supplementary = rules.final_submission_supplementary.output.supplementary,
    output:
        zip = manuscript_dir / "out" / "final-submission.zip"
    params:
        submission_dir = str(final_submission_dir),
        # LaTeX's working files are not part of the submission.
        excluded = "'*.aux' '*.log' '*.blg' '*.out'",
        output = lambda _, output: str(pathlib.Path(output.zip).resolve()),
    shell:
        """
        rm -f {params.output}
        cd {params.submission_dir}
        zip -r {params.output} . -x {params.excluded}
        """

rule ms_submission:
    input: rules.final_submission_zip.output.zip
