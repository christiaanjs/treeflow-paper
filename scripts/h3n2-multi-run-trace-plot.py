"""Accessory figure: parameter-trace convergence for each of the H3N2 multi-run
VI seeds, one subplot per run.

Each run's trace was recorded with --max-trace-coords (see
multi_run_variational_fit in workflow/data.smk), so the variational location
and scale variables already hold a fixed random sample of coordinates rather
than their full (D and D x D) values.

Run through the workflow (`snakemake -s workflow/ms.smk
manuscript/figures/h3n2-multi-run-traces.png`) rather than directly, so the
approximation and seeds match the runs the figure is meant to describe.
"""
import pathlib
import pickle
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from treeflow.vi.plotting import plot_parameter_traces

trace_paths = list(snakemake.input.traces)
seeds = list(snakemake.params.seeds)
approx = snakemake.params.approx
output_path = pathlib.Path(snakemake.output[0])

n_cols = 2
n_rows = -(-len(seeds) // n_cols)
fig, axes = plt.subplots(
    n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows), sharex=True, squeeze=False
)

for seed, trace_path, ax in zip(seeds, trace_paths, axes.ravel()):
    with open(trace_path, "rb") as f:
        trace = pickle.load(f)
    plot_parameter_traces(
        trace.parameters,
        sample=True,
        ax=ax,
        coords_per_var=10,
        title=f"Run {seed} (seed {seed})",
    )

for ax in axes.ravel()[len(seeds):]:
    ax.set_visible(False)

fig.suptitle(f"H3N2 multi-run VI ({approx}): parameter-trace convergence")
fig.tight_layout()
output_path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output_path, dpi=150)
print(f"Saved to {output_path}")
