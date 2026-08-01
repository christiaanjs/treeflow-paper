"""Accessory figure: parameter-trace convergence for each of the H3N2 multi-run
VI seeds, one subplot per run.

Each run's trace was recorded with --max-trace-coords (see
flu_variational_fit_run in workflow/data.smk), so the full_rank_loc/
full_rank_scale_raw variables already hold a fixed random sample of
coordinates rather than their full (D and D x D) values.
"""
import pathlib
import pickle
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from treeflow.vi.plotting import plot_parameter_traces

trace_dir = pathlib.Path("out/h3n2/variational-multi-run")
seeds = [1, 2, 3, 4]
output_path = pathlib.Path("manuscript/figures/h3n2-multi-run-traces.png")

fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

for seed, ax in zip(seeds, axes.ravel()):
    with open(trace_dir / f"trace-run{seed}.pickle", "rb") as f:
        trace = pickle.load(f)
    plot_parameter_traces(
        trace.parameters,
        sample=True,
        ax=ax,
        coords_per_var=10,
        title=f"Run {seed} (seed {seed})",
    )

fig.suptitle("H3N2 multi-run VI: parameter-trace convergence")
fig.tight_layout()
output_path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output_path, dpi=150)
print(f"Saved to {output_path}")
