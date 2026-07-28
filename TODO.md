# TreeFlow paper TODO

Handoff for continuing the Syst Biol revision (USYB-2024-141) on a machine that
can actually run TreeFlow. Written 2026-07-27 on an Intel Mac, which cannot
(see [Environment](#environment)).

## Where the work lives

Two repos, cloned as **sibling directories** (`workflow/ms.smk` resolves the
benchmark data and carnivores example outputs through the `treeflow` package
location, so this layout is required):

```
<parent>/treeflow-paper
<parent>/treeflow
```

Christiaan's branches, stacked:

| Repo | PR | Branch | Base |
|---|---|---|---|
| treeflow-paper | [#15](https://github.com/christiaanjs/treeflow-paper/pull/15) | `claude/treeflow-paper-review-l5dlxb` | `master` |
| treeflow-paper | [#17](https://github.com/christiaanjs/treeflow-paper/pull/17) | `claude/treeflow-benchmark-inline-21lm8p` | above |
| treeflow | [#113](https://github.com/christiaanjs/treeflow/pull/113) | `claude/treeflow-benchmark-inline-21lm8p` | `claude/treeflow-paper-review-l5dlxb` |

Local branch `revisions-fullrank` in **both** repos sits on top of PR #17's head.
In treeflow-paper it also has `origin/master` merged in, which restores
`out/h3n2/{beast.log,beast.trees,timing-data.csv}` — PR #17 predates the commit
that added them, so a two-dot diff makes them look deleted. They are not.

> **Not yet pushed.** Commit `2f8e465` (treeflow-paper) and the `full_rank.py`
> fix (treeflow) exist only on the Intel Mac. Push both before switching
> machines, or redo the work.

## Environment

TreeFlow's `requirements.txt` pins `tensorflow==2.20.0`. There are **no macOS
x86-64 TF wheels past 2.16.2**, and the 2.16.2 wheel is built without AVX2/FMA.
Measured on the Intel Mac: carnivores VI at **2.4 it/s**, against **62–70 it/s**
recorded in the committed `examples/carnivores.ipynb` outputs from Christiaan's
machine — roughly 26x slower. That is what blocks the H3N2 work below. An M4
should meet or beat Christiaan's numbers.

Setup:

```bash
cd treeflow       && python3 -m venv .venv && .venv/bin/pip install -r requirements.txt && .venv/bin/pip install -e .
cd treeflow-paper && pip install -e . && pip install -e ../treeflow && pip install arviz
```

The manuscript pipeline additionally needs R (with the CRAN/Bioconductor
packages listed in `README.md`), a LaTeX distribution, and BEAST 2 >= 2.7 with
Feast. Running the example notebooks needs `nbformat nbclient ipykernel
matplotlib pandas tqdm` on top of `requirements.txt`.

Confirm the environment is good before anything else:

```bash
cd treeflow/examples
TREEFLOW_EXAMPLE_NUM_STEPS=100 TREEFLOW_EXAMPLE_N_RUNS=1 \
  ../.venv/bin/python run_example.py carnivores --timeout 600
```

## Done (commit `2f8e465`, treeflow-paper)

- [x] **Variational inference / Model approximations subsections** rewritten
      around the full-rank Gaussian family that `treeflow_vi` now uses by
      default, with mean field as the contrast case. Covers the trainable
      affine `Lz + mu` applied *before* the split into model variables (which is
      why clock-rate/node-height correlations are representable), its quadratic
      parameter cost, and a recast IAF paragraph — the old "scales
      quadratically" justification no longer distinguished IAF from full rank.
- [x] **Carnivores interpretation** rewritten for the full-rank results:
      locations agree throughout; kappa, site rate shape, birth rate and all
      four base frequencies now match BEAST 2 in spread (the mean-field
      T-frequency artifact is gone); tree height and tree length are over- not
      under-dispersed.
- [x] **Dispersion ratios computed, not hardcoded.**
      `get_carnivores_dispersion_vars` in `treeflow_pipeline/manuscript.py`,
      wired through `workflow/ms.smk`, populates
      `\VAR{carnivores_tree_height_sd_ratio}` and
      `\VAR{carnivores_tree_length_sd_ratio}` from the samples. Cross-check when
      the pipeline first runs: tree height should come out near **1.8**
      (BEAST sd 0.01099 vs VI sd 0.01987, means 0.5022 vs 0.5036), derived from
      the `tree_height` histogram in notebook cell 16 and the BEAST log in
      `supplementary-data/carnivores/beast-results.zip`.
- [x] **Fig. 1 architecture diagram** tidied and rebuilt
      (`manuscript/out/architecture.pdf`). Ronquist's Developer API /
      TensorArray overlap was already fixed; the remaining problems were the
      native-op node stretching the figure to 2:1 with a dead region, and the
      TFP dependency hanging off the left with crossing arrows.

## To do

### 1. H3N2 rerun with the full-rank approximation

Decided: full rank throughout, for consistency with carnivores and with the CLI
default. The likelihood dominates at 980 taxa, so the speed cost should be
small — but measure it rather than asserting it, since the manuscript now
claims exactly that (`treeflow.j2.tex`, Model approximations, "In the analyses
presented here this cost is small relative to evaluating the phylogenetic
likelihood and its gradient").

`workflow/data.smk:177` (`rule variational_fit`) passes **no `-va` flag**, so it
already picks up the `full_rank` default. What needs changing:

- [ ] `-n 30000` is hardcoded in the rule. The last submission converged at
      22,000 of 30,000 iterations, which Christiaan thinks was short for some
      situations. Raise it, and consider lifting it into `config`.
- [ ] Single `config[seed]`. Christiaan wants multiple runs; mirror what
      `examples/carnivores.ipynb` does (loop over seeds, annotate samples with a
      `run` column) so the H3N2 marginals figure can carry the same inter-run
      Monte Carlo band as carnivores.
- [ ] Rerun, then update `out/h3n2/timing-data.csv` and the interpretation.

### 2. H3N2 prose — all still says "mean field"

Blocked on item 1; deliberately left untouched rather than rewritten
speculatively, since the rerun decides which way the uncertainty error goes.
In `manuscript/tex/treeflow.j2.tex`:

- [ ] **L220** — Fig. 4 caption: "underestimated by the mean-field variational
      approximation".
- [ ] **L231** — tree height / pop size "substantially underestimated";
      T-frequency passage and its simplex explanation, which the carnivores
      full-rank result has already falsified.
- [ ] **L233** — "TreeFlow's mean field variational approximation assumes ...
      independent Normal distributions".
- [ ] **L285** — Discussion: "would be ignored by variational inference using a
      mean-field posterior approximation", and the suggestion to adopt
      block-diagonal covariance, which full rank supersedes.

### 3. Verify the speed comparison is apples-to-apples

Ronquist's L440-444 comment turns entirely on this, and the numbers are
awkward: BEAST 22,899 s (3e7 iterations, min ESS 103) against VI 37,844 s to
convergence — VI is currently ~1.65x *slower*.

- [ ] Confirm the H3N2 BEAST and VI timings in `out/h3n2/timing-data.csv` came
      from the same machine under comparable load. The carnivores BEAST screen
      output in `beast-results.zip` carries the warning *"concurrent with other
      workloads — timing not valid for manuscript"*; if the H3N2 timings have
      the same provenance, the comparison needs redoing.
- [ ] Consider whether the native C++ likelihood op should be used for the H3N2
      VI run. At 1.72 s/iteration the traversal likelihood dominates completely,
      so this could change the timing story rather than just improve it. Note
      `treeflow_vi` has no `--native` flag; `use_native="auto"` picks it up if
      it has been built.

### 4. Figure captions and page fit

- [ ] Christiaan's item: captions grew during the revision and the figures no
      longer sit well on a page. Needs a LaTeX build to judge.

### 5. Response letter

- [ ] Write it. `manuscript/reviews/response-to-decision-letter.txt` has
      point-by-point responses from the earlier pass; the full-rank switch
      changes the answers to at least L412-413 (uncertainty summary) and
      L438-444 (convergence and speed).
- [ ] `manuscript/reviews/revision-checklist.md` Section 2 lists the items that
      needed the maintainer's toolchain. Fig. 1 is now done; the rest (Fig. 3/4/5
      re-render, `treeflow.pdf` / `supplementary.pdf` recompile, Sphinx docs
      build) still need doing and can be ticked off on the M4.

## Loose ends

- **`treeflow`, uncommitted:** `treeflow/model/approximation/full_rank.py`
  collected variables via `distribution.trainable_variables`, which walks the
  whole bijector chain and dies with `Error processing property
  '_bijectors_trackable'` on TF 2.16 / Python 3.12. Taking `loc_var` / `raw_var`
  directly is equivalent and version-proof. May not reproduce on TF 2.20, but
  the fix is worth keeping. **`mean_field.py:158` has the identical pattern**
  and will fail the same way — fix both.
- **Manuscript never built here.** No snakemake/jinja2 environment on the Intel
  Mac and `demo-out/` is not committed, so the LaTeX edits above are unverified
  beyond confirming every `\VAR{}` name resolves to something
  `get_treeflow_manuscript_vars` provides. Build early on the M4.
- **`examples/demo-out/` is not committed** — the carnivores samples CSV, tree
  samples and marginal-likelihood YAML that `workflow/ms.smk` consumes are
  regenerated by running the notebook (~55 min on Christiaan's machine for the
  base model's 4 runs, plus ~18 min for the alt model).
- Marginal likelihoods from the last notebook run, for reference: base
  **-193607.64**, lineage variation **-192436.75**.
