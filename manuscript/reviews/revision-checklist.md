# Revision checklist — Ronquist final review (2026-05-29)

Tracks each point in `final-review.txt` and the software/manuscript changes made
in response. Point-by-point wording for the reviewer is in
`response-to-decision-letter.txt`.

## Section 1 — Completed (source edited in this environment)

Manuscript (`manuscript/tex/treeflow.j2.tex` unless noted):

- [x] **L14** — citation order made consistent (ascending year) across `\cite{}` groups.
- [x] **L74–78** — "In other work…" sentence rewritten and split in two.
- [x] **L95–99** — Koptagel et al. (2022, VaiPhy) cited; entry added to `manuscript/tex/main.bib`.
- [x] **L120ff** — added pointer to online manual, YAML format reference and the tutorial as a starting point.
- [x] **Fig. 1** — `manuscript/tex/architecture.tex`: moved "Developer API" label to box corner and raised the `tensorflow.TensorArray` node so they no longer overlap.
- [x] **L212–222** — added an intuitive explanation of the quadratic cost (immutable Tensors ⇒ N copies of size N).
- [x] **Fig. 2 caption** — specified grey band = BEAST 2 Monte Carlo error (bootstrap); noted TreeFlow VI error negligible.
- [x] **L362** — clarified "every lineage" = every branch (one independent kappa per branch).
- [x] **L372–374** — added explanation of the importance-sampling → marginal-likelihood step.
- [x] **Fig. 3 caption + model** — rewrote captions; clarified independent per-branch kappa model and defined "Age" (lineage height in substitutions, not branch length).
- [x] **L383–386** — corrected: uncertainty grows deeper in the tree for *age* estimates; kappa uncertainty shows the opposite pattern.
- [x] **Fig. 4 caption** — moved in-graph text to caption (method comparison, grey band, mean/SD panels).
- [x] **L412–413** — revised uncertainty summary (tree height & pop size substantially, clock rate slightly); added T-frequency observation + intuition.
- [x] **L438–440** — expanded VI convergence-monitoring instructions; pointed to new docs page.
- [x] **L440–444** — reframed the speed conclusion (removed "substantial speed improvements"; justified via BEAST preliminary tuning runs and scaling).
- [x] **Fig. 5** — `scripts/improved-benchmark-plot.R` rewritten to facet by task with software as lines; caption updated.
- [x] **L533–534** — revised scaling conclusion re Table 1 slopes (1.36 vs 1.25) with justification.

Software (`treeflow` repo, branch `claude/treeflow-paper-review-l5dlxb`):

- [x] `docs/source/model-definition.md` — completed reference for all tree/clock/site/substitution options and priors; added per-branch kappa note.
- [x] `docs/source/convergence.md` — new VI convergence-monitoring guide; added to `docs/source/index.rst` toctree.
- [x] `docs/source/installation.md` — added a "Where to start" section.
- [x] `examples/README.md` — fixed `rates-and-dates-model.yaml` link (and a typo).

Deliverables:

- [x] `manuscript/reviews/response-to-decision-letter.txt` — plain-text, ScholarOne-compatible, review copied verbatim with responses inserted after each point.
- [x] This checklist.

## Section 2 — Cannot be fully resolved in this environment (needs maintainer action)

R and LaTeX are not installed in the session used to make these edits, so the
following require running the maintainer's toolchain:

- [ ] **Re-render Fig. 3, Fig. 4, Fig. 5 PNGs** — needs R (`snakemake -s workflow/ms.smk`). Fig. 5's script was restructured; the PNG must be regenerated to reflect the new layout. Fig. 3/4 changes are caption-only, so their PNGs need re-rendering only if the plotted axis labels are also changed.
- [ ] **Recompile `treeflow.pdf` / `supplementary.pdf`** — needs LaTeX; confirm new line numbers, that captions fit, and there is no page overrun.
- [ ] **Rebuild `architecture.pdf`** — needs LaTeX; visually confirm the Fig. 1 label no longer overlaps the TensorArray node (the fix was made blind to rendering).
- [ ] **Build the Sphinx docs** (`cd treeflow/docs && make html`) — confirm the new `convergence.md` toctree entry and cross-references resolve.
- [ ] **Author sign-off** on the reframed scientific wordings (speed comparison, JC-gradient scaling, uncertainty summary, T-frequency intuition, kappa-uncertainty correction) — these change interpretation and should be reviewed by the authors before submission.
