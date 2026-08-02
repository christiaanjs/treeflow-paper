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
- [x] **Fig. 2 caption** — specified band = BEAST 2 Monte Carlo error (bootstrap); added a TreeFlow VI band from four independent variational runs, with the curve showing the pooled posterior.
- [x] **L362** — clarified "every lineage" = every branch (one independent kappa per branch).
- [x] **L372–374** — added explanation of the importance-sampling → marginal-likelihood step.
- [x] **Fig. 3 caption + model** — rewrote captions; clarified independent per-branch kappa model. The reviewer's suspicion about the x-axis was correct: `scripts/carnivores-kappa-plot.R` computed `parentHeight - height` (the branch length) while labelling it "Age". Fixed to the branch midpoint height `(parentHeight + height) / 2` and regenerated. Caption also corrected to posterior *mean* (not median) and notes the log scale.
- [x] **Fig. 3(b)** — `scripts/carnivores-tree-plot.R` summarised each node's *parent* height, so every internal node was plotted twice; now summarises each internal node's own height (61 per model). Caption names the axes instead of relying on "against".
- [x] **L383–386** — the reviewer was right about the figure as it stood, but the axis fix reverses the conclusion: with age on the x-axis, kappa uncertainty is *greatest* for the oldest branches, while remaining substantial on the youngest. Passage rewritten to describe both, which fits Duchene et al. (2015) better than either earlier version.
- [x] **Fig. 4 caption** — moved in-graph text to caption (method comparison, grey band, mean/SD panels).
- [x] **L412–413** — the T-frequency intuition was tested and confirmed. New `root_full_rank` variational family places the base-frequency coordinates and the root height in a covariance block; it resolves the T frequency on both datasets while leaving the other three unchanged, and the clock rate is now well matched. Tree height and pop size remain underestimated, attributed to skew plus the mode-seeking KL objective rather than missing correlation. All main-text analyses switched to this family; Supplementary Appendix Figures S3/S4 and Tables S1/S2 compare all three.
- [x] **L438–440** — expanded VI convergence-monitoring instructions; pointed to new docs page.
- [x] **L440–444** — speed conclusion restored and now supported: native accelerated modules make VI converge around eleven times faster than the BEAST 2 run, with the full optimization taking about a quarter of the MCMC time. Also made explicit that BEAST 2 was run with BEAGLE, so both sides use their accelerated likelihoods.
- [x] **Fig. 5** — `scripts/improved-benchmark-plot.R` rewritten to facet by task with software as lines; caption updated. Note the figure now has five series, not three: `TreeFlow (native)` and `JAX (JIT)` were added.
- [x] **L533–534** — benchmarks re-run on a new harness, so the reviewer's quoted slopes (1.36 vs 1.25) no longer appear; they are now 1.48 vs 1.27 and the gap is wider. Rather than call it noise, the conclusion is now anchored to a measured quadratic reference: JIT-compiled JAX has slopes above 2 throughout, so TreeFlow and bito/BEAGLE are both firmly sub-quadratic, and TreeFlow's native op has a shallower slope than bito/BEAGLE. Wording throughout changed from "near-linear"/"close to one" to "sub-quadratic".

Software (`treeflow` repo, branch `claude/treeflow-paper-review-l5dlxb`):

- [x] `docs/source/model-definition.md` — completed reference for all tree/clock/site/substitution options and priors; added per-branch kappa note.
- [x] `docs/source/convergence.md` — new VI convergence-monitoring guide; added to `docs/source/index.rst` toctree.
- [x] `docs/source/installation.md` — added a "Where to start" section.
- [x] `examples/README.md` — fixed `rates-and-dates-model.yaml` link (and a typo).

Deliverables:

- [x] `manuscript/reviews/response-to-decision-letter.txt` — plain-text, ScholarOne-compatible, review copied verbatim with responses inserted after each point.
- [x] This checklist.

## Section 2 — Additions not prompted by a review comment

Disclosed to the reviewer under "ADDITIONS SINCE THE REVIEWED VERSION" in
`response-to-decision-letter.txt`:

- [x] **Native accelerated modules** — compiled C++ TensorFlow custom operations for the pruning likelihood and the node height ratio transform, each with a hand-derived analytic gradient. Reported in the Benchmarks section as `TreeFlow (native)`; also what makes the H3N2 variational timings favourable.
- [x] **Benchmark harness rewritten** — per-evaluation timings (minimum of 100 repeats, 10 replicate datasets), a `JAX (JIT)` series, and simulation moved from BEAST 2 to TreeFlow itself (constant-size coalescent, strict clock, JC, no site rate variation).
- [x] **`root_full_rank` variational family** — see the L412–413 entry above.

## Section 3 — Remaining

- [x] **Re-render figures** — all seven regenerated via `snakemake -s workflow/ms.smk`.
- [x] **Recompile `treeflow.pdf` / `supplementary.pdf` / `treeflow-diff.pdf`**.
- [x] **Rebuild `architecture.pdf`** — Fig. 1 label no longer overlaps the TensorArray node.
- [ ] **Build the Sphinx docs** (`cd treeflow/docs && make html`) — confirm the new `convergence.md` toctree entry and cross-references resolve.
- [ ] **Author sign-off** on the reframed scientific wordings (speed comparison, JC-gradient scaling, uncertainty summary, T-frequency result, kappa-uncertainty reversal, carnivores over-/under-dispersion) — these change interpretation and should be reviewed by the authors before submission.
