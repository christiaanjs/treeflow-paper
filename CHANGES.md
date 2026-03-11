# Change Log for Revision 2

Tracking all changes made for the response to the editor and reviewers.

## Text corrections (R3 detailed comments)

- L65: "of discrete part" → "of the discrete part"
- L82: Added "automatically" before "generate" re LinguaPhylo description
- L82: Fixed "user interace" → "user interface" (additional typo found)
- L92-94: "phylogenetic inference" → "phylogenetic problems"
- L98: Fixed BibTeX rendering of "Matsen IV" by bracing `{IV}` in all 5 bib entries
- L100: "represention" → "representation"
- L109: "fixed-topology phylogenetic inference" → "parameter inference in fixed-topology models of a standard form"
- L153: "probabilsitic" → "probabilistic"
- L161: "on a marginalizing" → "on marginalizing"
- L189: "might be need" → "might need"
- L474ff: Wrote out numbers under ten ("4 distinct" → "four distinct", "4 rate categories" → "four rate categories", "6 GTR" → "six GTR", "4 base" → "four base")
- L562: "it's" → "its"

## References added (R3)

- L13-14: Replaced RAxML with PhyloBayes (Lartillot et al. 2009) in opening paragraph; reframed sentence to stress Bayesian MCMC focus
- L19-20: Added Lakner et al. (2008, Syst Bio) to MCMC efficiency citations
- L28: Added Höhna et al. (2014, Syst Bio) on probabilistic graphical models for phylogenetics
- L78: Added Senderov et al. (2024, TreePPL) alongside Ronquist et al. (2021)
- L369: Added Burda et al. (2015) for importance sampling correcting VI approximation error

## References updated to published versions

- Drummond et al. LinguaPhylo: bioRxiv 2022 → PLOS Computational Biology 2023
- Ji et al. ratio transformations: arXiv 2021 → Systematic Biology 72(5), 2023
- Fisher et al. shrinkage clocks: arXiv 2021 → Molecular Biology and Evolution 40(11), 2023
- Senderov et al. TreePPL: updated to 2024 v2 with correct title and authors

## Structural reorganization (R3)

- [x] Reorder Software Description: new "User interface" subsection (CLI + model definition) now leads, followed by "Developer API" subsection introducing the library internals, then the existing technical subsections
- [x] Reorder Discussion: moved comparison with RevBayes/Blang (TreeFlow advantages) to lead the Discussion; limitations follow with forward-looking additions
- [x] Refocus example conclusions: added TreeFlow-focused concluding sentences to both carnivores (ease of model extension) and influenza (scalability, YAML-only specification) examples
- [x] Moved Figs 3 (code listing) & 5 (YAML model spec) to new supplementary appendix (`manuscript/tex/supplementary.tex`) with introductory text on model specification philosophy; main text now references Supplementary Appendix Figures S1 and S2

## Figure modifications (R3)

- [x] Fig 1 (architecture.tex): added distinct fill colors (light blue for User Interface, beige for Developer API), thick borders, bold labels
- [x] Fig 2 (marginals): regenerated with bootstrap density bands for BEAST MC error; fixed snakemake R input accessor (`[` → `[[`) in `improved-marginals-plot.R`
- [ ] Fig 6a (kappa): add batch-means MC SE error bars — R script updated (`carnivores-kappa-plot.R`); awaiting `carnivores.ipynb` notebook execution to produce `carnivores-alt-trees.nexus`
- [x] Fig 7 (benchmark): regenerated with reversed legend symbol order

## Clarifications & expansions (R3)

- [x] L236-239: simplified bijector chaining explanation with concrete example (real values → ratios → node heights → full tree)
- [x] L274-276: simplified sentence about topology inference support
- [x] L369: added intuitive explanation of why importance sampling corrects VI error (importance weights account for discrepancy between approximation and true posterior)
- [x] L369-371: removed redundant marginal likelihood numbers from first paragraph, moved to second paragraph with interpretation
- [x] L428: clarified BEAST used built-in autotuning plus manual operator weight adjustment; added ELBO monitoring as VI convergence diagnostic
- [x] L529-538: expanded structured approximations with concrete examples (block-diagonal covariance for clock rate/tree height correlation, normalizing flows)
- [x] L539-545: added forward-looking direction for topology inference (subsplit Bayesian networks combined with TreeFlow's existing components)

## Manuscript version specificity

- Specified BEAST version as 2.7.7 at first analysis mention (carnivores example)

## Pipeline compatibility fixes (treeflow-paper repo)

- Updated all BEAST XML templates to 2.7 namespaces (`beast.core.*` → `beast.base.*`, `beast.evolution.*` → `beast.base.evolution.*`, `beast.app.seqgen` → `beastfx.app.seqgen`): `beast-analysis.j2.xml`, `sim-seq.j2.xml`, `tree-sim.j2.xml`, `rate-sim.j2.xml`
- Updated fully qualified class names in `templating.py` (`YuleModel`, `WeibullSiteModel`)
- Updated `supplementary-data/*/beast-2.7.xml` version attributes from `required="BEAST v2.5.2"` to `required="BEAST v2.7.0"`
- Fixed `data.smk`: filtered starting values to only include free parameters in `variational_fit` and `ml_fit` rules (prevents `Unknown parameters in initial values: {'clock_rate'}` for strict clock models)
- Fixed `setup.py`: removed `importlib` dependency (Python 2 backport, not needed)
- Fixed `results.py`: arviz 1.0 compatibility — `compute_beast_ess` now passes `(chain, draw)` shaped arrays via `np.expand_dims`
- Fixed `treeflow.j2.tex`: added `backgrounds` tikz library and `apibg`/`uibg` color definitions so `\includestandalone{architecture}` works in main document

## Documentation

- Wrote comprehensive README for treeflow-paper repo (reproduction guide, requirements, workflows)
- Wrote README for treeflow-benchmarks repo (dependencies, usage, configuration, outputs)
- Wrote README for phylojax repo (features, installation, usage example)

## Software / installation fixes (R1 & R2) — in ~/Git/treeflow repo

- [x] Fix install docs: removed `--user` from venv instructions in `docs/source/installation.md`
- [x] Add `tf-keras` as explicit dependency in `setup.cfg` and `requirements.txt`
- [x] Fix TF/TFP version incompatibility: replaced private `nest._get_attrs_items` API with `attr.fields` in `treeflow/tf_util/attrs.py`
- [x] Fix Docker build: updated `silence-tensorflow` from 1.2.1 to 1.2.3; fixed `FROM...as` casing; fixed CMD to JSON format
- [x] Pin `numpy>=1.19,<2.1` in `setup.cfg` to avoid TFP incompatibility with numpy 2.1+ (`reshape(newshape=)` removed)
- [x] Fix Keras 3 compatibility: removed `name=` kwarg from `apply_gradients()` in `treeflow/vi/optimizers/robust_optimizer.py`
- [x] Removed unused `numpy.core.fromnumeric` import causing deprecation warning in `treeflow/tree/topology/tensorflow_tree_topology.py`
- [x] Removed debug "Temporary hotfix" warning from `attrs.py`
- [x] Verified end-to-end: `pip install .` + `treeflow_vi` CLI runs cleanly with Python 3.12, TF 2.20, TFP 0.25, numpy 2.0.2
- [x] Test suite: 255/260 pass (4 failures: 1 optional `bito` dep, 1 float precision, 2 incomplete cascading flows)
- [x] Updated installation docs: Python 3.9+ (3.12 recommended), fixed Docker flag order
- [x] CI docs build: updated Python 3.8→3.12, actions/setup-python v4→v5
- [x] Fixed Sphinx docs: replaced broken `sphinxcontrib.napoleon` with built-in `sphinx.ext.napoleon`
- [x] Fixed dev/requirements.txt: resolved Sphinx/docutils/nbsphinx-link version conflicts (Sphinx 7.4.7, docutils 0.20.1, myst-parser 4.0.1)
- [x] All changes in treeflow PR #75 (`fix/reviewer-installation-issues` branch), CI passing
