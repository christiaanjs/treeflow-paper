# Revision Plan — TreeFlow Manuscript (USYB-2024-141, Round 2)

Decision: Minor Revisions (Editor-in-Chief) / Major Revision (AE recommendation)
Reviews received: 2024-12-31
Reviewers: R1 (anonymous), R2 (anonymous), R3 (Fredrik Ronquist)

## Phase 1: Mechanical text fixes ✅
Unambiguous edits to `manuscript/tex/treeflow.j2.tex`:

- [x] L65: "of discrete part" → "of the discrete part"
- [x] L100: "represention" → "representation"
- [x] L153: "probabilsitic" → "probabilistic"
- [x] L161: "on a marginalizing" → "on marginalizing"
- [x] L189: "might be need" → "might need"
- [x] L562 (Discussion): "it's" → "its"
- [x] L98: "Zhang and IV" → "Zhang and Matsen IV"
- [x] L82: add "automatically" before "generate"
- [x] L92-94: "phylogenetic inference" → "phylogenetic problems"
- [x] L109: rephrase "fixed-topology phylogenetic inference"
- [x] L307/L474ff: write out numbers under ten ("We benchmark 4" → "We benchmark four")

## Phase 2: Missing references ✅
Add to `manuscript/tex/main.bib` and cite in text:

- [x] L13-14: Added PhyloBayes (Lartillot 2009), reframed opening for Bayesian MCMC focus, replaced RAxML with PhyloBayes
- [x] L19-20: Added Lakner et al. (2008, Syst Bio)
- [x] L28: Added Höhna et al. (2014, Syst Bio) on phylogenetic graphical models
- [x] L78: Added Senderov et al. (2023) preprint
- [x] L369: Added Burda et al. (2015) for importance sampling correcting VI approximation error

## Phase 3: Structural reorganization
Needs author input on scope:

- [ ] Reorder Software Description: lead with User Interface (command line + model spec), then Developer API
- [ ] Reorder Discussion: lead with TreeFlow advantages (current L553ff), then limitations
- [ ] Refocus example conclusions on TreeFlow capabilities, not biology (carnivores L384-387, influenza)
- [ ] Move Fig 3 (code listing) and Fig 5 (YAML model spec) to supplementary appendix
- [ ] Create supplementary appendix with simpler introductory examples of model spec language and API

## Phase 4: Figure modifications

- [ ] Fig 1 (architecture.tex): add color/shading to distinguish Developer API vs User Interface; reposition User Interface to top or side
- [ ] Figs 2 & 6a: add Monte Carlo error estimates for BEAST (needs ESS data)
- [ ] Fig 7: reverse legend symbol order (needs plot regeneration)

## Phase 5: Clarifications & expansions

- [ ] L236-239: simplify explanation of bijector chaining and fixed-topology distributions
- [ ] L274-276: rewrite for clarity
- [ ] L369: add intuitive explanation of why importance sampling corrects VI error
- [ ] L369-371: remove redundant sentence, move numbers to L375 with interpretation guide
- [ ] L428: add detail on BEAST tuning time, autotuning, convergence diagnostics
- [ ] L529-538: expand on structured posterior approximations
- [ ] L539-545: outline how topology inference could be addressed in future

## Phase 6: Software / installation fixes (R1 & R2 — blocking)
Fixes needed in the **treeflow** repo (not this paper repo):

- [ ] Fix install docs: remove `--user` from venv instructions
- [ ] Add `tf_keras` as explicit dependency
- [ ] Fix TF/TFP version incompatibility (`_get_attrs_items` error)
- [ ] Fix Docker build (`silence-tensorflow`/`support_developer` dependency)

## Phase 7: Response letter
- [ ] Draft point-by-point responses to all three reviewers + editor/AE
