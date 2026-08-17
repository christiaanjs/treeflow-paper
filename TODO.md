# TreeFlow paper TODO

## Submission status

Revision expired at Syst Biol (USYB-2024-141) — Christiaan requesting
reopening from editor.

### Remaining before submission

- [ ] **Supplementary PDF**: compile `manuscript/tex/supplementary.tex` and commit
- [ ] Christiaan to resubmit once editor reopens

### Final submission files (editor's checklist)

`snakemake -s workflow/ms.smk ms_submission` assembles everything the editor
asked for into `manuscript/out/final-submission/` (and zips it to
`manuscript/out/final-submission.zip`): the main text as LaTeX plus PDF with the
figures removed and their legends after the reference list, the `.bib`, class
and bibliography style files, and one file per figure. `MANIFEST.txt` in the
bundle says what each file is.

- [ ] Run the target on a machine with R (Figure 5 is re-rendered as a vector
      PDF by `scripts/improved-benchmark-plot.R`, to meet the journal's 600 dpi
      minimum; the 300 dpi PNG does not)
- [ ] Figures 2-4 are still the committed 300 dpi PNGs (Figures 3 and 4 are
      PDFs with those PNGs embedded). If the editor asks for 600 dpi on those
      too, re-render the panels from the analysis outputs at a higher `dpi=` in
      the plotting scripts.

### Completed

- [x] treeflow PRs #75, #76, #79 merged (Keras 3 fixes, installation fixes)
- [x] H3N2 and Carnivores BEAST re-runs
- [x] Figures regenerated
- [x] Response letter drafted (`manuscript/response-letter.tex`)
- [x] All manuscript text changes verified
- [x] Software fixes on master (`--user`, `tf-keras`, `attrs.py`, Docker, numpy, READMEs)
- [x] Non-repo files sent to Christiaan (12 Mar)
