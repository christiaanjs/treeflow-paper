# TreeFlow paper TODO

## H3N2 BEAST re-run (clean timing)

Previous run had correct speed (13.7 min/M with BEAGLE CPU arm64) but laptop
went to sleep after ~17M iterations, inflating the cumulative average to 26 min/M.

```bash
# Clean outputs and run with caffeinate to prevent sleep (expects venv + ~/bin/beast on PATH)
cd ~/Git/treeflow-paper
source .venv/bin/activate
export PATH="$HOME/bin:$PATH"
caffeinate -s nohup snakemake -s workflow/data.smk out/h3n2/beast.log out/h3n2/beast.trees --rerun-triggers mtime --cores 1 > out/h3n2/snakemake.log 2>&1 &
```

Expected: ~13.7 min/M, total ~6.8 hours for 30M iterations.

Check progress with: `tail -5 out/h3n2/beast-log.txt`
Check timing with: `python scripts/check-beast-timing.py out/h3n2/beast-log.txt`

## Carnivores BEAST re-run

Same approach, after H3N2 is done:

```bash
rm -f out/carnivores/beast.log out/carnivores/beast.trees out/carnivores/beast-log.txt out/carnivores/beast-benchmark.txt
caffeinate -s nohup snakemake -s workflow/data.smk out/carnivores/beast.log out/carnivores/beast.trees --rerun-triggers mtime --cores 1 > out/carnivores/snakemake.log 2>&1 &
```

## After both runs complete

- Generate timing CSVs: `snakemake -s workflow/data.smk out/h3n2/timing-data.csv out/carnivores/timing-data.csv --rerun-triggers mtime --cores 1`
- Compile manuscript: `snakemake -s workflow/ms.smk --cores 1`
- ~~Write response letter~~ Draft: `manuscript/response-letter.tex`

## Response letter verification (11 Mar 2026)

### Blockers before submission

- [ ] **treeflow PR #76** (draft): not merged — need to run VI first to verify results/speed match previous submission
- [ ] **Keras 3 fix** (`fix/keras3-robust-optimizer` branch): not on master yet — merge into PR #76 or separately
- [ ] **H3N2 BEAST re-run**: running on laptop (ETA ~4:45pm 11 Mar)
- [ ] **Carnivores BEAST re-run**: not started — run after H3N2 or on second machine
- [ ] **Regenerate figures** (Figs 2, 6a, 7): blocked on BEAST re-runs
- [ ] **Run treeflow test suite** on Apple Silicon (TF 2.20+) to get current pass/fail counts

### Claims to verify in manuscript

- [ ] RAxML in opening paragraph (L13-14): was it replaced with PhyloBayes there, or only added elsewhere? RAxML still at L170 (different context — ML topology — probably fine)
- [ ] CLI docs: were example invocations actually added, or just auto-generated sphinx-click?
- [ ] "it's" → "its" fix: line numbers shifted, confirm it's in the current text

### Claims verified ✓

All other manuscript text changes (typos, references, structural reorganisation, Discussion reorder, supplementary appendix, Burda citation, bijector explanation, structured approximations, subsplit BNs) confirmed present in `treeflow.j2.tex` and `main.bib`.

Software fixes confirmed on master: `--user` removed, `tf-keras` added, `attrs.py` fixed, Docker fixed, numpy pinned, READMEs written.
