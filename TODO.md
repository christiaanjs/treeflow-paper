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
- Write response letter
