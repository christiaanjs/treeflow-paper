module sim:
    snakefile: "workflow/sim.smk"

use rule * from sim as sim_*
