#!/bin/bash




name=$1        # Can be CTnI_hmr CTnI_runs CTnT_hmr CTnT_runs with _run1 _run2 ...
                                # Example: CTnI_hmr-run3-S1P
                                # Example: WT-run3
prmtop=$2
trajs=$3


cpptraj.OMP <<- EOF
    parm ${prmtop}
    trajin ${trajs} 1 last 10
    rms first @CA,C,O,N
    average crdset average_structure @CA,C,O,N
    run
    # Step 2 - rms fit to average and calculate atomic flucts
    rms ref average_structure @CA,C,O,N
    atomicfluct out rmsf_${name}.dat @CA,C,O,N byres
    run
EOF
