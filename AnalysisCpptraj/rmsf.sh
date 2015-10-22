#!/bin/bash



simtime=$1        # Can be 000-050, 000-0500, 100-150 ...
name=$2        # Can be CTnI_hmr CTnI_runs CTnT_hmr CTnT_runs with _run1 _run2 ...
                                # Example: CTnI_hmr-run3-S1P
                                # Example: WT-run3
prmtop=$3
trajs=$4


cpptraj <<- EOF
    parm ${prmtop}
    trajin ${trajs}
    rms first @CA,C,O,N,H
    average crdset average_structure @CA,C,O,N,H
    run
    # Step 2 - rms fit to average and calculate atomic flucts
    rms ref average_structure @CA,C,O,N,H
    atomicfluct out rmsf_${name}.${simtime}ns.dat @CA,C,O,N,H byres
    run
EOF
