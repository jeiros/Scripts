#!/bin/bash



name=$1        # Can be CTnI_hmr CTnI_runs CTnT_hmr CTnT_runs with _run1 _run2 ...
                                # Example: CTnI_hmr-run3-S1P
                                # Example: WT-run3
prmtop=$2
trajs=$3


cpptraj.OMP <<- EOF
	parm ${prmtop}
	trajin ${trajs}
	rms rmsd @CA,C,O,N first out ./rmsd_${name}ns.dat mass
	run
EOF
