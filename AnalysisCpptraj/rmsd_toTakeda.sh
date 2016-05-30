#!/bin/bash




name=$1        # Can be CTnI_hmr CTnI_runs CTnT_hmr CTnT_runs with _run1 _run2 ...
                                # Example: CTnI_hmr-run3-S1P
                                # Example: WT-run3
prmtop=$2
trajs=$3


cpptraj.OMP <<- EOF
    parm ${prmtop} [traj_top]
    parm /Users/je714/Troponin/IAN_Troponin/completehowarthcut/salted/ff99SB/1j1d.prmtop [pdb_top]  # Topology matching the reference
    reference /Users/je714/Troponin/IAN_Troponin/completehowarthcut/salted/ff99SB/1j1d.inpcrd parm [pdb_top] [pdb_ref] # The reference coordinates
    trajin ${trajs} parm [traj_top]
    rms :1-89,92-161,162-231,283-384,393-408@CA,C,O,N @CA,C,O,N out rmsd_${name}_toTakeda.dat mass ref [pdb_ref]
    run
EOF
