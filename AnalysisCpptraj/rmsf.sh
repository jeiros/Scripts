#!/bin/bash


WORKDIR=$PWD

simtime=$1        # Can be 000-050, 000-0500, 100-150 ...
name=$2        # Can be CTnI_hmr CTnI_runs CTnT_hmr CTnT_runs with _run1 _run2 ...
                                # Example: CTnI_hmr-run3-S1P
                                # Example: WT-run3
prmtop=$3
trajs=$4





# OPTION 1
# cpptraj <<- EOF
#     parm ${prmtop}
#     trajin ${trajs}
#     loadtraj name loaded_trajs
#     crdaction loaded_trajs average crdset average_structure @CA,C,O,N,H
#     crdaction loaded_trajs rms ref average_structure @CA,C,O,N,H
#     # crdaction loaded_trajs atomicfluct out rmsf_OPT1_${name}.${simtime}ns.dat @CA,C,O,N,H byres
#     crdaction loaded_trajs atomicfluct out rmsf_OPT1.dat @CA,C,O,N,H byres
#     run
# EOF

# OPTION 1

# cpptraj <<- EOF
#     parm ${prmtop}
#     trajin ${trajs}
    # rms first
    # average crdset average_structure
    # createcrd loaded_trajs
    # run
    # crdaction loaded_trajs rms ref average_structure @CA,C,O,N,H
    # atomicfluct out rmsf_OPT.dat @CA,C,O,N,H byres
    # run
# EOF

# OPTION 2
cpptraj <<- EOF
    parm ${prmtop}
    loadcrd ${trajs}
    crdaction ${trajs} average avg.pdb @CA,C,O,N,H
    parm avg.pdb
    reference avg.pdb parm avg.pdb
    crdaction ${trajs} rms reference @CA,C,O,N,H
    crdaction ${trajs} atomicfluct out rmsf_OPT2.dat @CA,C,O,N,H byres
    run
EOF


# Dan Roe


cpptraj <<- EOF
    # Step 1 - create average structure, rms-fit to first frame
    parm ${prmtop}
    trajin ${trajs}
    rms first @CA,C,O,N,H
    average crdset average_structure @CA,C,O,N,H
    run
    # Step 2 - rms fit to average and calculate atomic flucts
    rms ref average_structure @CA,C,O,N,H
    atomicfluct out rmsf_danroe.dat @CA,C,O,N,H byres
    run
EOF
