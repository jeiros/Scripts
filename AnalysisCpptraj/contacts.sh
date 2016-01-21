#!/bin/bash


simtime=$1        # Can be 000-050, 000-0500, 100-150 ...
name=$2        # Can be CTnI_hmr CTnI_runs CTnT_hmr CTnT_runs with _run1 _run2 ...
                                # Example: CTnI_hmr-run3-S1P
                                # Example: WT-run3
prmtop=$3

trajs=$4


printf "\nThe selected arguments are:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
printf "%s\t\t\t\t\tSimulation time\n" ${simtime}
printf "%s\t\t\t\t\tName\n" ${name}
printf "%s\t\t\tTopology file\n" ${prmtop}
printf "%s\t Trajectories\n\n" ${trajs}


cpptraj.OMP <<- EOF
    parm ${prmtop}
    trajin ${trajs}
    rms first
    average crdset average_structure
    run
    rms ref average_structure
    nativecontacts :1-89    :396-412 ref average_structure byresidue resout cmap_${name}.${simtime}ns_NcTnC-switch.dat
    nativecontacts :1-161   :386-395 ref average_structure byresidue resout cmap_${name}.${simtime}ns_cTnC-inhib.dat
    nativecontacts :1-89    :249-289 ref average_structure byresidue resout cmap_${name}.${simtime}ns_NcTnC-NcTnI.dat

    nativecontacts :249-289 :386-395 ref average_structure byresidue resout cmap_${name}.${simtime}ns_NcTnI-inhib.dat

    nativecontacts :232-248 :249-289 ref average_structure byresidue resout cmap_${name}.${simtime}ns_CcTnT-NcTnI.dat
    nativecontacts :232-248 :386-395 ref average_structure byresidue resout cmap_${name}.${simtime}ns_CcTnT-inhib.dat
 	run
EOF
