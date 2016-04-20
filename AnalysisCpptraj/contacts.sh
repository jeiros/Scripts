#!/bin/bash



name=$1        # Can be CTnI_hmr CTnI_runs CTnT_hmr CTnT_runs with _run1 _run2 ...
                                # Example: CTnI_hmr-run3-S1P
                                # Example: WT-run3
prmtop=$2

if [ -z "$3" ]; then
    stride=1
else
    stride=$3
fi

trajs=$4

printf "\nThe selected arguments are:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
printf "%s\t\t\t\t\tSimulation time\n" ${simtime}
printf "%s\t\t\t\t\tName\n" ${name}
printf "%s\t\t\tTopology file\n" ${prmtop}
printf "%s\t Trajectories\n\n" ${trajs}


cpptraj.OMP <<- EOF
    parm ${prmtop}
    trajin ${trajs} 1 last ${stride}
    rms first
    average crdset average_structure
    run
    rms ref average_structure
    nativecontacts :1-89    :249-289 ref average_structure byresidue resout cmap_${name}_NcTnC-NcTnI.dat
    nativecontacts :1-89    :396-412 ref average_structure byresidue resout cmap_${name}_NcTnC-switch.dat

    nativecontacts :1-161   :386-395 ref average_structure byresidue resout cmap_${name}_cTnC-inhib.dat

    nativecontacts :249-289 :386-395 ref average_structure byresidue resout cmap_${name}_NcTnI-inhib.dat

    nativecontacts :232-248 :249-289 ref average_structure byresidue resout cmap_${name}_CcTnT-NcTnI.dat
    nativecontacts :232-248 :386-395 ref average_structure byresidue resout cmap_${name}_CcTnT-inhib.dat
 	run
EOF

for file in cmap*; do
    sort -k1 -n ${file} > ${file}_sorted
    mv ${file}_sorted ${file}
done