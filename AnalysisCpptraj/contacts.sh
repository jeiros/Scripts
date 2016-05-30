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


/usr/local/amber15/bin/cpptraj.OMP <<- EOF
    parm ${prmtop}
    trajin ${trajs} 1 last ${stride}
    rms first
    average crdset average_structure
    run
    rms ref average_structure
    nativecontacts :1-89    :249-289 ref average_structure byresidue resout cmap_NcTnC-NcTnI_${name}.dat
    nativecontacts :1-89    :396-412 ref average_structure byresidue resout cmap_NcTnC-switch_${name}.dat

    nativecontacts :1-161   :386-395 ref average_structure byresidue resout cmap_cTnC-inhib_${name}.dat

    nativecontacts :249-289 :386-395 ref average_structure byresidue resout cmap_NcTnI-inhib_${name}.dat

    nativecontacts :232-248 :249-289 ref average_structure byresidue resout cmap_CcTnT-NcTnI_${name}.dat
    nativecontacts :232-248 :386-395 ref average_structure byresidue resout cmap_CcTnT-inhib_${name}.dat
 	run
EOF

# cpptraj.OMP <<- EOF
#     parm ${prmtop} [traj_top]
    
#     parm /Users/je714/Troponin/IAN_Troponin/completehowarthcut/salted/completehowcut_correct.prmtop [pdb_top]  # Topology matching the reference
#     reference /Users/je714/Troponin/IAN_Troponin/completehowarthcut/salted/completehowcut_correct.inpcrd parm [pdb_top] [pdb_ref] # The reference coordinates
    
#     trajin ${trajs} 1 last ${stride} parm [traj_top]

#     nativecontacts :1-89    :249-289 ref [pdb_ref] writecontacts atom_contacts.dat\
#         map mapout resmap.gnu contactpdb NcTnC-NcTnI.pdb \
#         byresidue resout cmap_NcTnC-NcTnI_${name}.dat


#     nativecontacts :1-89    :396-412 ref [pdb_ref] byresidue resout cmap_NcTnC-switch_${name}.dat

#     nativecontacts :1-161   :386-395 ref [pdb_ref] byresidue resout cmap_cTnC-inhib_${name}.dat

#     nativecontacts :249-289 :386-395 ref [pdb_ref] byresidue resout cmap_NcTnI-inhib_${name}.dat

#     nativecontacts :232-248 :249-289 ref [pdb_ref] byresidue resout cmap_CcTnT-NcTnI_${name}.dat
#     nativecontacts :232-248 :386-395 ref [pdb_ref] byresidue resout cmap_CcTnT-inhib_${name}.dat
#     run
# EOF

for file in cmap*; do
    sort -k1 -n ${file} > ${file}_sorted
    mv ${file}_sorted ${file}
done