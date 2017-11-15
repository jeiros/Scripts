#!/bin/bash

name=$1        # Can be CTnI_hmr CTnI_runs CTnT_hmr CTnT_runs with _run1 _run2 ...
								# Example: CTnI_hmr-run3-S1P
								# Example: WT-run3
prmtop=$2
trajs=$3

printf "\nThe selected arguments are:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
printf "%s\t\t\t\t\tName\n" ${name}
printf "%s\t\t\tTopology file\n" ${prmtop}
printf "%s\t Trajectories\n\n" ${trajs}


mkdir -p dist
cpptraj.OMP <<- EOF
	parm ${prmtop}
	trajin ${trajs}

	# Coordination sphere catalytic calcium
	distance CAL420D65OD1 :420 :65@OD1 out ./dist/420-D65OD1_${name}.dat
	distance CAL420D65OD2 :420 :65@OD2 out ./dist/420-D65OD2_${name}.dat
	distance CAL420D65O   :420 :65@O   out ./dist/420-D65O_${name}.dat

	distance CAL420E66OE1 :420 :66@OE1 out ./dist/420-E66OE1_${name}.dat
	distance CAL420E66OE2 :420 :66@OE2 out ./dist/420-E66OE2_${name}.dat
	distance CAL420E66O   :420 :66@O   out ./dist/420-E66O_${name}.dat

	distance CAL420D67OD1 :420 :67@OD1 out ./dist/420-D67OD1_${name}.dat
	distance CAL420D67OD2 :420 :67@OD2 out ./dist/420-D67OD2_${name}.dat
	distance CAL420D67O   :420 :67@O   out ./dist/420-D67O_${name}.dat

	distance CAL420S69OG  :420 :69@OG  out ./dist/420-S69OG_${name}.dat
	distance CAL420S69O   :420 :69@O  out ./dist/420-S69O_${name}.dat

	distance CAL420T71    :420 :71@OG1 out ./dist/420-T71OG1_${name}.dat
	distance CAL4240T71O  :420 :71@O   out ./dist/420-T71O_${name}.dat

	distance CAL420D73OD1 :420 :73@OD1 out ./dist/420-D73OD1_${name}.dat
	distance CAL420D73OD2 :420 :73@OD2 out ./dist/420-D73OD2_${name}.dat
	distance CAL420D73O   :420 :73@O   out ./dist/420-D73O_${name}.dat

	distance CAL420D75OD1 :420 :75@OD1 out ./dist/420-D75OD1_${name}.dat
	distance CAL420D75OD2 :420 :75@OD2 out ./dist/420-D75OD2_${name}.dat
	distance CAL420D75O   :420 :75@O   out ./dist/420-D75O_${name}.dat
	
	distance CAL420E76OE1 :420 :76@OE1 out ./dist/420-E76OE1_${name}.dat
	distance CAL420E76OE2 :420 :76@OE2 out ./dist/420-E76OE2_${name}.dat
	distance CAL420E76O   :420 :76@O   out ./dist/420-E76O_${name}.dat

	distance CAL420E280OE1 :420 :280@OE1 out ./dist/420-E280OE1_${name}.dat
	distance CAL420E280OE2 :420 :280@OE2 out ./dist/420-E280OE2_${name}.dat
	distance CAL420E280O   :420 :280@O   out ./dist/420-E280O_${name}.dat


	# distance K233-421 :233@NZ :421 out ./dist/K233-421_${name}.dat
	# distance K233-422 :233@NZ :422 out ./dist/K233-422_${name}.dat
	# distance K236-421 :236@NZ :421 out ./dist/K236-421_${name}.dat
	# distance K236-422 :236@NZ :422 out ./dist/K236-422_${name}.dat
	# distance K240D3 :240@NZ :3@CG out ./dist/K240-D3_${name}.dat
	# distance K240D2 :240@NZ :2@CG out ./dist/K240-D2_${name}.dat
	# distance K240CAL421 :240@NZ :421 out ./dist/K240-421_${name}.dat
	# distance K240CAL422 :240@NZ :422 out ./dist/K240-422_${name}.dat
	# distance K242D3 :242@NZ :3@CG out ./dist/K242-D3_${name}.dat
	# distance K242D2 :242@NZ :2@CG out ./dist/K242-D2_${name}.dat
	# distance K242CAL421 :242@NZ :421 out ./dist/K242-421_${name}.dat
	# distance K242CAL422 :242@NZ :422 out ./dist/K242-422_${name}.dat
	# distance 30-38_271-272 :30-38 :271-272 out ./dist/30-38_271-272_${name}.dat geom
	# distance 271-272_CAL420 :271-272 :420 out ./dist/271-272_420_${name}.dat geom
	# distance 284-286_CAL420 :284-286 :420 out ./dist/284-286_420_${name}.dat geom
	# distance 420-271O2P :420 :271@O2P out ./dist/420-271@O2P_${name}.dat
	# distance 420-271O3P :420 :271@O3P out ./dist/420-271@O3P_${name}.dat
	# distance 420-272O2P :420 :272@O2P out ./dist/420-272@O2P_${name}.dat
	# distance 420-272O3P :420 :272@O3P out ./dist/420-272@O3P_${name}.dat
	run
EOF
