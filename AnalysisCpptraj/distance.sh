#!/bin/bash

simtime=$1		  # Can be 000-050, 000-0500, 100-150 ...
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

	# Coordination sphere catalytic calcium
	distance CAL420D65OD1 :420 :65@OD1 out ./420-D65OD1_${name}.${simtime}ns.dat
	distance CAL420D65OD2 :420 :65@OD2 out ./420-D65OD2_${name}.${simtime}ns.dat

	distance CAL420E66OE1 :420 :66@OE1 out ./420-E66OE1_${name}.${simtime}ns.dat
	distance CAL420E66OE2 :420 :66@OE2 out ./420-E66OE2_${name}.${simtime}ns.dat

	distance CAL420D67OD1 :420 :67@OD1 out ./420-D67OD1_${name}.${simtime}ns.dat
	distance CAL420D67OD2 :420 :67@OD2 out ./420-D67OD2_${name}.${simtime}ns.dat

	distance CAL420S69 	  :420 :69@OG  out ./420-S69OG_${name}.${simtime}ns.dat

	distance CAL420T71    :420 :71@OG1 out ./420-T71OG1_${name}.${simtime}ns.dat

	distance CAL420D73OD1 :420 :73@OD1 out ./420-D73OD1_${name}.${simtime}ns.dat
	distance CAL420D73OD2 :420 :73@OD2 out ./420-D73OD2_${name}.${simtime}ns.dat

	distance CAL420D75OD1 :420 :75@OD1 out ./420-D75OD1_${name}.${simtime}ns.dat
	distance CAL420D75OD2 :420 :75@OD2 out ./420-D75OD2_${name}.${simtime}ns.dat
	
	distance CAL420E76OE1 :420 :76@OE1 out ./420-E76OE1_${name}.${simtime}ns.dat
	distance CAL420E76OE2 :420 :76@OE2 out ./420-E76OE2_${name}.${simtime}ns.dat


	# distance K233-421 :233@NZ :421 out ./K233-421_${name}.${simtime}ns.dat
	# distance K233-422 :233@NZ :422 out ./K233-422_${name}.${simtime}ns.dat
	# distance K236-421 :236@NZ :421 out ./K236-421_${name}.${simtime}ns.dat
	# distance K236-422 :236@NZ :422 out ./K236-422_${name}.${simtime}ns.dat
	# distance K240D3 :240@NZ :3@CG out ./K240-D3_${name}.${simtime}ns.dat
	# distance K240D2 :240@NZ :2@CG out ./K240-D2_${name}.${simtime}ns.dat
	# distance K240CAL421 :240@NZ :421 out ./K240-421_${name}.${simtime}ns.dat
	# distance K240CAL422 :240@NZ :422 out ./K240-422_${name}.${simtime}ns.dat
	# distance K242D3 :242@NZ :3@CG out ./K242-D3_${name}.${simtime}ns.dat
	# distance K242D2 :242@NZ :2@CG out ./K242-D2_${name}.${simtime}ns.dat
	# distance K242CAL421 :242@NZ :421 out ./K242-421_${name}.${simtime}ns.dat
	# distance K242CAL422 :242@NZ :422 out ./K242-422_${name}.${simtime}ns.dat
	# distance 30-38_271-272 :30-38 :271-272 out ./30-38_271-272_${name}.${simtime}ns.dat geom
	# distance 271-272_CAL420 :271-272 :420 out ./271-272_420_${name}.${simtime}ns.dat geom
	# distance 284-286_CAL420 :284-286 :420 out ./284-286_420_${name}.${simtime}ns.dat geom
	# distance 420-271O2P :420 :271@O2P out ./420-271@O2P_${name}.${simtime}ns.dat
	# distance 420-271O3P :420 :271@O3P out ./420-271@O3P_${name}.${simtime}ns.dat
	# distance 420-272O2P :420 :272@O2P out ./420-272@O2P_${name}.${simtime}ns.dat
	# distance 420-272O3P :420 :272@O3P out ./420-272@O3P_${name}.${simtime}ns.dat
	run
EOF
