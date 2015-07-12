#!/bin/bash


WORKDIR=$PWD

cluster=$1
simtime=000-250


echo $cluster
echo $simtime

cd ${WORKDIR}/${cluster}/run1/S1P/

cpptraj <<- EOF
	parm ./repstr.c0_phosS1P_nowat.prmtop
	trajin ./05_Prod*.nc
	distance CAL420S69 :420 :69@OG out ${WORKDIR}/420-S69OG_${cluster}.${simtime}ns.dat
	distance CAL420T71 :420 :71@OG1 out ${WORKDIR}/420-T71OG1_${cluster}.${simtime}ns.dat
	distance CAL420D65OD1 :420 :65@OD1 out ${WORKDIR}/420-D65OD1_${cluster}.${simtime}ns.dat
	distance CAL420D65OD2 :420 :65@OD2 out ${WORKDIR}/420-D65OD2_${cluster}.${simtime}ns.dat
	distance CAL420D67OD1 :420 :67@OD1 out ${WORKDIR}/420-D67OD1_${cluster}.${simtime}ns.dat
	distance CAL420D67OD2 :420 :67@OD2 out ${WORKDIR}/420-D67OD2_${cluster}.${simtime}ns.dat
	distance CAL420D73OD1 :420 :73@OD1 out ${WORKDIR}/420-D73OD1_${cluster}.${simtime}ns.dat
	distance CAL420D73OD2 :420 :73@OD2 out ${WORKDIR}/420-D73OD2_${cluster}.${simtime}ns.dat
	distance CAL420E76OE1 :420 :76@OE1 out ${WORKDIR}/420-E76OE1_${cluster}.${simtime}ns.dat
	distance CAL420E76OE2 :420 :76@OE2 out ${WORKDIR}/420-E76OE2_${cluster}.${simtime}ns.dat
	distance K233-421 :233@NZ :421 out ${WORKDIR}/K233-421_${cluster}.${simtime}ns.dat
	distance K233-422 :233@NZ :422 out ${WORKDIR}/K233-422_${cluster}.${simtime}ns.dat
	distance K236-421 :236@NZ :421 out ${WORKDIR}/K236-421_${cluster}.${simtime}ns.dat
	distance K236-422 :236@NZ :422 out ${WORKDIR}/K236-422_${cluster}.${simtime}ns.dat
	distance K240D3 :240@NZ :3@CG out ${WORKDIR}/K240-D3_${cluster}.${simtime}ns.dat
	distance K240D2 :240@NZ :2@CG out ${WORKDIR}/K240-D2_${cluster}.${simtime}ns.dat
	distance K240CAL421 :240@NZ :421 out ${WORKDIR}/K240-421_${cluster}.${simtime}ns.dat
	distance K240CAL422 :240@NZ :422 out ${WORKDIR}/K240-422_${cluster}.${simtime}ns.dat
	distance K242D3 :242@NZ :3@CG out ${WORKDIR}/K242-D3_${cluster}.${simtime}ns.dat
	distance K242D2 :242@NZ :2@CG out ${WORKDIR}/K242-D2_${cluster}.${simtime}ns.dat
	distance K242CAL421 :242@NZ :421 out ${WORKDIR}/K242-421_${cluster}.${simtime}ns.dat
	distance K242CAL422 :242@NZ :422 out ${WORKDIR}/K242-422_${cluster}.${simtime}ns.dat
	distance 30-38_271-272 :30-38 :271-272 out ${WORKDIR}/30-38_271-272_${cluster}.${simtime}ns.dat geom
	distance 271-272_CAL420 :271-272 :420 out ${WORKDIR}/271-272_420_${cluster}.${simtime}ns.dat geom
	distance 284-286_CAL420 :284-286 :420 out ${WORKDIR}/284-286_420_${cluster}.${simtime}ns.dat geom
	distance 420-271O2P :420 :271@O2P out ${WORKDIR}/420-271@O2P_${cluster}.${simtime}ns.dat
	distance 420-271O3P :420 :271@O3P out ${WORKDIR}/420-271@O3P_${cluster}.${simtime}ns.dat
	distance 420-272O2P :420 :272@O2P out ${WORKDIR}/420-272@O2P_${cluster}.${simtime}ns.dat
	distance 420-272O3P :420 :272@O3P out ${WORKDIR}/420-272@O3P_${cluster}.${simtime}ns.dat
	run
EOF