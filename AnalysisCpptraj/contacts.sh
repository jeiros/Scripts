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
	rms first 
	average avg_${cluster}_${simtime}ns.nc
	run
	reference avg_${cluster}_${simtime}ns.nc [ref1]
	nativecontacts :1-161@CA :249-419@CA  distance 5.4 ref [ref1] byresidue resout ${WORKDIR}/TnCI_${cluster}.${simtime}ns.dat
	nativecontacts :1-161@CA :162-248@CA  distance 5.4 ref [ref1] byresidue resout ${WORKDIR}/TnCT_${cluster}.${simtime}ns.dat
 	nativecontacts :162-248@CA :249-419@C distance 5.4 ref [ref1] byresidue resout ${WORKDIR}/TnIT_${cluster}.${simtime}ns.dat
 	run
EOF
