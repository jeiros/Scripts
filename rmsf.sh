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
	rms ref [ref1] mass @CA,C,O,N,H
	atomicfluct out ${WORKDIR}/rmsf_${cluster}.${simtime}ns.dat @CA,C,O,N,H byres
	run
EOF
cd ${WORKDIR}

