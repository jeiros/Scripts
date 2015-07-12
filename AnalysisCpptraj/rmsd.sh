#!/bin/bash


WORKDIR=$PWD

cluster=$1


simtime=000-250

cd ${WORKDIR}/${cluster}/run1/S1P/

cpptraj <<- EOF
	parm ./repstr.c0_phosS1P_nowat.prmtop
	trajin ./05_Prod*.nc
	rms rmsd @CA,C,O,N,H first out ${WORKDIR}/rmsd_${cluster}.${simtime}ns.dat mass
	run
EOF
cd $WORKDIR