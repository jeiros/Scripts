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
	# Fit to first frame, create average structure, save coordinates
	rms first
	average avg_${cluster}.${simtime}ns.nc
	createcrd CRD_${cluster}_${simtime}
	run
	# Load average as reference, fit saved coordinates to reference
	reference avg_${cluster}.${simtime}ns.nc [ref1]
	crdaction CRD_${cluster}_${simtime} rms @CA ref [ref1]

	# Calculate covariance matrix from rms-fit coords
	crdaction CRD_${cluster}_${simtime} matrix covar name matrixdat @CA out ${WORKDIR}/covmat-ca_${cluster}.${simtime}ns.dat

	# Diagonalize matrix, save eigenvectors as MyEvecsmhr7
	runanalysis diagmatrix matrixdat out ${WORKDIR}/evecs-ca_${cluster}.${simtime}ns.dat vecs 174 name MyEvecs${cluster}

	# Project fit and saved coordinates along eigenvectors
	crdaction CRD_${cluster}_${simtime} projection PROJECT modes MyEvecs${cluster} beg 1 end 2 @CA out ${WORKDIR}/myevecs_${cluster}.${simtime}ns.dat
	run
EOF

cd ${WORKDIR}