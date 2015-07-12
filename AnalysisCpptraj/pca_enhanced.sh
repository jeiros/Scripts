#!/bin/bash


WORKDIR=$PWD

cluster=$1
run=$2
simtime=$3


DESTINATION=${WORKDIR}/${cluster}/${run}/S1P/


if [ -d "$DESTINATION" ]; then 
	printf "The directory exists\n\n"
else
	printf "The destination directory does not exist\n"
	printf "\nPlease check the arguments
--------------------------------

--------------------------------\n" 
	printf "Exiting now...\n\n"
	exit 1
fi



cd ${WORKDIR}/${cluster}/${run}/S1P/


cpptraj <<- EOF
	parm ./repstr.c0_phosS1P_nowat.prmtop
	trajin ./05_Prod*.nc
	# Fit to first frame, create average structure, save coordinates
	rms first
	average crdset average_structure
	createcrd loaded_trajs
	run
	#Fit the frames to the previously averaged structure, using the backbone
	crdaction loaded_trajs rms ref average_structure @CA,C,O,N,H
	#Calculate the coordinate covariance matrix
	crdaction loaded_trajs matrix covar name matrix_covar @CA,C,O,N,H
	#Diagonalize the cov matrix and get first 3 eigenvectors
	runanalysis diagmatrix matrix_covar out ${WORKDIR}/evecs-ca_${cluster}.${simtime}ns.dat \
		vecs 3 name myEvecs @CA,C,O,N,H
	# Project fit and saved coordinates along eigenvectors
	crdaction loaded_trajs projection PROJECT modes Evecs \
		beg 1 end 3 @CA,C,O,N,H
	# Make normalized histogram of the 3 calculated projections
	hist PROJECT:1 bins 100 out ${cluster}_${simtime}-hist.dat norm name PROJECT:1
	hist PROJECT:2 bins 100 out ${cluster}_${simtime}-hist.dat norm name PROJECT:2
	hist PROJECT:3 bins 100 out ${cluster}_${simtime}-hist.dat norm name PROJECT:3
	# Run and clear
	run
	clear all
	# Read the file with the eigenvectors
	readdata ${WORKDIR}/evecs-ca_${cluster}.${simtime}ns.dat

	# Create a topology for the backbone only
	parm ./repstr.c0_phosS1P_nowat.prmtop
	parmstrip !(@CA,C,O,N,H)
	parmwrite out ./repstr.c0_phosS1P_nowat_striped.prmtop 

	#Create a .nc trajectory with the modes of motion of the 1st PC
	runanalysis modes name Evecs trajout ./${cluster}_${simtime}_1PCA.nc \
		pcmin -100 pcmax 100 tmode 1 trajoutmask @CA,C,O,N,H trajoutfmt netcdf

EOF



