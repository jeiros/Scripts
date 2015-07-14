#!/bin/bash

# Advanced analysis of PCA with cpptraj
# Must run on the directory where the trajectories
# and the prmtop are. 





simtime=$1		  # Can be 000-050, 000-0500, 100-150 ...
name=$2        # Can be CTnI_hmr CTnI_runs CTnT_hmr CTnT_runs with _run1 _run2 ...
								# Example: CTnI_hmr-run3-S1P
								# Example: WT-run3

prmtop=$3
trajs=$4



cpptraj <<- EOF
	parm ./${prmtop}
	trajin ./${trajs}
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
	runanalysis diagmatrix matrix_covar out \
		./evecs-ca_${name}.${simtime}ns.dat vecs 200 name \
		myEvecs

	# Project fit and saved coordinates along eigenvectors
	crdaction loaded_trajs projection PROJECT modes myEvecs beg 1 end 3 \
		@CA,C,O,N,H out ./myevecs_${name}.${simtime}ns.dat

	# Make normalized histogram of the 3 calculated projections
	hist PROJECT:1 bins 100 out \
		./${name}_${simtime}-hist.dat norm name PROJECT-1

	hist PROJECT:2 bins 100 out \
		./${name}_${simtime}-hist.dat norm name PROJECT-2

	hist PROJECT:3 bins 100 out \
		./${name}_${simtime}-hist.dat norm name PROJECT-3

	# Run and clear
	run
	clear all

	# Read the file with the eigenvectors
	readdata ./evecs-ca_${name}.${simtime}ns.dat \
		name Evecs

	# Create a topology for the backbone only
	parm ./${prmtop}
	parmstrip !(@CA,C,O,N,H)
	parmwrite out ./repstr.c0_phosS1P_nowat_backbone.prmtop 

	#Create a .nc trajectory with the modes of motion of the 1st PC
	runanalysis modes name Evecs trajout ./${name}_${simtime}_1PCA.nc \
		pcmin -100 pcmax 100 tmode 1 trajoutmask @CA,C,O,N,H trajoutfmt netcdf

EOF