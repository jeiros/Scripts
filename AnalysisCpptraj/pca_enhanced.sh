#!/bin/bash
# Advanced analysis of PCA with cpptraj.OMP
name=$1        
prmtop=$2
stride=$3
mask=$4
name_mask=$5
trajs=$6

printf "Stride set to %s\n" ${stride}
cpptraj.OMP <<- EOF
	parm ${prmtop}
	trajin ${trajs} 1 last ${stride}
	# Fit to first frame, create average structure, save coordinates
	rms first
	average crdset average_structure
	createcrd loaded_trajs
	run
	#Fit the frames to the previously averaged structure, using the backbone
	crdaction loaded_trajs rms ref average_structure ${mask}

	# Calculate the coordinate covariance matrix
	# The input to PCA will be a coordinate covariance matrix. The entries to this matrix are the
	# covariance between the X, Y, and Z components of each atom, so the final matrix will have a size
	# of [3 * N_selected_atoms] X [3 * N_selected_atoms]. This means that in order to properly populate this
	# matrix we will need at least as many input frames to calculate the coordinate covariance matrix as we
	# have rows/columns (i.e. 3 * N_selected_atoms).
	crdaction loaded_trajs matrix covar name matrix_covar ${mask}

	#Diagonalize the cov matrix and get first 20 eigenvectors
	runanalysis diagmatrix matrix_covar out \
		./evecs-ca_${name}_${name_mask}.dat vecs 20 name \
		myEvecs nmwiz nmwizfile nmwiz.nmd nmwizvecs 20 nmwizmask ${mask}

	# Project fit and saved coordinates along eigenvectors
	crdaction loaded_trajs projection PROJECT modes myEvecs beg 1 end 20 \
		${mask} out ./myevecs_${name}_${name_mask}.dat

	# Histogram the 3 calculated projections with a Gaussian KDE
	
	kde PROJECT:1 out ./${name}_${name_mask}-kde.dat bins 400
	kde PROJECT:2 out ./${name}_${name_mask}-kde.dat bins 400
	kde PROJECT:3 out ./${name}_${name_mask}-kde.dat bins 400

	# Run and clear
	run
	clear all

	# Read the file with the eigenvectors
	readdata ./evecs-ca_${name}_${name_mask}.dat \
		name Evecs

	# Create a topology for the backbone only
	parm ./${prmtop}
	parmstrip !(${mask})
	parmwrite out ./${name_mask}.prmtop 

	# Create a .nc trajectory with the modes of motion of the 1st PC
	runanalysis modes name Evecs trajout ./${name}_${name_mask}_1PCA.nc \
		pcmin -100 pcmax 100 tmode 1 trajoutmask ${mask} trajoutfmt netcdf
	# Create a .nc trajectory with the modes of motion of the 2nd PC
	runanalysis modes name Evecs trajout ./${name}_${name_mask}_2PCA.nc \
		pcmin -100 pcmax 100 tmode 2 trajoutmask ${mask} trajoutfmt netcdf
	# Create a .nc trajectory with the modes of motion of the 3rd PC
	runanalysis modes name Evecs trajout ./${name}_${name_mask}_3PCA.nc \
		pcmin -100 pcmax 100 tmode 3 trajoutmask ${mask} trajoutfmt netcdf
	# Create a file with the eigenvalues
	runanalysis modes name Evecs eigenval out evalues.gnu
EOF
