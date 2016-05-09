#!/bin/bash

WORKDIR=$PWD
topology=${WORKDIR}/WT-ff14SB_clean.prmtop

traj=~/wt_data/*/05*nc


function create_dir {
	# Create a directory where all the clustering data is stored
	# Only 3 different masks can be given
	# :232-248 for C-terminus region of cTnT
	# :249-294 for N-terminus region of cTnI
	# :383-419 for C-terminus region of CTnI		
	mask=$1

	if [[ "$mask" == ":232-248" ]]; then
		echo "Mask set to CTnT"
		mkdir -p ./cluster_CTnT.runs_all
		cd ./cluster_CTnT.runs_all
	elif [[ "$mask" == ":249-294" ]]; then
		echo "Mask set to NTnI"
		mkdir -p ./cluster_NTnI.runs_all
		cd ./cluster_NTnI.runs_all
	elif [[ "$mask" == ":383-419" ]]; then
		echo "Mask set to CTnI"
		mkdir -p ./cluster_CTnI.runs_all
		cd ./cluster_CTnI.runs_all
	else
		printf "Mask was set to %s\n" $mask
		region=`echo $mask | tr -cd '[[:alnum:]]._-'`
		mkdir -p ./cluster_${region}.runs_all
		cd ./cluster_${region}.runs_all
	fi
}

function get_Kdists {
	# Quickly get the K-dists plots for k values ranging from 
	# 4 to 10. Use a sieve of 20 to load the data faster (doing
	# the plots does not require to be that accurate)
	for k in {4..10}
	do
		cpptraj.OMP <<- EOF
			parm $topology
			trajin ${traj} 1 last 10
			rms first
			cluster dbscan kdist $k rms $mask sieve 20
			run
		EOF
	done
}

function plotKdists {
	# Use gnuplot to graph the K dist plots obtained
	# with the get_Kdists function. Save the graph 
	# in Kplots.png		
	gnuplot <<- EOF
		set term png
		set output "Kplots.png"
		plot for [i=4:10] 'Kdist.'.i.'.dat' u 1:2 w l title "Kdist ".i
	EOF
}

function get_matrix {
	# Get the pairwise distance matrix so it doesn't have to be 
	# calculated every time in the get_clusters_dbscan function
	# The clustering metric is the all-atom RMSD				
	cpptraj.OMP <<- EOF
		parm $topology
		trajin $traj 1 last 10
		rms first
		cluster dbscan minpoints 4 epsilon 5.0 rms $mask savepairdist sieve 10 random
		run
	EOF
}



function get_clusters_dbscan {
	# Do the actual clustering on the previously save pairwise distance matrix
	# on RMSD space. k values are 4, 5 and 6, and eps goes from 2.5 to 4.0 in
	# 0.2 steps. A new directory is created for every combination
	cd $WORKDIR
	mask=$1

	# Check that mask is one of the 3 that have to be used
	# and cd into it's directory
	if [[ "$mask" == ":232-248" ]]; then
		region=CTnT
		cd ./cluster_CTnT.runs_all
	elif [[ "$mask" == ":249-294" ]]; then
		region=NTnI
		cd ./cluster_NTnI.runs_all
	elif [[ "$mask" == ":383-419" ]]; then
		region=CTnI
		cd ./cluster_CTnI.runs_all
	else
		region=`echo $mask | tr -cd '[[:alnum:]]._-'`
		cd ./cluster_${region}.runs_all
	fi
		
	for k in {4..6}
	do
		for eps in `seq 2.5 .2 4.0`
		do
			mkdir -p ./K.${k}.eps.${eps}
			cd ./K.${k}.eps.${eps}
			cpptraj.OMP <<- EOF
				parm $topology
				trajin $traj 1 last 10
				rms first
				cluster dbscan minpoints $k epsilon $eps rms mass $mask \
					sieve 10 random \
					loadpairdist pairdist ../CpptrajPairDist out cnumvtime.dat \
					summary avg.summary.dat info info.dat sil silhouette \
					summarysplit avg.summarysplit.dat \
					clusterout $region.nc clusterfmt netcdf \
					repout ${region}_repout.pdb repfmt pdb
				run
				EOF
				cd ../
			done
		done
}

function get_clusters_kmeans {

	cd $WORKDIR
	mask=$1
	# Check that mask is one of the 3 that have to be used
	if [[ "$mask" == ":232-248" ]]; then
		region=CTnT
		cd ./cluster_CTnT.runs_all
	elif [[ "$mask" == ":249-294" ]]; then
		region=NTnI
		cd ./cluster_NTnI.runs_all
	elif [[ "$mask" == ":383-419" ]]; then
		region=CTnI
		cd ./cluster_CTnI.runs_all
	else
		region=`echo $mask | tr -cd '[[:alnum:]]._-'`
		cd ./cluster_${region}.runs_all
	fi

	for cluster_number in `seq 1 10 100`
	do
		mkdir -p ./kmeans_clusters_${cluster_number}
		cd ./kmeans_clusters_${cluster_number}
		cpptraj.OMP <<- EOF
			parm $topology
			trajin $traj
			rms first
			cluster kmeans clusters ${cluster_number} randompoint rms mass $mask \
				# msieve 10 random \
				loadpairdist pairdist ../CpptrajPairDist out cnumvtime.dat \
				summary avg.summary.dat info info.dat sil silhouette \
				summarysplit avg.summarysplit.dat \
				clusterout $region.nc clusterfmt netcdf \
				repout ${region}_repout.pdb repfmt pdb
			run
		EOF
		cd ../
	done
}




# Call the functions
# The $1 argument is the mask that's going to be clustered
# has to be one of the options in the create_dir function

create_dir $1
# get_Kdists
# plotKdists
# get_matrix
# get_clusters_dbscan $1
# get_metrics_dbscan $1
get_clusters_kmeans $1



#cluster dbscan minpoints $k epsilon $eps rms mass $mask loadpairdist pairdist ../CpptrajPairDist out cnumvtime.dat summary avg.summary.dat info info.dat clusterout clusters clusterfmt netcdf sil silhouette summarysplit avg.summarysplit.dat sieve 10


##PLOT IN R

#ggplot(clustersdbipsf_sorted, aes (x = Clusters, y = DBI)) + geom_line() + ylab("DBI index") + xlab("Number of clusters") + scale_x_continuous(breaks = seq(2,60,by=2))
#ggsave(filename, path = )
