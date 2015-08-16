#!/bin/bash

WORKDIR=$PWD
topology=${WORKDIR}/striped.prmtop

traj=${WORKDIR}/runs_all.nc



function get_Kdists {
for k in {4..10}
do
	cpptraj <<- EOF
		parm $topology
		trajin $traj 1 last 10
		rms first
		cluster dbscan kdist $k rms $mask sieve 20
		run
	EOF
done
}

function plotKdists {
gnuplot <<- EOF
	set term png
	set output "Kplots.png"
	plot for [i=4:10] 'Kdist.'.i.'.dat' u 1:2 w l title "Kdist ".i
EOF
}

function get_matrix {
cpptraj <<- EOF
	parm $topology
	trajin $traj 1 last 10
	rms first
	cluster dbscan minpoints 4 epsilon 5.0 rms $mask savepairdist sieve 10
	run
EOF
}



function get_clusters {
for k in {4..6}
do
for eps in `seq 3.0 .2 6.0`
do
mkdir -p ./K.${k}.eps.${eps}
cd ./K.${k}.eps.${eps}
cpptraj <<- EOF
	parm $topology
	trajin $traj 1 last 10
	rms first
	cluster dbscan minpoints $k epsilon $eps rms mass $mask loadpairdist pairdist ../CpptrajPairDist out cnumvtime.dat summary avg.summary.dat info info.dat sil silhouette summarysplit avg.summarysplit.dat sieve 10
	run
EOF
cd ../
done
done
}

function create_dir {

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
		echo "Please enter a valid mask"
		exit
	fi
}

function get_metrics {
cd $WORKDIR
mask=$1

if [[ "$mask" == ":232-248" ]]; then
	region=CTnT
	cd ./cluster_CTnT.runs_all
	WORKDIR=$PWD
	echo "Working directory set to $WORKDIR"
elif [[ "$mask" == ":249-294" ]]; then
	region=NTnI
	cd ./cluster_NTnI.runs_all
	WORKDIR=$PWD
	echo "Working directory set to $WORKDIR"
elif [[ "$mask" == ":383-419" ]]; then
	region=CTnI
	cd ./cluster_CTnI.runs_all
	WORKDIR=$PWD
	echo "Working directory set to $WORKDIR"
else
	echo "Please enter a valid mask"
	exit
fi

cat > MetricsFile.dat <<EOF
DBI and PSF metrics for ${region}
Mask is set to ${mask}

EOF

for k in {4..6}
do
	for eps in `seq 3.0 .2 6.0`
	do
		printf "Minpoints:\t%s\n" $k >> MetricsFile.dat
		printf "Epsilon:\t%s\n" $eps >> MetricsFile.dat
		nlines=`wc -l < ./K.${k}.eps.${eps}/avg.summary.dat`
		nclusters=$((nlines-1))
		printf "Clusters:\t%s\n" $nclusters >> MetricsFile.dat

		DBI=`grep -E 'DBI' ./K.${k}.eps.${eps}/info.dat | grep -o '[0-9]*\.[0-9]*'`
		pSF=`grep -E 'pSF' ./K.${k}.eps.${eps}/info.dat | grep -o '[0-9]*\.[0-9]*'`

		if [[ "$pSF" == "" ]]; then
			pSF="NaN"
		else
			pSF=`grep -E 'pSF' ./K.${k}.eps.${eps}/info.dat | grep -o '[0-9]*\.[0-9]*'`
		fi



		printf "DBI index:\t%s\n" $DBI >> MetricsFile.dat
		printf "pSF index:\t%s\n\n" $pSF >> MetricsFile.dat
	done
done

grep -E 'DBI' MetricsFile.dat | grep -o '[0-9]*\.[0-9]*' > dbi.dat
grep -E 'pSF' MetricsFile.dat | cut -b 12-22 > psf.dat
grep -E 'Clusters' MetricsFile.dat | grep -o '[0-9]*' > clusters.dat

paste clusters.dat dbi.dat psf.dat > clusters_dbipsf.dat
sort -s -n -k 1,1 clusters_dbipsf.dat > clusters_dbipsf_sorted.dat


R --vanilla<<EOF
source("../plot_DBI_pSF.R")
clusters <- read.delim("./clusters_dbipsf_sorted.dat", header = FALSE, na.strings = "NaN")
data <- reduce_matrix(clusters)
plotdbipsf(data)
q()
EOF
}





# Call the functions
create_dir $1
get_Kdists
plotKdists
get_matrix
get_clusters
get_metrics $1


#cluster dbscan minpoints $k epsilon $eps rms mass $mask loadpairdist pairdist ../CpptrajPairDist out cnumvtime.dat summary avg.summary.dat info info.dat clusterout clusters clusterfmt netcdf sil silhouette summarysplit avg.summarysplit.dat sieve 10


##PLOT IN R

#ggplot(clustersdbipsf_sorted, aes (x = Clusters, y = DBI)) + geom_line() + ylab("DBI index") + xlab("Number of clusters") + scale_x_continuous(breaks = seq(2,60,by=2))
#ggsave(filename, path = )