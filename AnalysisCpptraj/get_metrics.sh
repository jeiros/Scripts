#!/bin/bash

directories=$1


cat > MetricsFile.dat <<-EOF
DBI and PSF metrics
Dirs are ${directories}

EOF

for dir in ${directories}; do
    echo $dir
    nlines=`wc -l < ./${dir}/avg.summary.dat`
    # rc=$?
    # if [[ $rc != 0 ]]; then continue; fi
    nclusters=$((nlines-1))
    echo $nclusters
    if [[ $nclusters -eq -1 ]]; then
        $nclusters="NaN"
        echo "derp"
    fi
    printf "Clusters:%s\n" $nclusters >> MetricsFile.dat
    DBI=`grep -E 'DBI' ./${dir}/info.dat | grep -o '[0-9]*\.[0-9]*'`
    echo $DBI
    pSF=`grep -E 'pSF' ./${dir}/info.dat | grep -o '[0-9]*\.[0-9]*'`
    echo $pSF
    if [[ "$pSF" == "" ]]; then
        pSF="NaN"
    fi
    if [[ "$DBI" == "" ]]; then
        DBI="NaN"
    fi
    printf "DBI index:\t%s\n" $DBI >> MetricsFile.dat
    printf "pSF index:\t%s\n\n" $pSF >> MetricsFile.dat
done

grep -E 'DBI' MetricsFile.dat | grep -o '[0-9]*\.[0-9]*' > dbi.dat
grep -E 'pSF' MetricsFile.dat | cut -b 12-22 > psf.dat
grep -E 'Clusters' MetricsFile.dat | cut -d ":" -f 2 > clusters.dat

paste clusters.dat dbi.dat psf.dat > clusters_dbipsf.dat
sort -s -n -k 1,1 clusters_dbipsf.dat > clusters_dbipsf_sorted.dat

R --vanilla<<-EOF
source("~/Scripts/Rplots/plot_DBI_pSF.R")
clusters <- read.table("./clusters_dbipsf_sorted.dat")
data <- reduce_matrix(clusters)
plotdbipsf(data)
q()
EOF