#!/bin/bash
for count in {3..9}
do
cluster=$((2**$count))
msmb KCenters -i ../05_Prod_comphowcut_000-050ns_run3.nc -o assignments/CA_run3_KCenters_$cluster -t assignments/CA_run3_KCenters_$cluster --metric "rmsd" --top ../../comphowcut_c_20Anowat.pdb --n_clusters $cluster
done

