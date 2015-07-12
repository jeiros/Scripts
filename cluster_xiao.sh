#!/bin/bash

#./script.sh TotalFrames NoiseFrames Mask(Default is @CA)


 echo -e "\n~~~~~~~Cpptraj Cluster Script written by X. Hu April 2015. Please make sure it fit your requirement.~~~~~~~"
 echo "**Please ensure home directory name consistency and regularity as a requirement for using this script !!" 
################################################Debug Setting####################################################

trap 'echo ~~~~~~~Exit normally ...~~~~~~~' EXIT
trap 'echo ~~~~~~~Exit with Error !! Please verify..~~~~~~~~' ERR

##################################################Set Up#########################################################

offset=10
sieve=10
k=3
frame=$1 #Please specify total frame number

##############################################Initial value check###############################################
#Please specify topology file as $Top

EXCUDIR=$PWD
TrajDIR=/home/je714/Troponin/IAN_Troponin/completehowarthcut/salted/IAN30Asalt/Everything
#Directory where the Trajs are put in
Top=/home/je714/Troponin/IAN_Troponin/completehowarthcut/salted/IAN30Asalt/nowatjuancorrect30nohmr.prmtop
#SPECIFY THE .prmtop file

  echo -e "\nPlease check !!\nWorking directory: $EXCUDIR \nTopplogy: $Top \nTrajectory Path: $TrajDIR"

noise=$2 #noise population
mask=$3
j=`echo "$frame/$offset" | bc` #Total population

  if [[ -z "$frame" ]]; then
    echo -e "\tWARNING: total frame number is NULL !! \nPlease use script as ./cluster.sh <Total frames> <noise population> <res/atom mask> !!"
    exit
  else
    echo "Noise percentage set to `echo  "100*$noise/$j" | bc` % with a total population ${j}."
  fi

  if [[ -z "$noise" ]]; then
    echo -e "\tWARNING: Noise population is NULL !! \nPlease use script as ./cluster.sh <Total frames> <noise population> <res/atom mask> !!"
    exit
  else
    echo "Noise percentage set to `echo  "100*$noise/$j" | bc` % with a total population ${j}."
  fi

  if [[ -z "$mask" ]]; then
    echo -e "\tWARNING: Residue/atom mask is NULL !! \nMask set to default: @CA !!"
    mask=@CA
    echo -e "\tMask set to default: @CA !!"
  fi

  echo -e "Please check set up values. \n\tSieve number set to ${sieve}. \n\tTrajectory offset set to $offset. \n\tTotal population $j. \n\tMask: $mask (cpptraj format) \nk value: $k\n"

#################################################Kdist plot##########################################################

for i in `ls $TrajDIR`;
do

  cd $EXCUDIR
    if [[ ! -d ./cluster_${i} ]]; then
      echo -e "\nMaking cluster directory ..."
      mkdir ./cluster_${i}
    else 
      echo -e "\nCluster directory already exists !! Clearing ..."
      rm  -rf ./cluster_${i}/*
    fi

  cd ./cluster_${i}

  echo -e "\tPerforming $k.th kdist plotting for ${i} ... (trajectory offset set to $offset)"

  cat > kdist.in << EOF
trajin ../${i} 1 $frame $offset
rmsd first @CA
cluster dbscan rms $mask mass nofit kdist $k sieve ${sieve} loadpairdist pairdist CpptrajPairDist
EOF

  cpptraj $Top < kdist.in > kdist.o

  echo "Done $k.th kdist analysis for trajectory ${i} !! Extracting DBSCAN parameters ..."

##################################################Clustering#########################################################

m=`echo "$noise/${sieve}-1" | bc`

epsilon=`grep \ ${m}\  $EXCUDIR/cluster_${i}/Kdist.$k.dat | awk '{print $2}'`

  echo -e "\n\tPerforming cluster analysis for ${i} ..."

  cat > cluster.in << EOF
trajin ../${i} 1 $frame $offset
rmsd first @CA
cluster dbscan rms $mask mass nofit minpoints ${noise} epsilon ${epsilon} out cluster.dat summary cluster_summary.dat clusterout ${i} cpopvtime pop_vs_time.dat normframe
EOF

  cpptraj $Top < cluster.in > cluster.o

 echo -e "Done cluster analysis for trajectory ${i} !!"

done