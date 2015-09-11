#!/bin/bash

# Use script as:
# $ ./strip_autoimage.sh run1 050-100 CTnI_hmr S1P



run=$1 		  # Can be run1 run2 run3 run4 ...
sim=$2		  # Can be 000-050, 050-100, 100-150 ...
cluster=$3        # Can be CTnI_hmr CTnI_runs CTnT_hmr CTnT_runs
phosphotype=$4    # Can be S1P or SEP

if [[ -z ${run+x} ]]; then
	printf "Variable %s is missing" ${run}
	printf "Exiting now... \n\n"
	exit 1
else
	:
fi

if [[ -z ${sim+x} ]]; then
	printf "Variable %s is missing" ${sim}
	printf "Exiting now... \n\n"
	exit 1
else
	:
fi

if [[ -z ${cluster+x} ]]; then
	printf "Variable %s is missing" ${cluster}
	printf "Exiting now... \n\n"
	exit 1
else
	:
fi

if [[ -z ${phosphotype+x} ]]; then
	printf "Variable %s is missing" ${phosphotype}
	printf "Exiting now... \n\n"
	exit 1
else
	:
fi

WORKDIR=$PWD
printf "\nCurrent directory is %s\n" $WORKDIR
DESTINATION=${WORKDIR}/${cluster}/${run}/${phosphotype}/${sim}
printf "The destination directory is %s\n\n" $DESTINATION

if [ -d "$DESTINATION" ]; then 
	printf "The directory exists\n\n"
else
	printf "The destination directory does not exist\n"
	printf "\nPlease check the arguments
--------------------------------
Run:\t\t\t%s\nSim:\t\t\t%s\nCluster:\t\t%s\nPhosphorylation type:\t%s
--------------------------------\n" ${run} ${sim} ${cluster} ${phosphotype}
	printf "Exiting now...\n\n"
	exit 1
fi

hard_drive=/mnt/ntfs

if [[ -d  "$hard_drive" ]]; then
	printf "The hard drive is mounted\n\n"
else
	printf "The hard drive is not mounted\n"
	printf "Please mount it and rerun script"
	printf "Exiting now...\n\n"
	exit 1
fi


cat << "EOF"

 _____  _   _  _____  _____  _____  _____  _____  _ 
/  ___|| | | |/  __ \/  __ \|  ___|/  ___|/  ___|| |
\ `--. | | | || /  \/| /  \/| |__  \ `--. \ `--. | |
 `--. \| | | || |    | |    |  __|  `--. \ `--. \| |
/\__/ /| |_| || \__/\| \__/\| |___ /\__/ //\__/ /|_|
\____/  \___/  \____/ \____/\____/ \____/ \____/ (_)
                                                    
                                                    
EOF


cd ${DESTINATION}/
printf "\n\nUntaring the simulation\nRun:\t\t\t%s\nSim:\t\t\t%s\nCluster:\t\t%s\nPhosphorylation type:\t%s\n" ${run} ${sim} ${cluster} ${phosphotype}
tar zxvf ${run}_${cluster}_${sim}.tgz



cd $WORKDIR
printf "\n\nStripping waters out...\n\n"
cpptraj <<- EOF
	parm /home/je714/Troponin/IAN_Troponin/completehowarthcut/phospho/structures/${cluster}/repstr.c0_phos${phosphotype}_watsalthmr.prmtop
	trajin ${DESTINATION}/05_Prod_${cluster}.phos${phosphotype}.${sim}_${run}.nc
	strip :WAT
	trajout ${DESTINATION}/${sim}.nc
	run
EOF
printf "\n\nAutoimaging...\n\n"

cpptraj <<- EOF
	parm ./${cluster}/${run}/${phosphotype}/repstr.c0_phos${phosphotype}_nowat.prmtop
	trajin ${DESTINATION}/${sim}.nc
	autoimage
	trajout ./${cluster}/${run}/${phosphotype}/05_Prod_${cluster}_${sim}ns_${run}.nc
	run
EOF


cd ${DESTINATION}/
printf "\n\nErasing everything that was untarred\n\n"
find ! -name '*.tgz' -type f -exec rm -f {} +


cd $WORKDIR
printf "\n\nMoving everything to the hard drive\n\n"
mv ${DESTINATION}/ /mnt/ntfs/completehowarthcut/phospho/hmr_runs/${cluster}/${run}/${phosphotype}/${sim}/
cp ./${cluster}/${run}/${phosphotype}/05_Prod_${cluster}_${sim}ns_${run}.nc /mnt/ntfs/completehowarthcut/phospho/hmr_runs/${cluster}/${run}/${phosphotype}/

cat << "EOF"
______  _____  _   _  _____  _ 
|  _  \|  _  || \ | ||  ___|| |
| | | || | | ||  \| || |__  | |
| | | || | | || . ` ||  __| | |
| |/ / \ \_/ /| |\  || |___ |_|
|___/   \___/ \_| \_/\____/ (_)
                               
EOF

