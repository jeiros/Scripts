#!/bin/bash

# Use script as:
# $ ./strip_autoimage.sh run1 050-100 CTnI_hmr S1P



run=$1 		  # Can be run1 run2 run3 run4 ...
sim=$2		  # Can be 000-050, 050-100, 100-150 ...
cluster=$3        # Can be CTnI_hmr CTnI_runs CTnT_hmr CTnT_runs
phosphotype=$4    # Can be S1P or SEP


# Check all the inputs are there
if [[ $# -ne 4 ]]; then
	printf "Arguments missing\nExiting...\n"
	exit 1
fi

WORKDIR=$PWD
printf "\nCurrent directory is %s\n" $WORKDIR

mkdir -p ${WORKDIR}/${cluster}/${run}/${phosphotype}/${sim}
DESTINATION=${WORKDIR}/${cluster}/${run}/${phosphotype}/${sim}
cd ${DESTINATION}
scp je714@login.cx1.hpc.ic.ac.uk:/work/je714/phosphoHMR/${cluster}/${run}/${phosphotype}/results/${run}_${cluster}_${sim}.tgz .
printf "The destination directory is %s\n\n" $DESTINATION

# if [ -d "$DESTINATION" ]; then 
# 	printf "The directory exists\n\n"
# else
# 	printf "The destination directory does not exist\n"
# 	printf "\nPlease check the arguments
# --------------------------------
# Run:\t\t\t%s\nSim:\t\t\t%s\nCluster:\t\t%s\nPhosphorylation type:\t%s
# --------------------------------\n" ${run} ${sim} ${cluster} ${phosphotype}
# 	printf "Exiting now...\n\n"
# 	exit 1
# fi

# hard_drive=/Volumes/Seagate_MD/completehowarthcut/phospho/hmr_runs/${cluster}/${run}/${phosphotype}/

# if [[ -d  "$hard_drive" ]]; then
# 	printf "The hard drive is mounted\n"
# 	printf "The destination directory in the hard drive exists\n"
# else
# 	printf "The hard drive is not mounted\n"
# 	printf "Please mount it and rerun script"
# 	printf "Exiting now...\n\n"
# 	exit 1
# fi


cat << "EOF"

 _____  _   _  _____  _____  _____  _____  _____  _ 
/  ___|| | | |/  __ \/  __ \|  ___|/  ___|/  ___|| |
\ `--. | | | || /  \/| /  \/| |__  \ `--. \ `--. | |
 `--. \| | | || |    | |    |  __|  `--. \ `--. \| |
/\__/ /| |_| || \__/\| \__/\| |___ /\__/ //\__/ /|_|
\____/  \___/  \____/ \____/\____/ \____/ \____/ (_)
                                                    
                                                    
EOF


printf "\n\nUntaring the simulation\nRun:\t\t\t%s\nSim:\t\t\t%s\nCluster:\t\t%s\nPhosphorylation type:\t%s\n" ${run} ${sim} ${cluster} ${phosphotype}
tar zxvf ${run}_${cluster}_${sim}.tgz



cd $WORKDIR
printf "\n\nStripping waters out...\n\n"
printf "\n\nAutoimaging...\n\n"
cpptraj <<- EOF
	parm ${DESTINATION}/repstr.c0_phos${phosphotype}_watsalthmr.prmtop
	trajin ${DESTINATION}/05_Prod_${cluster}.phos${phosphotype}.${sim}_${run}.nc
	strip :WAT
	autoimage
	trajout ./${cluster}/${run}/${phosphotype}/05_Prod_${cluster}_${sim}ns_${run}.nc
	run
EOF

# printf "Move the restart file to the Hard Drive"
# scp ${DESTINATION}/*.rst /Volumes/Seagate_MD/completehowarthcut/phospho/hmr_runs/${cluster}/${run}/${phosphotype}/


cd ${DESTINATION}/
printf "\n\nErasing everything that was untarred\n\n"
ls | grep -v ${run}_${cluster}_${sim}.tgz | xargs rm  #  Erase all files except the tarred trajectory

cd $WORKDIR
# printf "\n\nMoving everything to the hard drive\n\n"
# # Move the tared file inside (and its containing directory) to the hard drive
# mv ${DESTINATION}/ /Volumes/Seagate_MD/completehowarthcut/phospho/hmr_runs/${cluster}/${run}/${phosphotype}/${sim}/

# # Copy the stripped and autoimaged trajectory to the hard drive
# scp ./${cluster}/${run}/${phosphotype}/05_Prod_${cluster}_${sim}ns_${run}.nc /Volumes/Seagate_MD/completehowarthcut/phospho/hmr_runs/${cluster}/${run}/${phosphotype}/

cat << "EOF"
______  _____  _   _  _____  _ 
|  _  \|  _  || \ | ||  ___|| |
| | | || | | ||  \| || |__  | |
| | | || | | || . ` ||  __| | |
| |/ / \ \_/ /| |\  || |___ |_|
|___/   \___/ \_| \_/\____/ (_)
                               
EOF

