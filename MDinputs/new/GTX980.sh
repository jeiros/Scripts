#PBS -lselect=1:ncpus=1:ngpus=1:mem=2000mb:host=cx1-51-3-1
#PBS -lwalltime=192:00:0
#PBS -q pqigould

# This shell script is specific to running jobs on GTX 980's which are on host cx1-51-3-1 or cx1-51-6-1
# Here we load all the libraries needed for the GPU/CUDA
module load cuda/6.5.19

# count is a really simple way to have the job use a sensible naming structure you need to set a value of count yourself

phosptype=S1P
sim=000-050
run=run1
cluster=CTnI_hmr

count=repstr.c0_phos${phosptype}_watsalt
prmtop=${count}hmr.prmtop
inpcrd=${count}.inpcrd

#This is just a bit of housekeeping
echo Working directory is $PBS_O_WORKDIR

#The job is run on a "backend" server with GPU's so first thing to do is to cd into the PBS temporary directory where the job is run from 
cd /tmp/pbs.$PBS_JOBID

#This block copies all the required input and control files that are neeeded to start a set of calculations to the appropriate directory
cp /work/je714/phosphoHMR/${cluster}/${run}/${phosptype}/*.in .
cp /work/je714/phosphoHMR/${cluster}/${run}/${phosptype}/*.inpcrd .
cp /work/je714/phosphoHMR/${cluster}/${run}/${phosptype}/*.prmtop .

# first calculation performs an initial minimisation to clean up bad contacts if you have a solvent box around your protein/solute etc
# it will be quite quick of the order of minutes if has harmonic retraints on protein/solute so that only counterions and water minimise
#
/home/igould/pmemd.cuda_SPFP -O -i premin.in -o premin_${count}_${cluster}.out -c $inpcrd -p $prmtop -r premin_${count}_${cluster}.rst  -ref $inpcrd
#
# The next job is to perform a short minimisation on the whole system with no constraints again this should be of the order of minutes
/home/igould/pmemd.cuda_SPFP -O -i sandermin1.in -o sandermin_${count}_${cluster}.out -c premin_${count}_${cluster}.rst -p $prmtop -r sandermin1_${count}_${cluster}.rst
#
# The next job performs a NVE heating up of the system with restraints on the protein/solute to 100K it does this in a small number of timesteps and 
# should usually be of the order of 10 minutes
#
/home/igould/pmemd.cuda_SPFP -O -i 02_Heat.in -o 02_Heat_${count}_${cluster}.out -c sandermin1_${count}_${cluster}.rst -p $prmtop -r 02_Heat_${count}_${cluster}.rst -x 02_Heat_${count}_${cluster}.nc -ref sandermin1_${count}_${cluster}.rst
#
# This job will "heat up" the system to the target temp of 300K, you may need to adjust your desired temp, this is performed as an NPT simulation with weak 
# restraints on the protein/solute usually takes of the order of an hour or two to perform
#
/home/igould/pmemd.cuda_SPFP -O -i 03_Heat2.in -o 03_Heat2_${count}_${cluster}.out -c 02_Heat_${count}_${cluster}.rst -p $prmtop -r 03_Heat2_${count}_${cluster}.rst -x 03_Heat2_${count}_${cluster}.nc -ref 02_Heat_${count}_${cluster}.rst
#
# This job is a "production" run and so you need to set the write frequency of the print, restart and trajectory to appropriate values. The pbsexec command 
# means that it will stop the calculation within the last 15 minutes of the 192 hours if it has not already finished so that you do not lose the output
#
#
pbsexec -grace 15 /home/igould/pmemd.cuda_SPFP -O -i 05_Prod.in -o 05_Prod_${cluster}.phos${phosptype}.${sim}_${run}.out -c 03_Heat2_${count}_${cluster}.rst -p $prmtop -r 05_Prod_${cluster}.phos${phosptype}.${sim}_${run}.rst -x 05_Prod_${cluster}.phos${phosptype}.${sim}_${run}.nc
#
#
#
# tar up everything in the temporary directory
tar -zcvf /tmp/pbs.$PBS_JOBID/${run}_${cluster}_${sim}.tgz *
#
#
# Copy output back to a results directoy in the submission directory you need to create the "results" directory before you submit the job. 
cp *tgz /work/je714/phosphoHMR/${cluster}/${run}/${phosptype}/results/
#
#
