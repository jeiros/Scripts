#PBS -lselect=1:ncpus=1:ngpus=1:mem=2000mb:host=cx1-51-3-1
#PBS -lwalltime=192:00:0
#PBS -q pqigould

# This shell script is specific to running jobs on GTX 980's which are on host cx1-51-3-2
# Here we load all the libraries needed for the GPU/CUDA
module load cuda/6.5.19

# count is a really simple way to have the job use a sensible naming structure you need to set a value of count yourself


phosptype=S1P
prevsim=000-050
sim=050-100
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
cp /work/je714/phosphoHMR/${cluster}/${run}/${phosptype}/*.rst .
cp /work/je714/phosphoHMR/${cluster}/${run}/${phosptype}/*.prmtop .



# This job is a "production" run and so you need to set the write frequency of the print, restart and trajectory to appropriate values. The pbsexec command 
# means that it will stop the calculation within the last 15 minutes of the 192 hours if it has not already finished so that you do not lose the output
#
#
pbsexec -grace 15 /home/igould/pmemd.cuda_SPFP -O -i 05_Prod.in -o 05_Prod_${cluster}.phos${phosptype}.${sim}_${run}.out -c 05_Prod_${cluster}.phos${phosptype}.${prevsim}_${run}.rst -p $prmtop -r 05_Prod_${cluster}.phos${phosptype}.${sim}_${run}.rst -x 05_Prod_${cluster}.phos${phosptype}.${sim}_${run}.nc

#
# tar up everything in the temporary directory
cp /tmp/pbs.$PBS_JOBID/05_Prod_${cluster}.phos${phosptype}.${sim}_${run}.rst /work/je714/phosphoHMR/${cluster}/${run}/${phosptype}/
rm /tmp/pbs.$PBS_JOBID/05_Prod_${cluster}.phos${phosptype}.${prevsim}_${run}.rst
tar -zcvf /tmp/pbs.$PBS_JOBID/${run}_${cluster}_${sim}.tgz *
#
#
# Copy output back to a results directoy in the submission directory you need to create the "results" directory before you submit the job. 
cp *tgz /work/je714/phosphoHMR/${cluster}/${run}/${phosptype}/results/

#
#

