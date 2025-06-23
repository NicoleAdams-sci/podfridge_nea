#!/bin/bash
#SBATCH --job-name=STR_sims
#SBATCH --output=/home/%u/PODFRIDGE/slurm/%x-%j.log
#SBATCH --time=15:00
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4g
#SBATCH --array=1-10
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=$(whoami)@umich.edu

########################################################################
# Reminder: Please set your GitHub PAT using the following command     #
# before running this script:                                          #
# export GITHUB_PAT='your_personal_access_token'                       #
# export GITHUB_USER='your_github_username'                            #
########################################################################

# Get the unique username
UNIQNAME=$(whoami)

# Inputs
RELATED=$1
UNRELATED=$2

echo "Starting job $SLURM_JOB_ID, array task $SLURM_ARRAY_TASK_ID"
echo "Inputting $RELATED related and $UNRELATED unrelated pairs"

# Ensure the directory structure exists
mkdir -p /home/$UNIQNAME/PODFRIDGE/slurm
mkdir -p /home/$UNIQNAME/PODFRIDGE/logfiles

# Change to the project directory
cd /home/$UNIQNAME/PODFRIDGE

# Activate conda
#CONDA_HOME="/home/$UNIQNAME/miniconda3"
#source "$CONDA_HOME/etc/profile.d/conda.sh"
source ~/.bashrc

# Activate the environment
conda activate rstats

# Log memory and CPU usage every 15 minutes
#log_resource_usage() {
#    squeue --job=$SLURM_JOB_ID --format="%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R %C %m %N" | tee -a /home/$UNIQNAME/PODFRIDGE/slurm/resource_usage_$SLURM_JOB_ID.log
#}

# Run the logging function every 15 minutes in the background
#while true; do
#    log_resource_usage
#    sleep 900  # Pause for 15 minutes
#done &

# Run the R script with command line arguments
#Rscript code/STR_sims_NEA.R $RELATED $UNRELATED
Rscript code/STR_sims_allPopLR_NEA.R $RELATED $UNRELATED

# Configure Git to use HTTPS and PAT
#git remote set-url origin https://github.com/lasisilab/PODFRIDGE.git

# Log memory and CPU usage after job completion
#sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,AllocCPUs,Elapsed,MaxRSS,MaxVMSize,State > /home/$UNIQNAME/PODFRIDGE/slurm/usage_$SLURM_JOB_ID.log

# Commit and push changes to GitHub
#git add .
#git commit -m "Automated commit of new results $(date)"
#echo $GITHUB_PAT | git push https://$GITHUB_USER:$GITHUB_PAT@github.com/lasisilab/PODFRIDGE.git
