#!/bin/bash
settingsfile=${1}
#settingsfile="settingsfile-zooplankton.txt"
RunAnalysisScript="Run_DNA-metabarcoding-MinION_recluster_without_settings.sh"
dos2unix $settingsfile
echo $settingsfile
RunName=$(grep '^RunName=' $settingsfile | awk -F '=' '{print $2}'| tr -d '"' | tr -d "'")
GPU=$(grep '^GPU=' $settingsfile | awk -F '=' '{print $2}'| tr -d '"' | tr -d "'")
OutDir=$(grep '^OutDir=' $settingsfile | awk -F '=' '{print $2}'| tr -d '"' | tr -d "'")
THREADS=$(grep '^THREADS=' $settingsfile | awk -F '=' '{print $2}'| tr -d '"' | tr -d "'")
RunModules=$(grep '^RunModules=' $settingsfile | awk -F '=' '{print $2}'| tr -d '"' | tr -d "'")
mail_user=$(grep '^mail_user=' $settingsfile | awk -F '=' '{print $2}'| tr -d '"' | tr -d "'")


echo $RunName
cat $settingsfile | cat - $RunAnalysisScript > Run_DNA-metabarcoding-MINION.$RunName.sh

### Create folder for LOGS of each sample and create folder for job_scheduler scripts
mkdir -p $OutDir/$RunName/job_scheduler/LOGS/
if ! command -v sinfo &> /dev/null
then
        work_env='local'
else
        work_env='HPC'
fi


function job_scheduler_script () {
echo "Generating job scheduler script for DNA metabarcoding..."
if  [[ "$work_env" == "HPC" ]] ; then
        echo "results will be sent to:  $mail_user"
fi
jobScheduler="${OutDir}/jobScheduler.${RunName}.sh"
printf "#!/bin/bash \n" > "$jobScheduler"
printf "#SBATCH --comment=RIKILT \n" >> "$jobScheduler"
printf "#SBATCH --mincpus=$THREADS\n" >> "$jobScheduler"
printf "#SBATCH --mem-per-cpu=8000 \n" >> "$jobScheduler"
printf "#SBATCH --time=200 \n" >> "$jobScheduler"
printf "#SBATCH --ntasks=1 \n" >> "$jobScheduler"
printf %s '#SBATCH --output='"$OutDir/$RunName"'/job_scheduler/LOGS/output_%j.txt '>> "$jobScheduler"
printf "\n" >> "$jobScheduler"
printf %s '#SBATCH --error='"$OutDir/$RunName"'/job_scheduler/LOGS/error_output_%j.txt ' >> "$jobScheduler"
printf "\n" >> "$jobScheduler"
printf "#SBATCH --job-name=Metabarcoding_${RunName} \n" >> "$jobScheduler"
printf "#SBATCH --mail-type=FAIL \n" >> "$jobScheduler"
printf "#SBATCH --mail-user=$mail_user \n\n" >> "$jobScheduler"
### Run program per sample
printf "time bash Run_DNA-metabarcoding-MINION.$RunName.sh  \n" >> "$jobScheduler"

if  [[ "$work_env" == "HPC" ]] ; then
sbatch $jobScheduler
else
bash $jobScheduler
fi
}

job_scheduler_script
