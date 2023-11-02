#!/bin/bash
settingsfile=${1}
#settingsfile="settingsfile-zooplankton.txt"
RunAnalysisScript="Run_DNA-metabarcoding-MinION_without_settings.sh"
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
if ! command -v sinfo &> /dev/null
then
	work_env='local'
else
	work_env='HPC'
fi
mkdir -p $OutDir
mkdir -p $OutDir/../LOGS
function job_scheduler_script () {
echo "Generating job scheduler script for DNA metabarcoding..."
if  [[ "$work_env" == "HPC" ]] ; then
	echo "results will be sent to:	$mail_user"
fi
jobScheduler="${OutDir}/jobScheduler.${RunName}.sh"
printf "#!/bin/bash \n" > "$jobScheduler"
printf "#SBATCH --comment=RIKILT \n" >> "$jobScheduler"
printf "#SBATCH --mincpus=$THREADS\n" >> "$jobScheduler"
printf "#SBATCH --mem-per-cpu=4000 \n" >> "$jobScheduler"
if [[ $GPU = 1 ]] ; then
	printf "#SBATCH --time=500 \n" >> "$jobScheduler"
	printf "#SBATCH --gres=gpu:1 \n" >> "$jobScheduler"
	printf "#SBATCH --partition=gpu \n" >> "$jobScheduler"
else
	printf "#SBATCH --time=2000 \n" >> "$jobScheduler"
fi
printf "#SBATCH --ntasks=1 \n" >> "$jobScheduler"
printf %s '#SBATCH --output='"$OutDir"'/../LOGS/output_%j.txt '>> "$jobScheduler"
printf "\n" >> "$jobScheduler"
printf %s '#SBATCH --error='"$OutDir"'/../LOGS/error_output_%j.txt ' >> "$jobScheduler"
printf "\n" >> "$jobScheduler"
printf "#SBATCH --job-name=Metabarcoding_${RunName} \n" >> "$jobScheduler"
printf "#SBATCH --mail-type=FAIL \n" >> "$jobScheduler"
printf "#SBATCH --mail-user=$mail_user \n\n" >> "$jobScheduler"
### Run program per sample
printf "bash Run_DNA-metabarcoding-MINION.$RunName.sh  \n" >> "$jobScheduler"
if ([[ "$RunModules" == *"all"* ]] || [[ "$RunModules" == *"uppy"* ]] || [[ "$RunModules" == *"oretools"* ]]) && [[ "$work_env" == "HPC" ]] ; then
	sbatch $jobScheduler
else
	bash $jobScheduler
fi
}

job_scheduler_script
