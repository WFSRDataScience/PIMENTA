#conda activate pimenta
# Run script
## Modified 28 February 2023, V. van der Vorst
if [ "$FAST5Folder" = "" ] || [ "$SampleDescription" = "" ]; then
	echo "something went wrong, check the settings"
	printf "FAST5 folder: $FAST5Folder \nSample description: $SampleDescription"
	exit 1
fi

if [[ $THREADS -eq '' ]]
  then
    THREADS=15
fi

if [ "$SampleDescription" = "" ]
  then
    SampleDescription="None"
  fi

#retrieve taxids from excluded taxid list
cat $scripts/Excluded.NCBI.taxids.tsv |sed 1d |  awk -F ';' '{print $1}' > $ExTaxids

OutFolderName="$OutDir/$RunName"
mkdir -p $OutFolderName

function UserInput(){
DATE=$(date)
#### Generating settingsfile
printf "Date of analysis: $DATE \n" > $OutFolderName/$RunName.settings.all.txt
printf "Run name:  $RunName \n" >> $OutFolderName/$RunName.settings.all.txt
printf "FAST5 Input folder:  $FAST5Folder \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Output folder:  $OutFolderName \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Output taxonomy folder:  ${OutFolderName}/Taxonomy \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Minion Sequencing kit:  $MinionKit \n" >> $OutFolderName/$RunName.settings.all.txt
printf "MinionFlowCell:  $MinionFlowCell \n" >> $OutFolderName/$RunName.settings.all.txt
printf "ONT ExpansionKit: $ExpansionKit \n" >> $OutFolderName/$RunName.settings.all.txt
printf "PoreChop setting - adapter_threshold: $AdapterThreshold \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Prinseq setting -range_len : $MinMaxLength \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Prinseq setting -min_qual_mean : $MinQualMean \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Prinseq setting -trim_left : $TrimLeft \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Prinseq setting -trim_right : $TrimRight \n" >> $OutFolderName/$RunName.settings.all.txt
printf "CD-HIT-est setting -c : $Ident1 \n" >> $OutFolderName/$RunName.settings.all.txt
printf "CD-HIT-est filtering -minclustsize : $MinClustSize \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Cutadapt setting --error-rate : $Error \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Blast setting -evalue : $Evalue  \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Blast setting -max_target_seqs : $MaxTargetSeqs  \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Blast filtering -qcovs : $Qcov \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Blast filtering -pident : $Pident \n" >> $OutFolderName/$RunName.settings.all.txt
printf "Number of CPU threads: $THREADS \n" >> $OutFolderName/$RunName.settings.all.txt
if [[ "$GPU" == 'cpu' ]]; then
	printf "Number of GPUs: $GPU \n">> $OutFolderName/$RunName.settings.all.txt
fi
}

RUN1(){

### Variable settings
Guppy_demultiplexed="${FAST5Folder}/Guppy_demultiplexed"
echo $Guppy_demultiplexed
Guppy_Folder="${FAST5Folder}/Guppy"
JobScheduler="$OutFolderName/job_scheduler/job_scheduler.sh"

## Poretools checks the quality of the fast5 data, only runs when specified in $RunModules
if [[ "$RunModules" == *"oretools"* ]] ; then
	time bash $scripts/RunPoretools.sh $OutFolderName $FAST5Folder $RunName $THREADS &
fi


### check if GPU argument is given, if not, don't run Guppy
if [[ "$RunModules" == *"all"* ]] || [[ "$RunModules" == *"uppy"* ]] ; then
	start=$SECONDS
	time bash $scripts/ReadCalling_Demultiplexing_HPC.sh $RunName $FAST5Folder $OutDir $Barcoding $MinionKit $MinionFlowCell $THREADS $SampleDescription $GPU $Guppy_demultiplexed $ExpansionKit $mail_user
	duration=$(( SECONDS - start ))
fi

### Create folder for LOGS of each sample and create folder for job_scheduler scripts
mkdir -p $OutFolderName/job_scheduler/LOGS/
## Remove old Taxonomy results
rm -r $OutFolderName/Taxonomy

MIDfolderCount=$(find $Guppy_demultiplexed/ -maxdepth 1 -name "barcode*" -type d | wc -l)
echo $MIDfolderCount
#Check if there is slurm in the working env, if so run in parallel.
if ! command -v sinfo &> /dev/null
        then
        echo "slurm not found, running workflow consecutively"
	for SLURM_ARRAY_TASK_ID in $(eval echo {1..$MIDfolderCount}); do
	MIDfolder=$(find "$Guppy_demultiplexed" -maxdepth 1 -name "barcode*" -type d | sort | head -n $SLURM_ARRAY_TASK_ID  | tail -n 1)
	echo $MIDfolder
	MID=$(basename $MIDfolder)
	echo "Task ID: $SLURM_ARRAY_TASK_ID \t Barcode: $MID"
	time bash $scripts/Pipeline_HPC.sh $RunName $OutDir $OutFolderName $MinMaxLength $TrimLeft $TrimRight $MinQualMean $Ident1 $MinClustSize $SampleDescription $THREADS $RunModules $scripts $MID
	done
else
#Creating job_scheduler script for slurm
printf "#!/bin/bash \n" > $JobScheduler
printf "#SBATCH --comment=Metabarcoding \n" >> $JobScheduler
printf "#SBATCH --time=$timepersample \n" >> $JobScheduler
printf "#SBATCH --mem-per-cpu=8000 \n" >> $JobScheduler
printf "#SBATCH --ntasks=1 \n" >> $JobScheduler
printf "#SBATCH --tmp=20 \n" >> $JobScheduler
printf "#SBATCH --job-name=Metabarcoding_$RunName \n" >> $JobScheduler
printf "#SBATCH --partition=main \n" >> $JobScheduler
printf "#SBATCH --mail-type=FAIL \n" >> $JobScheduler
printf "#SBATCH --mail-user=$mail_user\n" >> $JobScheduler
if [[ "$MIDfolderCount" > 1 ]]; then
	printf %s '#SBATCH --array=1-'"${MIDfolderCount}"'%'"${MIDfolderCount}" >> $JobScheduler
	printf "\n#SBATCH --cpus-per-task=$THREADS \n" >> $JobScheduler
	printf %s '#SBATCH --output='"$OutFolderName"'/job_scheduler/LOGS/output_%A_%a.txt '>> $JobScheduler
	printf "\n" >> $JobScheduler
	printf %s '#SBATCH --error='"$OutFolderName"'/job_scheduler/LOGS/error_output_%A_%a.txt ' >> $JobScheduler
	printf "\n" >> $JobScheduler
else
	printf "#SBATCH --mincpus=$THREADS \n" >> $JobScheduler
	printf %s '#SBATCH --output='"$OutFolderName"'/job_scheduler/LOGS/output_%j.txt '>> $JobScheduler
	printf "\n" >> $JobScheduler
	printf %s '#SBATCH --error='"$OutFolderName"'/job_scheduler/LOGS/error_output_%j.txt ' >> $JobScheduler
	printf "\n" >> $JobScheduler
fi

### load modules
printf "\nsource ~/.bashrc\nconda activate pimenta\n" >> $JobScheduler

if [[ "$MIDfolderCount" > 1 ]]; then
	### Run program per sample
	printf 'MIDfolder=$(find '"$Guppy_demultiplexed"' -maxdepth 1 -name "barcode*" -type d | sort | head -n $SLURM_ARRAY_TASK_ID  | tail -n 1)'"\n"  >> $JobScheduler
	printf 'MID=$(basename $MIDfolder | sed '"'s/\.[^.]*$//'"')' >> $JobScheduler
	printf "\n"'echo "Task ID: $SLURM_ARRAY_TASK_ID \t Barcode: $MID"'"\n" >> $JobScheduler
else
	printf 'MID=$(cat '$SampleDescription' | grep -v "#" | awk -F\; '"'"'{print $1}'"'"')' >> $JobScheduler
	printf "\n"'echo "Barcode: $MID"'"\n" >> $JobScheduler
	if [[ $Barcoding == "FALSE" ]]; then
		Guppy_demultiplexed=$Guppy_Folder
	fi

fi

printf "time bash $scripts/Pipeline_HPC.sh $RunName $OutDir $OutFolderName $MinMaxLength $TrimLeft $TrimRight $MinQualMean $Ident1 $MinClustSize $SampleDescription $THREADS $RunModules $scripts $FastqFolderSize $Guppy_demultiplexed "'$MID'" \n" >> $JobScheduler
if [[ "$RunModules" == *"oldmode"* ]] ; then
printf "bash $scripts/SelectBarcodeBlast.sh $RunName $OutDir $OutFolderName $Ident1 $PrimerFile $Error $Evalue $Qcov $Pident $MaxTargetSeqs $SampleDescription $THREADS $RunModules $NT_dmp "'$MID'" $ExTaxids $ExSeqids \n"  >> $JobScheduler
printf "old mode\n"
fi
### running sbatch
module unload slurm
module load slurm
sbatch $JobScheduler
fi
}

##Run analysis
UserInput
RUN1

