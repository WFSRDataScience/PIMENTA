#!/bin/bash

##Modified 5 September 2023, V. van der Vorst

### Input variables
RunName=${1}
FAST5Folder=${2}
OutDir=${3}
Barcoding=${4}
MinionKit=${5}
MinionFlowCell=${6}
AdapterThreshold=${7}
THREADS=${8}
SampleDescription=${9}
GPU=${10}
Guppy_demultiplexed=${11}
ExpansionKit=${12}
mail_user=${13}

printf "$AdapterThreshold"

printf "$GPU"

if [ "$GPU" = "1" ] ; then
	module load  SHARED/.guppy/6.5.7
#	nvidia-smi
#	modinfo nvidia
#	cat /proc/driver/nvidia/version
elif [ "$GPU" = "cpu" ]  ; then
	module load SHARED/.guppy/6.5.7-CPU
else
	echo 'ERROR: Please specify if you want to run Guppy using the GPU or CPU!'
	echo "$GPU"
fi

##Additional variables
OutFolderName="$OutDir/$RunName"


### Convert FAST5 into FASTQ format using ONT Guppy
function PrepareFastqFilesGuppy(){
	printf "\nRunning ONT Guppy...\n"
	### Perform basecalling
	if [[ "$GPU" == '1' ]]; then
		guppy_basecaller -i $FAST5Folder/ -x 'cuda:0' --flowcell $MinionFlowCell --kit $MinionKit --recursive --chunk_size 3500 -s $FAST5Folder/Guppy
	else
		guppy_basecaller -i $FAST5Folder --flowcell $MinionFlowCell --kit $MinionKit --recursive --cpu_threads_per_caller 1 --chunk_size 3500 --num_callers $THREADS -s $FAST5Folder/Guppy
	fi
	printf "Guppy basecaller finished"
    if [[ $Barcoding == "TRUE" ]] ; then
	if [[ "$GPU" == '1' ]]; then
		if [[ $ExpansionKit == "none" ]] ; then
                       guppy_barcoder --barcode_kits "${MinionKit}" -x 'cuda:0' -i $FAST5Folder/Guppy/pass -t $THREADS -s $Guppy_demultiplexed --enable_trim_barcodes
 		else
			guppy_barcoder --barcode_kits "${ExpansionKit}" -x 'cuda:0' -i $FAST5Folder/Guppy/pass -t $THREADS -s $Guppy_demultiplexed --enable_trim_barcodes
		fi
	else
		if [[ $ExpansionKit == "none" ]] ; then
			guppy_barcoder --barcode_kits "${MinionKit}" -i $FAST5Folder/Guppy/pass -t $THREADS -s $Guppy_demultiplexed --enable_trim_barcodes
		else
			guppy_barcoder --barcode_kits "${ExpansionKit}" -i $FAST5Folder/Guppy/pass -t $THREADS -s $Guppy_demultiplexed --enable_trim_barcodes
		fi
	fi
	else
		echo "skipping demultiplexing"
	fi
	}


####Run workflow
PrepareFastqFilesGuppy
