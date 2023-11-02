#!/bin/bash

RunName=${1}
FAST5Folder=${2}
OutFolderName=${3}
AdapterThreshold=${4}
THREADS=${5}
SampleDescription=${6}
Guppy_demultiplexed=${7}
Barcoding=${8}
scripts=${9}
MID=${10}


##Additional variables
RandSize=100
FastqFolderSize=19

### combine the guppy fastq files into one per sample
RemoveMinIONAdaptersGuppy(){
if [[ $Barcoding == "TRUE" ]] ; then 
	##Only copy barcode folders with more than 19 reads
	ReadCount=$(cat ${Guppy_demultiplexed}/${MID}/*.fastq | grep "^@" | wc -l)
	printf $ReadCount
	if (($ReadCount > $FastqFolderSize)); then
		##add SampleName to folder structure
		SampleName=$(cat $SampleDescription | grep "^$MID" | awk -F";" '{print $2}')
		printf $SampleName
		mkdir -p $OutFolderName/$MID.$SampleName
		printf "\nWorking on $MID...\n"
		cat ${Guppy_demultiplexed}/${MID}/*.fastq > $OutFolderName/$MID.$SampleName/$MID.$SampleName.adapter_trim.fastq

	fi
else
	ReadCount=$(cat ${Guppy_demultiplexed}/pass/*.fastq | grep "^@" | wc -l)
	printf $ReadCount
		if (($ReadCount > $FastqFolderSize)); then
		##add SampleName to folder structure
		SampleName=$(cat $SampleDescription | grep -v "#" | awk -F";" '{print $2}' | tr -d '\n')
		printf $SampleName
		mkdir -p $OutFolderName/$SampleName
		printf "\nWorking on $SampleName...\n"
		cat ${Guppy_demultiplexed}/pass/*.fastq > $OutFolderName/$MID.$SampleName/$MID.$SampleName.adapter_trim.fastq
	fi
fi
}

#Run
RemoveMinIONAdaptersGuppy
