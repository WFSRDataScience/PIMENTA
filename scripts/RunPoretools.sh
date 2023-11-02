#!/bin/bash
module load SHARED/metabarcoding/.poretools-0.6.0

#variables
OutFolderName=${1}
FAST5Folder=${2}
RunName=${3}
THREADS=${4}
## get stats for each FAST5Folder
function RunPoretools(){
	mkdir -p $OutFolderName/poretools
	if [ -d "$FAST5Folder/fast5_pass" ] && [ ! -d "$FAST5Folder/fast5/" ]; then
		multi_to_single_fast5 -i $FAST5Folder/ -s $FAST5Folder/fast5/ -t $THREADS --recursive 
	fi
	printf "Stats of FAST5 files \n" > $OutFolderName/poretools/$RunName.poretools.stats
	printf "FAST5 folder: $FAST5Folder\n" >> $OutFolderName/poretools/$RunName.poretools.stats
	total_folders=$(find $FAST5Folder/fast5/ -mindepth 1 -maxdepth 1 -type d | wc -l)
	printf "total subfolders: $total_folders\n\n" >> $OutFolderName/poretools/$RunName.poretools.stats
	if [ $total_folders -gt 1 ] ; then 
		for folder in $FAST5Folder/fast5/* ; do
			if [ -d $folder ]; then 
				basename=$(basename -- "$folder")
				echo $folder
				printf "\nfolder number: $basename\n" >> $OutFolderName/poretools/$RunName.poretools.stats 
				poretools stats $folder >> $OutFolderName/poretools/$RunName.poretools.stats
				poretools hist --saveas $OutFolderName/poretools/$RunName.$basename.poretools.hist.lengths.png $folder
			fi
		done
	else
		poretools stats $FAST5Folder/fast5/ > $OutFolderName/poretools/$RunName.poretools.stats
		poretools hist --saveas $OutFolderName/poretools/$RunName.poretools.hist.lengths.png $FAST5Folder/fast5/
	fi 
	cat $OutFolderName/poretools/$RunName.poretools.stats
}

#RUN
RunPoretools
