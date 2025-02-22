# Run script
## Modified 16 april 2020, V. van der Vorst


if [[ $THREADS -eq '' ]]
  then
    THREADS=15
fi

#retrieve taxids from excluded taxid list
cat $scripts/Excluded.NCBI.taxids.tsv |sed 1d |  awk -F ';' '{print $1}' > $ExTaxids

OutFolderName="$OutDir/$RunName" 
mkdir -p $OutFolderName

#Activate conda environment
source ~/.bashrc
#conda info --envs
#. $CONDA_PREFIX/etc/profile.d/conda.sh
conda deactivate 
conda activate pimenta
module load R
if [ -z "$filtertax" ]; then
	filtertax="notaxfilter"
fi
printf "$scripts/Pipeline_HPC_recluster.sh $RunName $OutDir $OutFolderName $Ident2 $MinClustSize2 $PrimerFile $Error $Evalue $Qcov $Pident $MaxTargetSeqs $THREADS $RunModules $NT_dmp $ExTaxids $Targets $ExSeqids $modus $DATABASE $filtertax $scripts\n"

time bash $scripts/Pipeline_HPC_recluster.sh $RunName $OutDir $OutFolderName $Ident2 $MinClustSize2 $PrimerFile $Error $Evalue $Qcov $Pident $MaxTargetSeqs $THREADS $RunModules $NT_dmp $ExTaxids $Targets $ExSeqids $modus $DATABASE $filtertax $scripts


