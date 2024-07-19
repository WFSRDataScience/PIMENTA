#!/bin/bash

##Modified 1 March 2023, V. van der Vorst


### Input variables
RunName=${1}
OutDir=${2}
OutFolderName=${3}
MinMaxLength=${4}
TrimLeft=${5}
TrimRight=${6}
MinQualMean=${7}
Ident1=${8}
MinClustSize=${9}
SampleDescription=${10}
THREADS=${11}
Runmodules=${12}
scripts=${13}
FastqFolderSize=${14}
Guppy_demultiplexed=${15}
MID=${16}

printf "$ExSeqids"
printf "Running modules: $Runmodules"
if [[ "$MID" == "" ]]; then
	exit 1
fi

##Additional variables
if [[ -f "$SampleDescription" ]]; then
SampleName=$(cat $SampleDescription | grep "^$MID" | awk -F";" '{print $2}')
else
SampleName=$MID
fi
MIDFolderName="$OutFolderName/$MID.$SampleName"
RandSize=100
nt=true #Use nt database instead of custom barcode databases

### Database variables
export BLASTDB=$BLASTDB:/lustre/shared/wfsr-databases/BLASTdb
DATABASE=/lustre/shared/wfsr-databases/BLASTdb/nt

### Create directory and print user input
function UserInput(){
	# Get the variables and create the output folder
	mkdir -p $MIDFolderName/
	mkdir -p $OutFolderName/Taxonomy/
	DATE=$(date)

	# Save the settings in a file
	printf "Species Identification using the Demultiplexed MinION Data...\n\n" > $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Date of analysiss: $DATE \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Sample or analysis name:  $MID.$SampleName \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Input file:  $MIDFolderName/$MID.adapter_trim.fastq \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Output folder:  $MIDFolderName \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Prinseq setting -range_len : $MinMaxLength \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Prinseq setting -min_qual_mean : $MinQualMean \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Prinseq setting -trim_left : $TrimLeft \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Prinseq setting -trim_right : $TrimRight \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "CD-HIT-est setting -c : $Ident1 \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "CD-HIT-est filtering -minclustsize : $MinClustSize \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Cutadapt setting --error-rate (-e) : $Error \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Blast setting -evalue : $Evalue  \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Blast setting -max_target_seqs : $MaxTargetSeqs  \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Blast filtering -qcovs : $Qcov \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Blast filtering -pident : $Pident \n" >> $MIDFolderName/$MID.$SampleName.settings.txt
	printf "Number of CPU threads: $THREADS \n">> $MIDFolderName/$MID.$SampleName.settings.txt
}

### QC trimming using Prinseq
QCtrimming(){
        #check if guppy demultiplexed folder for $MID (barcode) has enough reads, if so copy to $MID output folder
        if [ -d "${Guppy_demultiplexed}/${MID}" ]; then  
                ReadCount=$(cat ${Guppy_demultiplexed}/${MID}/*.fastq | grep "^@" | wc -l)
                printf $ReadCount
                if (($ReadCount > $FastqFolderSize)); then
                   ##add SampleName to folder structure
                   SampleName=$(cat $SampleDescription | grep "^$MID" | awk -F";" '{print $2}')
                   printf $SampleName
                   mkdir -p $OutFolderName/$MID.$SampleName
                   printf "\nWorking on $MID...\n"
                   cat ${Guppy_demultiplexed}/${MID}/*.fastq > $OutFolderName/$MID.$SampleName/$MID.$SampleName.adapter_trim.fastq
        	else
		   exit 0
                fi
        else
                SampleName=$(cat $SampleDescription | grep -v "#" | awk -F";" '{print $2}' | tr -d '\n')
                printf $SampleName
                mkdir -p $OutFolderName/$SampleName
                printf "\nWorking on $SampleName...\n"
                cat ${Guppy_demultiplexed}/pass/*.fastq > $OutFolderName/$MID.$SampleName/$MID.$SampleName.adapter_trim.fastq

        fi

	printf "\nRunning Prinseq...\n"
	prinseq-lite.pl -fastq $MIDFolderName/$MID.$SampleName.adapter_trim.fastq -trim_left $TrimLeft -trim_right $TrimRight -range_len $MinMaxLength -min_qual_mean $MinQualMean -out_good $MIDFolderName/$MID.$SampleName.QC -out_bad $MIDFolderName/$MID.$SampleName.QC.bad
	printf "prinseq-lite.pl -fastq $MIDFolderName/$MID.$SampleName.adapter_trim.fastq -trim_left $TrimLeft -trim_right $TrimRight -range_len $MinMaxLength -min_qual_mean $MinQualMean -out_good $MIDFolderName/$MID.$SampleName.QC -out_bad null\n"
	prinseq-lite.pl -fastq $MIDFolderName/$MID.$SampleName.QC.fastq -graph_stats ld,qd,ns -graph_data $MIDFolderName/$MID.$SampleName.QC.gd -out_good null -out_bad null
	printf "prinseq-lite.pl -fastq $MIDFolderName/$MID.$SampleName.QC.fastq -graph_stats ld,qd,ns -graph_data $MIDFolderName/$MID.$SampleName.QC.gd -out_good null -out_bad null\n" 
	#prinseq-graphs.pl -i $MIDFolderName/$MID.$SampleName.QC.gd -html_all -o $MIDFolderName/$MID.$SampleName.QC
	#cp $MIDFolderName/$MID.$SampleName.QC.html $OutFolderName/Taxonomy/$MID.$SampleName.QC.html
	###Convert fastq to fasta
	cat $MIDFolderName/$MID.$SampleName.QC.fastq | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' > $MIDFolderName/$MID.$SampleName.QC.fasta

}

### CD-HIT cluster>&2 echo "error"
Clustering(){
	rm -r --force $MIDFolderName/ClustCons/
	mkdir -p $MIDFolderName/ClustCons/
	##Set appropriate word size $n
	##See https://github.com/weizhongli/cdhit/wiki/3.-User%27s-Guide#clstr_sort_bypl


	Int=$(awk -v n="$Ident1" 'BEGIN{printf("%.0f\n",n*100)}')
	if [ "$Int" -ge 95 -a "$Int" -le 100 ]; then
		n=10
	elif [ "$Int" -ge 90 -a "$Int" -le 94 ]; then
		n=8
	elif [ "$Int" -ge 88 -a "$Int" -le 89 ]; then
		n=7
	elif [ "$Int" -ge 85 -a "$Int" -le 87 ]; then
		n=6
	elif [ "$Int" -ge 80 -a "$Int" -le 84 ]; then
		n=5
	elif [ "$Int" -ge 75 -a "$Int" -le 79 ]; then
		n=4
	fi

	printf "\nRunning CD-HIT-est...RUN1 at 90%% identity\n"
	cd-hit-est -i $MIDFolderName/$MID.$SampleName.QC.fasta -o $MIDFolderName/ClustCons/$MID.$SampleName.CD-hit.90 -c 0.90 -G 1 -g 1 -gap -1 -gap-ext 0 -p 1 -r 1 -d 0 -T $THREADS -n 8 -M 20000
	printf "\nRunning CD-HIT-est...RUN2 at $Int%% identity\n"
	cd-hit-est -i $MIDFolderName/ClustCons/$MID.$SampleName.CD-hit.90 -o $MIDFolderName/ClustCons/$MID.$SampleName.CD-hit.$Ident1 -c $Ident1 -G 1 -g 1 -gap -1 -gap-ext 0 -p 1 -r 1 -d 0 -T $THREADS -n $n -M 20000
	###Combine the two clstr-output files to include all original sequences from the parent clstr-file from RUN1
	clstr_rev.pl $MIDFolderName/ClustCons/$MID.$SampleName.CD-hit.90.clstr $MIDFolderName/ClustCons/$MID.$SampleName.CD-hit.$Ident1.clstr > $MIDFolderName/ClustCons/$MID.$SampleName.CD-hit.90.$Ident1.clstr
	###Renumber clusters
	clstr_renumber.pl $MIDFolderName/ClustCons/$MID.$SampleName.CD-hit.90.$Ident1.clstr > $MIDFolderName/ClustCons/$MID.$SampleName.CD-hit.90.$Ident1.renumber.clstr
	###Calculate cluster sizes of merged clstr-files.
	plot_len1.pl $MIDFolderName/ClustCons/$MID.$SampleName.CD-hit.90.$Ident1.renumber.clstr 1,2-5,6-9,10-19,20-29,30-39,40-49,50-100,101-200,201-300,301-400,401-500,501-1000,1001-5000,5001-99999 > $MIDFolderName/ClustCons/$MID.$SampleName.CD-hit.90.$Ident1.renumber.plotlen
	###Extract FASTA sequences from clusters
	printf "\nGenerating fasta files for clusters with more than $MinClustSize member sequences\n"
	make_multi_seq.pl $MIDFolderName/$MID.$SampleName.QC.fasta $MIDFolderName/ClustCons/$MID.$SampleName.CD-hit.90.$Ident1.renumber.clstr $MIDFolderName/multi-seq $MinClustSize
	mv --force $MIDFolderName/multi-seq/ $MIDFolderName/ClustCons/ #############SIMPLIFY
	find $MIDFolderName/ClustCons/ -type f -not -name "*.*" -exec mv "{}" "{}".fasta \;
}



###Consensus building
Consensus(){
printf "\nRunning MAFFT MSA and consensus building....\n"
rm $MIDFolderName/ClustCons/multi-seq/*.cons.fasta
rm $MIDFolderName/ClustCons/multi-seq/*.msf
rm $MIDFolderName/ClustCons/multi-seq/*.Rand*.fasta
### Get a consensus sequence per cluster
for ClustFile in $MIDFolderName/ClustCons/multi-seq/*.fasta
	do
	### Get the clustername without extention
	Clust=$(echo "$ClustFile" | sed -e 's/\/.*\///g' | cut -d. -f1)
	ReadCount=$(cat $ClustFile | grep "^>" | wc -l)

	# Check the if the number of reads is too large
	if (($ReadCount > $RandSize)); then
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}END {printf("\n");}' < $ClustFile | awk 'NR>1{ printf("%s",$0); n++; if(n%2==0) { printf("\n");} else { printf("\t");} }' | awk -v k=$RandSize 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | awk -F "\t" '{print $1"\n"$2 > "'$MIDFolderName/ClustCons/multi-seq/$Clust.Rand$RandSize.fasta'"}'

	### Perform a MSA with setting: 'L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information)'
	mafft --adjustdirection --thread $THREADS --localpair --maxiterate 1000 $MIDFolderName/ClustCons/multi-seq/$Clust.Rand$RandSize.fasta > $MIDFolderName/ClustCons/multi-seq/$Clust.mafft.msf
	### Create consensus
	perl $scripts/consensus.pl -in $MIDFolderName/ClustCons/multi-seq/$Clust.mafft.msf -t 30 -out $MIDFolderName/ClustCons/multi-seq/$Clust.cons.msf
	###Remove question marks and reformat FASTA
	cat $MIDFolderName/ClustCons/multi-seq/$Clust.cons.msf | sed 's/?//g' | awk '{if (substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' > $MIDFolderName/ClustCons/multi-seq/$Clust.cons.fasta
	### Rename the file to select the correct one
	awk '/^>/{print ">'$Clust';size='$ReadCount'"; next}{print}' < $MIDFolderName/ClustCons/multi-seq/$Clust.cons.fasta > $MIDFolderName/ClustCons/multi-seq/$Clust.renamed.cons.fasta
		else
	### Perform a MSA with setting: 'L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information)'
	>&2 echo "MAFFT"
	mafft --adjustdirection --thread $THREADS --localpair --maxiterate 1000 $ClustFile > $MIDFolderName/ClustCons/multi-seq/$Clust.mafft.msf
	### Create consensus
	>&2 echo "Consensus"
	perl $scripts/consensus.pl -in $MIDFolderName/ClustCons/multi-seq/$Clust.mafft.msf -t 30 -out $MIDFolderName/ClustCons/multi-seq/$Clust.cons.msf
	###Remove question marks and reformat FASTA
	>&2 echo "Rename"
	cat $MIDFolderName/ClustCons/multi-seq/$Clust.cons.msf | sed 's/?//g' | awk '{if (substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' > $MIDFolderName/ClustCons/multi-seq/$Clust.cons.fasta
	### Rename the file to select the correct one
	awk '/^>/{print ">'$Clust';size='$ReadCount'"; next}{print}' < $MIDFolderName/ClustCons/multi-seq/$Clust.cons.fasta > $MIDFolderName/ClustCons/multi-seq/$Clust.renamed.cons.fasta

	fi
done

### Combine all the fasta files from the clusters, and remove all empty newlines
>&2 echo "Combine"
cat $MIDFolderName/ClustCons/multi-seq/*.renamed.cons.fasta | awk NF > $MIDFolderName/ClustCons/multi-seq/$MID.$SampleName.combined.fasta
}



####Run workflow
if [[ $SampleName != '' ]] ; then
	UserInput
	if [[ "$Runmodules" == *"all"* ]] || [[ "$Runmodules" == *"QCtrimming"* ]] ; then
		time QCtrimming
	fi
	if [[ "$Runmodules" == *"all"* ]] || [[ "$Runmodules" == *"lustering"* ]] ; then
		time Clustering
	fi
	if [[ "$Runmodules" == *"all"* ]] || [[ "$Runmodules" == *"onsensus"* ]] ; then
		time Consensus
	fi

fi
