
### Load modules
#module load WUR/RIKILT/cutadapt-1.16 SHARED/metabarcoding/.ncbi-blast-2.10.1 R/3.5.0
#module load WUR/RIKILT/cutadapt-1.16 R/4.1.2 WUR/RIKILT/openpyxl/1.0 WUR/RIKILT/blast-plus/2.12.0
### Input variables
RunName=${1}
OutDir=${2}
OutFolderName=${3}
Ident1=${4}
PrimerFile=${5}
Error=${6}
Evalue=${7}
Qcov=${8}
Pident=${9}
MaxTargetSeqs=${10}
SampleDescription=${11}
THREADS=${12}
Runmodules=${13}
NT_dmp=${14}
MID=${15}
ExTaxids=${16}
ExSeqids=${17}
##Additional variables
SampleName=$(cat $SampleDescription | grep "^$MID" | awk -F";" '{print $2}')
MIDFolderName="$OutFolderName/$MID.$SampleName"
RandSize=100
nt=true #Use nt database instead of custom barcode databases

### Database variables
export BLASTDB=$BLASTDB:/lustre/shared/wfsr-databases/BLASTdb
DATABASE=/lustre/shared/wfsr-databases/BLASTdb/nt


####Select reads containing the DNA forward and reverse barcode primers. Primers are also removed.
function SelectBarcodeBlast(){
	# Create the output directory
	rm -r --force $MIDFolderName/Barcodes/
	rm -r --force $MIDFolderName/BLAST/
	mkdir -p $MIDFolderName/Barcodes
	### Get all the targets from the SampleDescription-file
	for TARGET in $(cat $SampleDescription | grep "$MID" | awk -F\; '{print $3}' | tr , '\n') ; do
	###ADD what to do with EMPTY files
	printf "Runnning SelectBarcodeBlast"
	printf "DNA barcode target: $TARGET\n"

		# Create two strings for the cutadapt statement.
		Primers=$(cat $PrimerFile | grep "^$TARGET;" | awk -F\; '{print $2}' | tr "[Ii]" "[Nn]" | sed 's/,*$//' | sed 's/,/ -b  /g' | sed -e 's/^/ -b /' | tr -d '\n')
		PrimersRC=$(cat $PrimerFile | grep "^$TARGET;" | awk -F\; '{print $2}' | tr "[ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVvIi]" "[TACGAtacgaNnRrYySsWwMmKkVvHhDdBbNn]" | rev | sed 's/,*$//' | sed 's/,/ -b  /g' | sed -e 's/^/ -b /' | tr -d '\n')


		### Combine the primers
		Combined=$Primers" "$PrimersRC
		echo $Combined

		### Cutadapt: remove up to 2 adapters from each read (-n 2) anywhere in the sequence (-b)
		Overlap=25
		cutadapt $Combined --overlap $Overlap -n 2 -o $MIDFolderName/Barcodes/$MID.$SampleName.$TARGET.PS.fasta -e $Error $MIDFolderName/ClustCons/multi-seq/$MID.$SampleName.combined.fasta -m 10 -j $THREADS --untrimmed-output $MIDFolderName/Barcodes/$MID.$SampleName.$TARGET.untrimmed.fasta > $MIDFolderName/Barcodes/$MID.$SampleName.$TARGET.PS.out

		### BLAST
		mkdir -p $MIDFolderName/BLAST

		if [ -s "$MIDFolderName/Barcodes/$MID.$SampleName.$TARGET.PS.fasta" ]
			then
				blastn -query $MIDFolderName/Barcodes/$MID.$SampleName.$TARGET.PS.fasta -task megablast -db $DATABASE -out $MIDFolderName/BLAST/$MID.$SampleName.$TARGET.blast.tsv -outfmt '6 qseqid qlen qcovhsp pident bitscore evalue sallacc staxids stitle sscinames' -num_threads $THREADS -evalue $Evalue -max_target_seqs $MaxTargetSeqs  -negative_taxidlist $ExTaxids
				printf "qseqid\tqlen\tqcovhsp\tpident\tbitscore\tevalue\tsallacc\tstaxids\tstitle\tsscinames\n" | cat - $MIDFolderName/BLAST/$MID.$SampleName.$TARGET.blast.tsv >  $MIDFolderName/BLAST/temp && mv  $MIDFolderName/BLAST/temp $MIDFolderName/BLAST/$MID.$SampleName.$TARGET.blast.tsv
				### Filter BLAST output
				cat $MIDFolderName/BLAST/$MID.$SampleName.$TARGET.blast.tsv | awk -F "\t" '{ if(($3 >= '"$Qcov"') && ($4 >= '"$Pident"')) { print } }' > $MIDFolderName/BLAST/$MID.$SampleName.$TARGET.blast.filtered.tsv
			else
				printf "\n$MIDFolderName/Barcodes/$MID.$SampleName.$TARGET.PS.fasta does not exist, or is empty\n"
		fi
	done
}




function SummaryTax_reads(){
#        source /cm/shared/apps/WUR/RIKILT/openpyxl/bin/activate
        rm -r --force $MIDFolderName/Taxonomy
        mkdir -p $MIDFolderName/Taxonomy
	python scripts/Taxonomy.py -i  $MIDFolderName/BLAST/ -o $MIDFolderName/Taxonomy/ -p $NT_dmp/ -e $ExSeqids
	echo "python scripts/Taxonomy.py -i  $MIDFolderName/BLAST/ -o $MIDFolderName/Taxonomy/ -p $NT_dmp/ -e $ExSeqids"
#	deactivate
        #Creating one summary for all barcodes
        Rscript Taxonomic_summary.R $MIDFolderName


	printf "qseqid\tTAXID\tTaxonomic Level\tTaxonomic Name\tNumber of Top Hits\tLineage\tReads\tMarker\n" > $MIDFolderName/Taxonomy/$MID.$SampleName.Overview.summary
	for i in $MIDFolderName/Taxonomy/*.Taxonomy.summary; do
	sed '1d' $i >> $MIDFolderName/Taxonomy/$MID.$SampleName.Overview.summary
	python scripts/short_summary.py $i > $i.short
	echo "Common name" > $i.common
	for taxid in $(sed '1d' $i.short | awk -F '\t' '{print $3}'); do
		grep "genbank common name" $NT_dmp/names.dmp | grep -P "^${taxid}\t" | awk -F '\t|\t' '{print $3}' | tail -n1 >> $i.common
		if [[ $(grep "genbank common name" $NT_dmp/names.dmp | grep -P "^${taxid}\t" | awk -F '\t|\t' '{print $3}' | tail -n1) == '' ]]; then
			printf "None\n" >> $i.common
		fi
	done
	paste -d '\t' $i.short $i.common > $i.tmp; mv $i.tmp $i.short

	#rename $i to summary
	cp $i $OutFolderName/Taxonomy/
	cp $i.short $OutFolderName/Taxonomy/
	done
	rm  $MIDFolderName/Taxonomy/*.common
	cp $MIDFolderName/$MID.$SampleName.Taxonomic.summary.all_barcodes.tsv $OutFolderName/Taxonomy/

}

### Get a summary of each barcode within a sample, these summaries are collected in a directory with summaries from all samples within a dataset. This directory is used for mail-pie-plotter.r
function SummaryTax_reads_old(){
	printf "qseqid\tTAXID\tTaxonomic Level\tTaxonomic Name\tNumber of Top Hits\tLineage\tReads\tLength\tMarker\n" > $MIDFolderName/Taxonomy/$MID.$SampleName.Overview.summary

	for i in $MIDFolderName/Taxonomy/*.Taxonomy.id.tsv
	do
    	echo $i
	echo "Reads"  > $i.count_reads
	echo "Length" > $i.read_length
	echo "Marker" > $i.marker
	marker=$(basename $i|awk -F '.' '{print $3}')
	#retrieving reads total and consensus length per cluster/qseqid
	for cluster in $(cat $i | awk -F '\t' '{ print $1 }' | grep -v 'qseqid' | awk -F ';' '{print $1}')
		do
			#echo $cluster
			grep '^>' $MIDFolderName/ClustCons/multi-seq/$cluster.fasta | wc -l >> $i.count_reads
			cat $MIDFolderName/ClustCons/multi-seq/$cluster.renamed.cons.fasta | grep -v "^>" | tr -d "\n" | wc -m >> $i.read_length
			printf "$marker\n" >> $i.marker
		done

	#adding reads and length per cluster
	paste -d '\t' $i $i.count_reads > $i.tmp ; mv $i.tmp $i
	paste -d '\t' $i $i.read_length > $i.tmp ; mv $i.tmp $i

	sed '1d' $i | sort -t$'\t' -k4 > $i.tmp

	#Adding results to overview (used for xlsx file)
        paste -d '\t' $i $i.marker > $i.tmp2 ; mv $i.tmp2 $i
	sed '1d' $i >> $MIDFolderName/Taxonomy/$MID.$SampleName.Overview.summary

	#create short summary , giving total reads and clusters per organism
	python short_summary.py $i.tmp > $i.short
	echo "Common name" > $i.common
	for taxid in $(sed '1d' $i.short | awk -F '\t' '{print $3}'); do
		grep "genbank common name" $NT_dmp/names.dmp | grep -P "^${taxid}\t" | awk -F '\t|\t' '{print $3}' | tail -n1 >> $i.common
		if [[ $(grep "genbank common name" $NT_dmp/names.dmp | grep -P "^${taxid}\t" | awk -F '\t|\t' '{print $3}' | tail -n1) == '' ]]; then
			printf "None\n" >> $i.common
		fi
	done
	paste -d '\t' $i.short $i.common > $i.tmp; mv $i.tmp $i.short

	#rename $i to summary
	mv $i $i.summary
	cp $i.summary $OutFolderName/Taxonomy/
	cp $i.short $OutFolderName/Taxonomy/
	done


	rm *.count_reads *.common *.read_length *.marker
	#Creating one summary for all barcodes
	Rscript scripts/Taxonomic_summary.R $MIDFolderName
	cp $MIDFolderName/$MID.$SampleName.Taxonomic.summary.all_barcodes.tsv $OutFolderName/Taxonomy/
}


function Create_xslx(){
	source /cm/shared/apps/WUR/RIKILT/openpyxl/bin/activate
       #Creating xlsx file
        python scripts/xlsx.py --out $MIDFolderName/Taxonomy/$MID.$SampleName.Overview.xlsx --input $MIDFolderName/Taxonomy/$MID.$SampleName.Overview.summary
        for i in $MIDFolderName/Taxonomy/*.Taxonomy.summary; do
                python scripts/xlsx.py --out $MIDFolderName/Taxonomy/$MID.$SampleName.Overview.xlsx --input $i
        done
        cp  $MIDFolderName/Taxonomy/$MID.$SampleName.Overview.xlsx $OutFolderName/Taxonomy/
	deactivate
}

SelectBarcodeBlast 
printf "BLAST done\n"
SummaryTax_reads
printf "Tax done\n"
Create_xslx 
