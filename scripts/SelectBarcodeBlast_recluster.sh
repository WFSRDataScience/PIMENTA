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
THREADS=${11}
Runmodules=${12}
NT_dmp=${13}
Targets=${14}
ExTaxids=${15}
ExSeqids=${16}
DATABASE=${17}
scripts=${18}
##Additional variables
RandSize=100
nt=true #Use nt database instead of custom barcode databases

### Database variables
folderDB=$(dirname $DATABASE)
export BLASTDB=$BLASTDB:$folderDB


####Select reads containing the DNA forward and reverse barcode primers. Primers are also removed.
function SelectBarcodeBlast(){
	# Create the output directory
	rm -r --force $OutFolderName/Barcodes/
	rm -r --force $OutFolderName/BLAST/
	mkdir -p $OutFolderName/Barcodes
	### Get all the targets from the SampleDescription-file
	for TARGET in $(echo $Targets | tr , '\n') ; do
	###ADD what to do with EMPTY files
	printf "Runnning SelectBarcodeBlast"
	printf "DNA barcode target: $TARGET\n"

		# Create two strings for the cutadapt statement.
		Primers=$(cat $PrimerFile | grep "^$TARGET;" | grep "forward" | awk -F\; '{print $2}' | tr "[Ii]" "[Nn]" | sed 's/,*$//' | sed 's/,/ -b  /g' | sed -e 's/^/ -g /' | tr -d '\n')
		PrimersRC=$(cat $PrimerFile | grep "^$TARGET;" | grep "reverse" | awk -F\; '{print $2}' | tr "[ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVvIi]" "[TACGAtacgaNnRrYySsWwMmKkVvHhDdBbNn]" | rev | sed 's/,*$//' | tr -d '\n')
		#PrimersRC=$(cat $PrimerFile | grep "^$TARGET;" | grep "reverse" | awk -F\; '{print $2}' | tr "[Ii]" "[Nn]" | sed 's/,*$//' | tr -d '\n')

		### Combine the primers
		Combined="$Primers -b $PrimersRC"
		echo $Combined
	
		### Cutadapt: remove up to 2 adapters from each read (-n 2) anywhere in the sequence (-b), -m minimum length
		Overlap=10
		cutadapt $Combined --overlap $Overlap -n 2 -o $OutFolderName/Barcodes/$RunName.$TARGET.PS.fasta -e $Error $OutFolderName/ClustCons/multi-seq/$RunName.combined.fasta -m 10 -j $THREADS --untrimmed-output $OutFolderName/Barcodes/$RunName.$TARGET.untrimmed.fasta > $OutFolderName/Barcodes/$RunName.$TARGET.PS.out --rc
		sed -E '/^>/s/( +|\t).*//' $OutFolderName/Barcodes/$RunName.$TARGET.PS.fasta > tmp; mv tmp $OutFolderName/Barcodes/$RunName.$TARGET.PS.fasta
		### BLAST
		mkdir -p $OutFolderName/BLAST
		
		if [ -s "$OutFolderName/Barcodes/$RunName.$TARGET.PS.fasta" ]
			then
				cat $OutFolderName/Barcodes/$RunName.$TARGET.PS.fasta | perl -ne 's/[^\x00-\x7F]+/ /g; print;' > tmp.fasta; mv tmp.fasta $OutFolderName/Barcodes/$RunName.$TARGET.PS.fasta
				time blastn -query $OutFolderName/Barcodes/$RunName.$TARGET.PS.fasta -task megablast -db $DATABASE -out $OutFolderName/BLAST/$RunName.$TARGET.blast.tsv -outfmt '6 qseqid qlen qcovhsp pident bitscore evalue sallacc staxids stitle sscinames' -num_threads $THREADS -evalue $Evalue -max_target_seqs $MaxTargetSeqs  -negative_taxidlist $ExTaxids
				printf "qseqid\tqlen\tqcovhsp\tpident\tbitscore\tevalue\tsallacc\tstaxids\tstitle\tsscinames\n" | cat - $OutFolderName/BLAST/$RunName.$TARGET.blast.tsv >  $OutFolderName/BLAST/temp && mv  $OutFolderName/BLAST/temp $OutFolderName/BLAST/$RunName.$TARGET.blast.tsv
				### Filter BLAST output
				cat $OutFolderName/BLAST/$RunName.$TARGET.blast.tsv | awk -F "\t" '{ if(($3 >= '"$Qcov"') && ($4 >= '"$Pident"')) { print } }' > $OutFolderName/BLAST/$RunName.$TARGET.blast.filtered.tsv
			else
				printf "\n$OutFolderName/Barcodes/$RunName.$TARGET.PS.fasta does not exist, or is empty\n"
		fi
	done
}


function SummaryTax_reads(){
	rm -r --force $OutFolderName/Taxonomy_$Ident1
	mkdir -p $OutFolderName/Taxonomy_$Ident1

	python $scripts/Taxonomy.py -i  $OutFolderName/BLAST/ -o $OutFolderName/Taxonomy_$Ident1/ -p $NT_dmp/ -e $ExSeqids
	printf "qseqid\tTAXID\tTaxonomic Level\tTaxonomic Name\tNumber of Top Hits\tLineage\tReads\tMarker\n" > $OutFolderName/Taxonomy_$Ident1/$RunName.Overview.summary
	for i in $OutFolderName/Taxonomy_$Ident1/*.Taxonomy.summary; do
	sed '1d' $i >> $OutFolderName/Taxonomy_$Ident1/$RunName.Overview.summary
	python $scripts/short_summary.py $i > $i.short
	echo "Common name" > $i.common
	for taxid in $(sed '1d' $i.short | awk -F '\t' '{print $3}'); do
		grep "genbank common name" $NT_dmp/names.dmp | grep -P "^${taxid}\t" | awk -F '\t|\t' '{print $3}' | tail -n1 >> $i.common
		if [[ $(grep "genbank common name" $NT_dmp/names.dmp | grep -P "^${taxid}\t" | awk -F '\t|\t' '{print $3}' | tail -n1) == '' ]]; then
			printf "None\n" >> $i.common
		fi
	done
	paste -d '\t' $i.short $i.common > $i.tmp; mv $i.tmp $i.short

	done

	rm $OutFolderName/Taxonomy_$Ident1/*.common
}


function Create_xslx(){
       #Creating xlsx file
        python $scripts/xlsx.py --out $OutFolderName/Taxonomy/$RunName.Overview.xlsx --input $OutFolderName/Taxonomy_$Ident1/$RunName.Overview.summary
        for i in $OutFolderName/Taxonomy_$Ident1/*.Taxonomy.summary; do
                python $scripts/xlsx.py --out $OutFolderName/Taxonomy_$Ident1/$RunName.Overview.xlsx --input $i
        done
}

if [[ "$Runmodules" == *"all"* ]] || [[ "$Runmodules" == *"blast"* ]] || [[ "$Runmodules" == *"BLAST"* ]] || [[ "$Runmodules" == *"Blast"* ]] ; then
	SelectBarcodeBlast
fi
if [[ "$Runmodules" == *"all"* ]] || [[ "$Runmodules" == *"axonomy"* ]] ; then
	SummaryTax_reads
#	Create_xslx
fi 
