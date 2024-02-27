#!/bin/bash

##Modified 15 November 2022, V. van der Vorst


### Input variables
RunName=${1}
OutDir=${2}
OutFolderName=${3}
Ident1=${4}
MinClustSize=${5}
PrimerFile=${6}
Error=${7}
Evalue=${8}
Qcov=${9}
Pident=${10}
MaxTargetSeqs=${11}
THREADS=${12}
Runmodules=${13}
NT_dmp=${14}
ExTaxids=${15}
Targets=${16}
ExSeqids=${17}
modus=${18}
DATABASE=${19}
filtertax=${20}
scripts=${21}

if [[ "$modus" == "one" ]]; then
	files=$(find $OutFolderName/barcode*.*/ClustCons/multi-seq/*.combined.fasta)
else
	files=$(find $OutDir/*/barcode*.*/ClustCons/multi-seq/*.combined.fasta)
fi
printf "accessionTaxa $accessionTaxa"
printf "DB $DATABASE\n"
printf "scripts $scripts\n"

Prepare_cons_fasta(){
        for i in $files; do name=$(basename $i | awk -F '.PS.fasta' '{print $1}') ;perl -p -e 's/^(>.*)$/$1-'$name'/g' $i ; done > $OutFolderName/$RunName.PS.fasta
}




### CD-HIT clustering
Clustering(){
	rm -r --force $OutFolderName/ClustCons/multi-seq/
	rm -r --force $OutFolderName/ClustCons/
	mkdir -p $OutFolderName/ClustCons/
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
	cd-hit-est -i $OutFolderName/$RunName.PS.fasta -o $OutFolderName/ClustCons/$RunName.CD-hit.90 -c 0.90 -G 1 -g 1 -gap -1 -gap-ext 0 -p 1 -r 1 -d 0 -T $THREADS -n 8 -M 20000
	printf "\nRunning CD-HIT-est...RUN2 at $Int%% identity\n"
	cd-hit-est -i $OutFolderName/ClustCons/$RunName.CD-hit.90 -o $OutFolderName/ClustCons/$RunName.CD-hit.$Ident1 -c $Ident1 -G 1 -g 1 -gap -1 -gap-ext 0 -p 1 -r 1 -d 0 -T $THREADS -n $n -M 20000
	###Combine the two clstr-output files to include all original sequences from the parent clstr-file from RUN1
	clstr_rev.pl $OutFolderName/ClustCons/$RunName.CD-hit.90.clstr $OutFolderName/ClustCons/$RunName.CD-hit.$Ident1.clstr > $OutFolderName/ClustCons/$RunName.CD-hit.90.$Ident1.clstr
	###Renumber clusters
	clstr_renumber.pl $OutFolderName/ClustCons/$RunName.CD-hit.90.$Ident1.clstr > $OutFolderName/ClustCons/$RunName.CD-hit.90.$Ident1.renumber.clstr
	###Calculate cluster sizes of merged clstr-files.
	plot_len1.pl $OutFolderName/ClustCons/$RunName.CD-hit.90.$Ident1.renumber.clstr 1,2-5,6-9,10-19,20-29,30-39,40-49,50-100,101-200,201-300,301-400,401-500,501-1000,1001-5000,5001-99999 > $OutFolderName/ClustCons/$RunName.CD-hit.90.$Ident1.renumber.plotlen
	###Extract FASTA sequences from clusters
	printf "\nGenerating fasta files for clusters with more than $MinClustSize member sequences\n"
	make_multi_seq.pl $OutFolderName/$RunName.PS.fasta $OutFolderName/ClustCons/$RunName.CD-hit.90.$Ident1.renumber.clstr $OutFolderName/multi-seq $MinClustSize	
	mv --force $OutFolderName/multi-seq/ $OutFolderName/ClustCons/ #############SIMPLIFY
	find $OutFolderName/ClustCons/ -type f -not -name "*.*" -exec mv "{}" "{}".fasta \;
} 



###Consensus building
Consensus(){
printf "\nRunning MAFFT MSA and consensus building....\n"
### Get a consensus sequence per cluster
for ClustFile in $OutFolderName/ClustCons/multi-seq/*.fasta
	do
	### Get the clustername without extention
	Clust=$(echo "$ClustFile" | sed -e 's/\/.*\///g' | cut -d. -f1)
	ReadCount=$(cat $ClustFile | grep "^>" | wc -l)
	RandSize=100
	if (($ReadCount == 1)); then
		printf "$Clust;size=$ReadCount"
		awk '/^>/{print ">'$Clust'"; next}{print}' < $OutFolderName/ClustCons/multi-seq/$Clust.fasta > $OutFolderName/ClustCons/multi-seq/$Clust.renamed.cons.fasta	
        elif (($ReadCount > $RandSize)); then
        awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}END {printf("\n");}' < $ClustFile | awk 'NR>1{ printf("%s",$0); n++; if(n%2==0) { printf("\n");} else { printf("\t");} }' | awk -v k=$RandSize 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | awk -F "\t" '{print $1"\n"$2 > "'$OutFolderName/ClustCons/multi-seq/$Clust.Rand$RandSize.fasta'"}'

        ### Perform a MSA with setting: 'L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information)'
        mafft --adjustdirection --thread $THREADS --localpair --maxiterate 1000 $OutFolderName/ClustCons/multi-seq/$Clust.Rand$RandSize.fasta > $OutFolderName/ClustCons/multi-seq/$Clust.mafft.msf
        ### Create consensus
        perl $scripts/consensus.pl -in $OutFolderName/ClustCons/multi-seq/$Clust.mafft.msf -t 30 -out $OutFolderName/ClustCons/multi-seq/$Clust.cons.msf
        ###Remove question marks and reformat FASTA
        cat $OutFolderName/ClustCons/multi-seq/$Clust.cons.msf | sed 's/?//g' | awk '{if (substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' > $OutFolderName/ClustCons/multi-seq/$Clust.cons.fasta
        ### Rename the file to select the correct one
        awk '/^>/{print ">'$Clust'"; next}{print}' < $OutFolderName/ClustCons/multi-seq/$Clust.cons.fasta > $OutFolderName/ClustCons/multi-seq/$Clust.renamed.cons.fasta
	else
	### Perform a MSA with setting: 'L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information)'
	mafft --adjustdirection --thread $THREADS --localpair --maxiterate 1000 $ClustFile > $OutFolderName/ClustCons/multi-seq/$Clust.mafft.msf
	### Create consensus
	perl $scripts/consensus.pl -in $OutFolderName/ClustCons/multi-seq/$Clust.mafft.msf -t 30 -out $OutFolderName/ClustCons/multi-seq/$Clust.cons.msf &>/dev/null
	###Remove question marks and reformat FASTA
	cat $OutFolderName/ClustCons/multi-seq/$Clust.cons.msf | sed 's/?//g' | awk '{if (substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' > $OutFolderName/ClustCons/multi-seq/$Clust.cons.fasta
	### Rename the file to select the correct one
	awk '/^>/{print ">'$Clust'"; next}{print}' < $OutFolderName/ClustCons/multi-seq/$Clust.cons.fasta > $OutFolderName/ClustCons/multi-seq/$Clust.renamed.cons.fasta

	fi
done

### Combine all the fasta files from the clusters, and remove all empty newlines
cat $OutFolderName/ClustCons/multi-seq/*.renamed.cons.fasta | awk NF > $OutFolderName/ClustCons/multi-seq/$RunName.combined.fasta
}


###Count the total number of reads across clusters of each barcode
function CountClusterSize(){
        header="DNA_barcode\tClusterID\tSize\tTAXID\tTaxLevel\tTaxonomy\tLineage\tClusters_from_samples"
        printf "$header\n" > $OutFolderName/$RunName.ClusterContent.tsv
        for TARGET in $(echo $Targets | tr , '\n') ; do
                for CLUSTER in $(cat $OutFolderName/Barcodes/$RunName.$TARGET.PS.fasta | grep "^>" | sed 's/>//g') ; do
                        Taxonomy=$(grep -P "^$CLUSTER\t" $OutFolderName/Taxonomy_$Ident1/$RunName.Overview.summary | awk -F '\t' '{print $2"\t"$3"\t"$4"\t"$6}')
                        if [[ "$Taxonomy" == "" ]]; then
                                Taxonomy="NA\tNA\tNA\tNA"
                        fi
                        SIZE=$(cat $OutFolderName/ClustCons/multi-seq/$CLUSTER.fasta | grep "^>" | wc -l)

                        counts="$TARGET\t$CLUSTER\t$SIZE\t$Taxonomy"

			content=$(grep "^>" $OutFolderName/ClustCons/multi-seq/$CLUSTER.fasta | sed 's/>//g'| tr '\n' ':')
                        printf "$counts\t$content\n" >> $OutFolderName/$RunName.ClusterContent.tsv

                done
        done
}

function GetTaxPerSample(){
	rm -r $OutFolderName/Taxonomy_per_sample_$Ident1
	mkdir -p $OutFolderName/Taxonomy_per_sample_$Ident1

	#Go through each marker of each sample to retrieve sample specific taxonomy
	for Sample in $files; do
                clusterid=0
		SampleName=$(basename $Sample | awk -F '.' '{print $2}' | tr -d '\n')
		if  [ "$SampleName" == "combined" ]; then 
			SampleName=$(basename $Sample | awk -F '.' '{print $1}' | tr -d '\n')
		fi
		echo "Sample: $SampleName"
		#retrieve all within a marker 
		grep "$SampleName" $OutFolderName/$RunName.ClusterContent.tsv | tr ' ' ':' | grep -vP 'NA\tNA\tNA\tNA'| awk -F '\t' '{print $1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}'  > $OutFolderName/Taxonomy_per_sample_$Ident1/$SampleName.tmp
   		printf "#DNA_barcode\tTaxid\tLevel\tTaxonomy\tLineage\tReadcount\n"  > $OutFolderName/Taxonomy_per_sample_$Ident1/$SampleName.count.tmp
		cat $OutFolderName/Taxonomy_per_sample_$Ident1/$SampleName.tmp
		while read -r marker taxid level tax lineage barcodes; do
			let clusterid++
			readcount=0
			#retrieve sample specific count for each cluster
			for barcode in $(echo $barcodes | tr ';' '\n'); do
				if [[ "$barcode" == *"$SampleName"* ]]; then
					bcount=$(echo $barcode | awk -F '=' '{print $2}' | awk -F '-' '{print $1}')
					readcount=$((readcount+bcount))
				fi
			done	
			if [ "$filtertax" != "notaxfilter" ]; then
			
			if [ "$lineage" == *"$filtertax"* ] ; then
    	                	printf "${marker}.${clusterid}\t$taxid\t$level\t$tax\t$lineage\t$readcount\n" | tr ':' ' ' >> $OutFolderName/Taxonomy_per_sample_$Ident1/$SampleName.count.tmp
			fi
			else
				printf "${marker}.${clusterid}\t$taxid\t$level\t$tax\t$lineage\t$readcount\n" | tr ':' ' ' >> $OutFolderName/Taxonomy_per_sample_$Ident1/$SampleName.count.tmp
			fi
		done < $OutFolderName/Taxonomy_per_sample_$Ident1/$SampleName.tmp
	for TARGET in $(echo $Targets | tr , '\n') ; do
        	mkdir -p $OutFolderName/Taxonomy_per_sample_$Ident1/$TARGET
		targetoutput=$(grep "$TARGET" $OutFolderName/Taxonomy_per_sample_$Ident1/$SampleName.count.tmp)
		echo $targetoutput
		if [ "$targetoutput" != "" ]; then
      		#	printf "#DNA_barcode\tTaxid\tlevel\tTaxonomy\tLineage\tReadcount\n"  > $OutFolderName/Taxonomy_per_sample_$Ident1/$TARGET/$SampleName.$TARGET.tsv
		  	echo "$targetoutput" > $OutFolderName/Taxonomy_per_sample_$Ident1/$TARGET/$SampleName.$TARGET.tsv
    		fi
	done
	done
	rm $OutFolderName/Taxonomy_per_sample_$Ident1/*.tmp


}

function GetKronaPlots(){
	mkdir $OutFolderName/Taxonomy_per_sample_$Ident1/krona
        for TARGET in $(echo $Targets | tr , '\n') ; do
		mkdir -p $OutFolderName/Taxonomy_per_sample_$Ident1/krona/$TARGET	
		#printf "DNA_barcode\tlevel\tTaxonomy\tLineage\tReadcount\n"  > $OutFolderName/Taxonomy_per_sample_$Ident1/krona/$TARGET/$TARGET.tsv
		for i in  $OutFolderName/Taxonomy_per_sample_$Ident1/$TARGET/*.$TARGET.tsv ; do
			sed 1d $i >> $OutFolderName/Taxonomy_per_sample_$Ident1/krona/$TARGET/$TARGET.tsv
		done
		ktImportTaxonomy -q 1 -t 2 -m 6 -o  $OutFolderName/Taxonomy_per_sample_$Ident1/krona/$TARGET.Taxonomy.krona.html $OutFolderName/Taxonomy_per_sample_$Ident1/krona/$TARGET/$TARGET.tsv
		for ksample in $OutFolderName/Taxonomy_per_sample_$Ident1/$TARGET/*.$TARGET.tsv; do
			SName=$(basename $ksample | awk -F '.' '{print $1}' | tr -d '\n')
			echo $SName
			ktImportTaxonomy -q 1 -t 2 -m 6 -o  $OutFolderName/Taxonomy_per_sample_$Ident1/krona/$TARGET/$SName.$TARGET.Taxonomy.krona.html $ksample
		done
	done
}	

function GetStats(){
	printf "Statistics for dataset:$RunName\n" > $OutFolderName/$RunName.stats.txt
	printf "Total PS $RunName consensus sequences: $(cat $OutFolderName/$RunName.PS.fasta | grep "^>" | wc -l)\n" >> $OutFolderName/$RunName.stats.txt
	printf "Total Clusters: $(cat $OutFolderName/ClustCons/$RunName.CD-hit.90.$Ident1.renumber.plotlen | grep "Total" | awk '{print $3}')\n" >> $OutFolderName/$RunName.stats.txt
	printf "Total Clusters >= $MinClustSize: $(cat $OutFolderName/ClustCons/multi-seq/$RunName.combined.fasta | grep "^>" | wc -l)\n" >> $OutFolderName/$RunName.stats.txt
	for TARGET in $(echo $Targets | tr , '\n') ; do
		printf "Total $TARGET Clusters: $(cat $OutFolderName/Barcodes/$RunName.$TARGET.PS.fasta | grep "^>" | wc -l)\n" >> $OutFolderName/$RunName.stats.txt
                Rscript $scripts/Taxonomic_summary_recluster.R $OutFolderName/Taxonomy_per_sample_$Ident1/$TARGET
		Rscript $scripts/Prepare_phylosec.R  $OutFolderName $OutFolderName/$RunName.ClusterContent.tsv $OutFolderName/Barcodes/$RunName.$TARGET.PS.fasta $TARGET
	done
}


####Run workflow
	if [[ "$Runmodules" == *"all"* ]] || [[ "$Runmodules" == *"lustering"* ]] ; then
		Prepare_cons_fasta
		start=$SECONDS
                time Clustering
		duration=$(( SECONDS - start ))
                echo "Clustering $duration"
	fi
	if [[ "$Runmodules" == *"all"* ]] || [[ "$Runmodules" == *"onsensus"* ]] ; then
                start=$SECONDS
                time Consensus
                duration=$(( SECONDS - start ))
                echo "Consensus $duration"
	fi
	if [[ "$Runmodules" == *"all"* ]] || [[ "$Runmodules" == *"blast"* ]] || [[ "$Runmodules" == *"BLAST"* ]] || [[ "$Runmodules" == *"Blast"* ]] || [[ "$Runmodules" == *"axonomy"* ]]; then
		bash $scripts/SelectBarcodeBlast_recluster.sh $RunName $OutDir $OutFolderName $Ident1 $PrimerFile $Error $Evalue $Qcov $Pident $MaxTargetSeqs $THREADS $Runmodules $NT_dmp $Targets $ExTaxids $ExSeqids $DATABASE $scripts
	fi
	CountClusterSize
        GetTaxPerSample
	GetStats
	GetKronaPlots
	
