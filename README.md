# Pimenta DNA metabarcoding combined clustering
This pipeline combines cluster data from multiple datasets/samples analysed with the DNA Metabarcoding pipeline (https://git.wur.nl/vorst016/hpc-dna-metabarcoding-minion-gpu-parallel) and reclusters the consensus sequences from these datasets/samples. 
New consensus sequences are created from these new clusters and the taxonomy is identified with a BLAST analysis. 

![alt text](https://git.wur.nl/vorst016/pimenta/-/raw/master/DNA_metabarcoding.drawio.png)

<strong>DNA metabarcoding</strong> <br>
Pimenta is a pipeline for the rapid identification of species in forensic samples using MinION data. The pipeline consists of eight linked tools, and data analysis passes through 3 phases: 1) pre-processing the MinION data through read calling, demultiplexing, trimming sequencing adapters, quality trimming and filtering the reads, 2) clustering the reads, continued by MSA and consensus building per cluster, 3) Taxonomy identification with the use of a BLAST analysis. The DNAmetabarcoding pipeline makes use of the frequently-used software tools Cutadapt v1.16 (http://cutadapt.readthedocs.io/en/stable/guide.html), PRINSEQ v0.20.4 (http://prinseq.sourceforge.net/), Porechop (https://github.com/rrwick/Porechop), BLAST+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), Guppy 4.0.11, CD-HIT v4.6.7 (https://github.com/weizhongli/cdhit), Consensus v1.0, Krona v2.7.1 (https://github.com/marbl/Krona/wiki) and MAFFT v7.130 (https://mafft.cbrc.jp/alignment/software/).
 

<strong>Required programs</strong> <br>
Mininconda has to be installed and loaded before installing Pimenta. <br>

<strong>Installing the pipeline</strong> <br>
First clone the pipeline into a chosen directory:
```
git clone https://github.com/WFSRDataScience/Pimenta.git
```

To install Pimenta:<br>
```
cd pimenta
bash install.sh 
```


<strong>Required R packages</strong> <br>dplyr, readr, stringi, taxonomizr, stringr and data.table are required for Taxonomic_summary.R <br>
The required R packages are automatically installed.

<strong> Required perl packages</strong> <br> Bio::SeqIO, Bio::AlignIO <br>
These are automatically installed with install.sh

<string>Required databases</strong><br>
Taxdump database, download and extract the files in a desired location: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/ <br>
The variable NT_dmp in the settingsfile is where you specify where you extracted the dmp files to. <br>
NCBI NT database, specify the location in DATABASE variable in settingsfile. <br>
`update_blastdb.pl --decompress nt`

<strong>General usage</strong> </br>
The basic command to run the example analysis, after editing the settings in settingsfile.txt: 
```
bash START.sh settingsfile.txt
```
The previous command runs Guppy and the first clustering per sample.
To run the reclustering and Taxonomic identification:
```
bash START_recluster.sh settingsfile.txt
```
<strong>Modify the following parameters in settingsfile.txt to run the pipeline:</strong><br>
```<RunName>```             Run name of the analysis, this will be the name of the output folder <br>
```<mail-user>```           Email address to receive updates about the job and results from the analysis <br>
```<modus>```   "all" to analyse all datasets within OutDir or "one" to analyse dataset with same RunName and OutDir <br>
```<OutDir>```              Path where the output of the DNA metabarcoding analysis is placed <br>
```<workdir>```             Location of the Pimenta scripts on the HPC server <br>
```<THREADS>```             Total CPU threads <br>


<strong>Other options:</strong> <br>
```<RunModules>```      Modules to run, options: all (default), Guppy,Clustering, Consensus, Blast, Taxonomy. <br> Another module called 'oldmode' can be also be used, which runs the tax identification per sample (instead of per dataset). Needs to be used in combination with 'all' or other modules (e.g. "oldmode,clustering,Consensus,Blast,taxonomy") (experimental)<br> 

```<GPU>```             GPU usage (only for Guppy basecalling): '1' for using GPU or 'cpu' to run Guppy on CPU. <br>


<strong>Output files</strong> <br>
By default, the output of the DNAMetabarcoding pipeline is written to /```<OutDir>```/```<RunName>``` <br>
The taxonomy and QC output is written to /```<OutDir>```/```<RunName>```/Taxonomy_per_sample_```<ident2>```<br>


