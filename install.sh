#!/bin/bash
source ~/.bashrc
if ! command -v conda activate &> /dev/null
        then
        echo "conda not installed"
        exit 1
fi

conda install -y -c conda-forge mamba 


conda config --add channels bioconda && \
conda config --add channels anaconda && \
conda config --add channels r && \
conda config --add channels defaults && \
conda config --add channels conda-forge
mamba create -n pimenta -y -c bioconda python=3.10
mamba install -n pimenta -y -c bioconda cd-hit
mamba install -n pimenta -y -c conda-forge r-base
mamba install -n pimenta -y -c bioconda perl-cgi perl-bioperl perl-cairo perl-statistics-basic  perl-json
#Before installing blastdb, install this:
#mamba install -n pimenta -y -c conda-forge -c bioconda -c defaults abricate

mamba install -n pimenta -y blast prinseq=0.20.4 cutadapt
mamba install -n pimenta -y mafft=7.471 openpyxl pandas
mamba install -n pimenta -y -c bioconda taxopy krona
conda activate pimenta
ktUpdateTaxonomy.sh

