#!/bin/bash
source ~/.bashrc
if ! command -v conda activate &> /dev/null
        then
        echo "conda not installed"
        exit 1
fi

conda install -y -c conda-forge mamba


mamba config --add channels bioconda && \
mamba config --add channels anaconda && \
mamba config --add channels r && \
mamba config --add channels defaults && \
mamba config --add channels conda-forge
mamba create -n pimenta -y -c bioconda
mamba install -n pimenta -y -c bioconda cd-hit
mamba install -n pimenta -y -c conda-forge r-base
mamba install -n pimenta -y -c bioconda perl-cgi perl-bioperl perl-cairo perl-statistics-basic  perl-json
#Before installing blastdb, install this:
#mamba install -n pimenta -y -c conda-forge -c bioconda -c defaults abricate

mamba install -n pimenta -y blast prinseq=0.20.4 cutadapt
mamba install -n pimenta -y mafft=7.471 openpyxl pandas
mamba install -n pimenta -y taxopy krona
conda activate pimenta
ktUpdateTaxonomy.sh

cd scripts/Porechop/
make
