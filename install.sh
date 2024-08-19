#!/bin/bash
set -e

ENVNAME="cenv"
BLASTVERSION="2.16.0"
BLASTFILE=ncbi-blast-$BLASTVERSION+-x64-linux.tar.gz


echo "Downloading $BLASTFILE"
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BLASTVERSION/$BLASTFILE -O ~/$BLASTFILE
tar zxvpf ~/$BLASTFILE -C ~

export PATH=$PATH:$HOME/ncbi-blast-$BLASTVERSION+/bin
PATHCMD="export PATH=\$PATH:$HOME/ncbi-blast-${BLASTVERSION}+/bin"
grep -qxF "$PATHCMD" $HOME/.bashrc || echo $PATHCMD >> $HOME/.bashrc

[ -d $HOME/blastdb ] || mkdir $HOME/blastdb
export BLASTDB=$HOME/blastdb
BLASTDBCMD="export BLASTDB=\$HOME/blastdb"
grep -qxF "$BLASTDBCMD" $HOME/.bashrc || echo $BLASTDBCMD >> $HOME/.bashrc


if ! command -v conda &> /dev/null
then
    echo "conda not found, installing miniconda"
    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm -rf ~/miniconda3/miniconda.sh
else
    echo "conda is already installed"
fi

echo "Setting up conda environment"
if conda info --envs | grep -q $ENVNAME; then echo "$ENVNAME environment already exists"; else conda create -y -n $ENVNAME; fi

source activate base
conda init bash
conda activate $ENVNAME
echo "conda environment activated"

echo "Installing required packages"
conda config --add channels bioconda

conda install -y pip
conda install -y parsnp

pip install biopython pandas xlsxwriter

echo "Installation complete"



