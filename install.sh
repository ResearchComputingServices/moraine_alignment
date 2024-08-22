#!/bin/bash
set -ex

ENVNAME="cenv"
BLASTVERSION="2.16.0"
BLASTFILE=ncbi-blast-$BLASTVERSION+-x64-linux.tar.gz

# echo "Downloading $BLASTFILE"
# wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BLASTVERSION/$BLASTFILE -O ~/$BLASTFILE
# tar zxvpf ~/$BLASTFILE -C ~
# rm -rf ~/$BLASTFILE

# export PATH=$PATH:$HOME/ncbi-blast-$BLASTVERSION+/bin
# PATHCMD="export PATH=\$PATH:$HOME/ncbi-blast-${BLASTVERSION}+/bin"
# grep -qxF "$PATHCMD" $HOME/.bashrc || echo $PATHCMD >>$HOME/.bashrc

# [ -d $HOME/blastdb ] || mkdir $HOME/blastdb
# export BLASTDB=$HOME/blastdb
# BLASTDBCMD="export BLASTDB=\$HOME/blastdb"
# grep -qxF "$BLASTDBCMD" $HOME/.bashrc || echo $BLASTDBCMD >>$HOME/.bashrc

USERHOME="/home/$SUDO_USER"
if ! command -v conda &>/dev/null; then
    echo "conda not found, installing miniconda"
    mkdir -p $USERHOME/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $USERHOME/miniconda3/miniconda.sh
    bash $USERHOME/miniconda3/miniconda.sh -b -u -p $USERHOME/miniconda3
    rm -rf $USERHOME/miniconda3/miniconda.sh
    sudo -u $SUDO_USER $USERHOME/miniconda3/bin/conda init bash
else
    echo "conda is already installed"
fi

export PATH="$USERHOME/miniconda3/bin:$PATH"

source activate base
echo "Setting up conda environment"
conda info --envs
if conda info --envs | grep -q $USERHOME/.*/$ENVNAME; then echo "$ENVNAME environment already exists"; else conda create -y -n $ENVNAME; fi

conda activate $ENVNAME
echo "conda environment activated"

echo "Installing required packages"
conda config --add channels bioconda

conda install -y pip parsnp

pip install biopython pandas xlsxwriter

echo "Installation complete. Please reopen the terminal."
