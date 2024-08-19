# DNA sequence alignment research
This project runs parsnp to align a set of DNA sequences. 

# Requirements

This toolsuite requires both Parsnp and Blast to be installed in this machine. Information on how to install them individually can be found [here](#third-party-software).

# Install
You can install the required software with the following command:
```bash
sudo ./install.sh
```

# Run the application
The install script will create a conda environment. You can activate the environment with
```bash
conda activate cenv
```

You can find the list of arguments in the help menu of the main file:
```bash
python main.py -h
```

Make sure the genomes being tested are in the correct folders, then run the main file with the appropriate command:
```bash
python main.py ...
```

A detailed manual on how to use this software can be found in [here](#doc/manual.docx)

## <a name="thirdpartysoftware">Third party software</a>
Parsnp should be installed following the instructions in:
https://github.com/marbl/parsnp

It also requires Blast to run locally:
https://www.ncbi.nlm.nih.gov/books/NBK569861/