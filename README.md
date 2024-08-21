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

<details>
<summary><h2 style="display:inline-block">Running this code on a virtual machine (VM)?</h2></summary>

The time to process this code can become long. When running this code on a VM, you may stop the code from running if you disconnect from the server, either voluntarely or involutarely.

To run the application in the background and detach it from the current terminal session, you can use the `nohup` command. Here's how you can modify the command to use `nohup`:

```bash
nohup python main.py ...
```

This will prevent the application from being terminated when you close the terminal session. The output of the application will be redirected to a file named `nohup.out` in the current directory.

Remember to replace `...` with the appropriate command-line arguments for your application.

You can check the progress and any error messages by viewing the `nohup.out` file using the `tail` command:

```bash
tail -f nohup.out
```

Make sure to monitor the progress of the application and check the `nohup.out` file periodically.

For more information on using `nohup`, you can refer to the [man page](https://man7.org/linux/man-pages/man1/nohup.1.html).

</details>


## <a name="thirdpartysoftware">Third party software</a>
Parsnp should be installed following the instructions in:
https://github.com/marbl/parsnp

It also requires Blast to run locally:
https://www.ncbi.nlm.nih.gov/books/NBK569861/