# Install R and RStudio

* Install the latest version of R, you can download it from [here](https://www.r-project.org/).
* Install the latest version of RStudio Desktop, you can download it from [Posit](https://posit.co/download/rstudio-desktop/).

# Install Bioconductor

Open RStudio, and run the following in Console:

```r
install.packages("BiocManager");
library(BiocManager);
BiocManager::install();
```

You might be prompted with `Update all/some/none? [a/s/n]:`, type `a` and press `enter`.

# Install Windows Subsystem for Linux (WSL) - Windows User Only

Follow this [Microsoft Tutorial](https://learn.microsoft.com/en-us/windows/wsl/install) to install WSL. The `<Distribution Name>` should be `Ubuntu`. 

# Install Miniconda in WSL - Windows User Only 

Open Ubuntu. In the terminal, run the following:

```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
```

After download complete, run:

```sh
sh Miniconda3-latest-Linux-x86_64.sh
```

You will be prompted with the following, press `enter`, then press `space` multiple times to read through the license:

```
Welcome to Miniconda3 py312_24.11.1-0

In order to continue the installation process, please review the license
agreement.
Please, press ENTER to continue
>>>
```

You will be prompted with the following, type `yes` and press `enter`:

```
Do you accept the license terms? [yes|no]
>>>
```

You will be prompted with the following, press `enter`:

```
Miniconda3 will now be installed into this location:
/home/jiajia/miniconda3

  - Press ENTER to confirm the location
  - Press CTRL-C to abort the installation
  - Or specify a different location below

[/home/jiajia/miniconda3] >>>
```

You will be prompted with the following, type `yes` and press `enter`:

```
Do you wish to update your shell profile to automatically initialize conda?
This will activate conda on startup and change the command prompt when activated.
If you'd prefer that conda's base environment not be activated on startup,
   run the following command when conda is activated:

conda config --set auto_activate_base false

You can undo this by running `conda init --reverse $SHELL`? [yes|no]
[no] >>> 
```

You will see `Thank you for installing Miniconda3!`, this means you have installed Miniconda. 

# Install Miniconda in Mac - Mac User Only

Go to this address https://www.anaconda.com/download/success, scroll down and make sure you are in the section of __Miniconda Installers__ rather than __Anaconda Installers__. Pick the right installer for your CPU. You can choose either Graphical Installer or Command Line Installer.

* If you choose Graphical Installer, just follow the instructions to install.
* If you choose the Command Line Installer:
    * Download the right installer.
    * Open Terminal app on your Mac.
    * Go to the folder where the installer was downloaded:
        * `cd ~/Downloads`
    * Run the installer script:
        * `bash Miniconda3-latest-MacOSX-arm64.sh`
        * The name of the `.sh` file might be different, make sure you type the correct name. 
    * Follow the on-screen prompts:
        * Press `enter` to review the license, and press `space` to scroll down the license.
        * Type `yes` to accept the terms.
        * Confirm the installation path.
        * Type `yes` to update the shell profile. 
    * You will see `Thank you for installing Miniconda3!`, this means you have installed Miniconda. 


# Install Integrative Genomics Viewer (IGV)

You can download it from [here](https://igv.org/doc/desktop/#DownloadPage/), pick the one with Java included. 

# Install SSH client (optional)

# Install Loupe browser (check later)

# Download scRNAseq files (check later)

From [here](https://genomedata.org/rnaseq-tutorial/scrna/).

# Reading Materials

Please read these articles before coming to the workshop (???):

https://rnabio.org/module-00-setup/0000/05/01/Prerequisites/ 