# pyRAP
Pipeline to annotate Rend-seq data genomewide

The pyRAP pipeline take wig files as input and outputs annotation files in GFF or GTF format.

The pipeline is written in Python3 and has the following dependencies:

- python3
- 
- 

All dependencies can be installed using conda. To install conda go to
https://www.anaconda.com/products/distribution
and download the installer appropriate for your operating system. Choose default options and default locations.

Open terminal and initialise conda by writing

$ conda init name_of_your_shell

 
The $ sign signifies a terminal command and shall not be included in your command line.

On Mac your shell is either bash or zsh

After initialising conda, restart terminal either be closing and reopening the window or by using the command

$ exec zsh -l

Use bash instead of zsh if you use a bash shell

After restarting your shell you will now see (base) written next to your username. To prevent conda from automatically activating base when opening terminal type this command in terminal

$ conda config --set auto_activate_base false


Now download the pyRAP pipeline from the GitHub repository https://github.com/ALawaetz/pyRAP

Move the folder to your desired location. Open terminal and navigate to the pyRAP folder using the command

$ cd path_to_pyRAP_folder

Make sure you are in the correct folder by typing

$ ls

You are in the right folder if you see among other files
-pipeline.py
-terminal_input.py
-Scripts

You will now create a conda environment called pyRAP and install all necessary dependencies by typing in terminal:

$ python3 install_dependencies.py

If that doesn't work try

$ python3 install_dependencies_alternative.py

If that doesn't work either, install dependencies manually. Dependencies are listed in the file named pyrap_environment.yml

To activate your environment type in terminal

$ conda activate pyRAP

In terminal you will see (pyRAP) next to your username. Every time you use pyRAP you need to activate your pyRAP environment.

Still being in your pyRAP folder you can run the pyRAP pipeline by typing in terminal

$ python3 pipeline.py -f FOLDER -c FILE --output_folder FOLDER

The pyRAP pipeline has a number of options. To view them all and for help of usage type:

$ python3 pipeline.py --help

The pyRAP package also includes a script named terminal_input.py which makes the use of options more easy. Open the script in your favourite editor (I recommend atom https://atom.io) and change the path variable names to match the paths on your computer. Then type in Terminal:

$ python3 terminal_input.py

If you want to experiment with different option parameters it is recommended that you use the --genome_start and --genome_end flags to analyse only a small region of your genome. Once you have optimised the parameters you can remove those flags to analyse the entire transcriptome which takes about 6 hours to finish. Alternatively, use default parameters and make genome-wide annotation files straight away. To increase the speed of the analysis, pyRAP uses multiprocessing. The number of processors to be used is set with --processors flag.

Annotation files can be viewed using various genome browsers. I recommend JBrowse which has a user-friendly desktop version that can be downloaded from https://jbrowse.org/jb2/download/
To get the browser running upload your fasta genome file along with .fai file. An .fai file can be created using samtools faidx which can be installed using conda.

With JBrowse you have the option of uploading annotation files in GFF or GTF format. Furthermore, you can visualise wig files using the Multi-wiggle track option. To upload wig files they need to be in the bigwig format. UCSC provides a tool to convert wig to BigWig called wigToBigWig and can be found following the links at http://hgdownload.soe.ucsc.edu/admin/exe/. 
It is convenient to visualise 5' and 3' reads from the same strand on the same track. That is done in JBrowse by uploading both files with the Multi-wiggle track option and then choosing xy-plot under Renderer type.


