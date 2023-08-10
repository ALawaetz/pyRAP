# pyRAP
Pipeline to annotate Rend-seq data genomewide

The pyRAP pipeline take wig files as input and outputs annotation files in GFF3 or GTF format.

## Installation
The pipeline is written completely in Python 3. Dependencies can be found in the file pyrap_environment.yml. For installation of dependencies, see below.

All dependencies can be installed using conda. To install conda go to
<https://www.anaconda.com/products/distribution>
and download the installer appropriate for your operating system. Choose default options and default locations.

Put conda in your PATH by typing in terminal

```
export PATH="/home/username/anaconda/bin:$PATH"
```

If you want conda to be in your path each time you open terminal then put the command in the .bashrc or .zshrc file on your system.

Open terminal and initialise conda by writing

```
conda init name_of_your_shell
```

On Mac your shell is either bash or zsh

After initialising conda, restart terminal either be closing and reopening the window or by using the command

```
exec zsh -l
```

Use bash instead of zsh if you use a bash shell

After restarting your shell you will now see (base) written next to your username. To prevent conda from automatically activating base when opening terminal type this command in terminal

```
conda config --set auto_activate_base false
```

Now download the pyRAP pipeline from the GitHub repository <https://github.com/ALawaetz/pyRAP>

Move the folder to your desired location. Open terminal and navigate to the pyRAP folder using the command

```
cd path_to_pyRAP_folder
```

You are in the right folder if you see among other files ```pipeline.py``` and ```Example_Terminal_input.py```

You can now create a conda environment called pyRAP and install all necessary dependencies by typing in terminal:

```
python3 install_dependencies.py
```

Or alternatively

```
python3 install_dependencies_alternative.py
```

Dependencies can also be installed manually. Dependencies are listed in the file named ```pyrap_environment.yml```

## Tutorial
Running the pipeline can be done using the wig files found in the folder ```Wig_files```. For running the pipeline on other datasets replace the files with your wig files of interest. For manual entry of condition names using the ```-n``` flag, beware to also change the names in the txt file found in the folder ```Txt_files_examples```.

First activate your environment by typing in terminal

```
conda activate pyRAP
```

In terminal you will see (pyRAP) next to your username. Every time you use pyRAP you need to activate your pyRAP environment.

Still being in your pyRAP folder you can run the pyRAP pipeline by typing in terminal

```
python3 pipeline.py -f FOLDER -c FILE --output_folder FOLDER
```

Three inputs are mandatory:  
1. ```-f``` FOLDER (path to wig folder)  
2. ```-c``` FILE (chrom size file, see folder ```chrom_size_files```)  
3. ```--output_folder``` FOLDER (path to output folder)

The pyRAP pipeline has a number of options. To view them all and for help of usage type:

```
python3 pipeline.py --help
```

The pyRAP package also includes a script named ```EXAMPLE_Terminal_input.py``` which makes the use of options more easy. Open the script in your favorite editor (for example [Visual Studio Code] (https://code.visualstudio.com)) to carry out custom analysis.

```
python3 Example_Terminal_input.py
```

If you want to experiment with different option parameters it is recommended that you use the ```--genome_start``` and ```--genome_end``` flags to analyse only a small region of your genome. Once you have optimised the parameters you can remove those flags to analyse the entire transcriptome which takes about 1 hour to finish pr condition. Alternatively, use default parameters and make genome-wide annotation files straight away. To increase the speed of the analysis, pyRAP uses multiprocessing. The number of processors to be used is set with the ```--processors``` flag.

Annotation files can be viewed using various genome browsers. We recommend [JBrowse2] (https://jbrowse.org/jb2/download/) which has a user-friendly desktop version.
To get the browser running, upload your fasta genome file along with .fai file. An .fai file can be created using ```samtools faidx``` which can be installed using conda.

With JBrowse you have the option of uploading annotation files in GFF or GTF format. Furthermore, you can visualise wig files using the Multi-wiggle track option. To upload wig files they need to be in the bigwig format. UCSC provides a tool to convert wig to BigWig called wigToBigWig which can be downloaded following the links at <http://hgdownload.soe.ucsc.edu/admin/exe/>. The script ```wigToBigWigBatch.py``` in the main directory can then be used together with the wigToBigWig programme for batch conversion.
 
It is convenient to visualise 5' and 3' reads from the same strand on the same track. That is done in JBrowse by uploading both files with the Multi-wiggle track option and then choosing xy-plot under Renderer type.

Annotation files obtained from alternative sources (e.g. tiling array data, Subtiwiki, EcoCyc, AureoWiki) can be used to transfer features and names to Rend-seq annotation and for statistical analysis/comparison. Alternative annotation files used to verify pyRAP can be found in the folder ```alternative_annotation_files```.

Operon data can be generated by pyRAP automatically (see method section of publication), or it can be supplied as a GFF3 file. Operon files used to verify pyRAP are found in the folder ```alternative_operons```.

Genomic coordinates of start and stop sites and operon boundaries determined by different methods (e.g. Cappable-seq, Term-seq or tiling arrays) can be used to verify low-confidence peaks (see method section of publication). Alternative coordinates used in this study are found in the folder ```alternative_coordinates```.

Naming of files is done based on the names of the input files. If other names are to be used, those names can be provided with the ```-n``` flag which takes txt files as input. pyRAP also includes a statistics module where the ```--feature_distance``` flag defines the features on which statistics is to be carried out. Differentially expressed peaks can be found using the ```--compare_peaks_across_samples``` flag which takes a txt file as input denoting what conditions should be compared. Examples are found in the folder ```Txt_files_examples```.

pyRAP needs to know the name of the genome version and the size of the genome for the wig files analysed. This file is provided by the flag ```-c``` flag and the input is a txt file (see example in folder ```Chrom_size_files```). 

The output of pyRAP is annotation files in GFF3 or GTF format. The location of output files are set with ```---output_folder``` and the format is set with ```--output```.

