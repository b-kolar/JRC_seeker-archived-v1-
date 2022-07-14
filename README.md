

<img src="https://user-images.githubusercontent.com/52743495/173834089-526c540a-df4b-452f-964e-26104bf6f261.png" width="350" />

**JRC_seeker** is a tool for identifying jointly regulated CpGs (JRCs) in the human methylome. These regions can be classified into a number of genomic phenomenon, such as imprinting regions or methylation quantitative trait loci (mQTLs). Developed by the Genetic Identification Lab at Erasmus Medical Center, a BAM file of WGBS reads can be inputted and the result of this tool is a list of JRC locations and their associated p-values. JRC_seeker is built using a Snakemake pipeline that combines Python scripts with Linux shell commands. 

## Requirements
```
Operating system: tested on Ubuntu 18.04.6 LTS (Bionic Beaver)
R: tested on R version 3.6.1 (2019-07-05) -- "Action of the Toes"
Python: Python 3.9.12
RAM requirements: Do not attempt to run without at least 50 GB of RAM.
Runtime: Approximately 1 day for 98 GB BAM file and 3.1 GB reference genome, using 20 cores.
```

## Dependencies

### [conda](https://www.anaconda.com/products/individual)

Download and install conda if you do not have it already on your machine.
```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
bash Anaconda3-2021.11-Linux-x86_64.sh
```

### [mamba](https://github.com/mamba-org/mamba)

Install Mamba into your Conda-based python distribution
```
conda install -n base -c conda-forge mamba
```
Activate the Conda base environment (which now includes Mamba).
```
conda activate base
```


### [snakemake](https://snakemake.readthedocs.io/) (at least v4.3.1)

Create a new conda environment called ```jrc_seeker``` with snakemake and python 3.9 in it.
```
mamba create -c conda-forge -c bioconda -n jrc_seeker snakemake python=3.9
```

**Option:**
If conda-forge did not work for you, simply create a conda environment like this and use the pip command options below:
```
conda create -n jrc_seeker python=3.9 snakemake
```

Activate the ```jrc_seeker``` conda environment.
```
conda activate jrc_seeker
```
Check whether Snakemake is succesfully installed by running the following command:
```
snakemake --help
```

### [biopython-1.79](https://biopython.org/docs/1.79/api/Bio.html)

Install the following packages: biopython=1.79 .
```
conda install -c conda-forge biopython=1.79
```
or
```
pip install "biopython==1.79"
```

### [tabix](https://github.com/samtools/htslib)

```
conda install -c bioconda tabix
```
or
```
pip install tabix
```

### [pandas](https://pandas.pydata.org/)

```
conda install -c anaconda pandas
```
or
```
pip install pandas
```

### [bgzip](https://github.com/xbrianh/bgzip.git)

```
conda install -c bioconda bgzip
```
or
```
pip install bgzip
```

### [ChromHMM](http://compbio.mit.edu/ChromHMM/)

Quick instructions on running ChromHMM:

1. Install Java 1.5 or later if not already installed.
2. Download and unzip the ChromHMM.zip file using the following code snippit:

```
wget http://compbio.mit.edu/ChromHMM/ChromHMM.zip
unzip ChromHMM.zip
```

### [bedtools](https://bedtools.readthedocs.io/en/latest/)

```
conda install -c bioconda bedtools
```
or
```
pip install bedtools
```

### [samtools](http://www.htslib.org/doc/samtools.html)

Recommended version: 1.14

```
conda install -c bioconda samtools
```
or
```
pip install "samtools==1.14"
```

### [samblaster](https://github.com/GregoryFaust/samblaster)

```
conda install -c bioconda samblaster
```
or
```
pip install samblaster
```

### [BISCUIT](https://huishenlab.github.io/biscuit/)

```
conda install -c bioconda biscuit
```

### R

If not installed already, be sure you have R version 4.1.2 or greater.

### R dependencies

To install R dependencies, open R by running:
```bash
R
```
And run the following R commands:

```r
# To install
install.packages('data.table')
install.packages('parallel')
install.packages('MASS')

# To verify that installation was succesful
library(data.table)
library(parallel)
library(MASS)
```

## Downloading JRC_Seeker

Clone the repository with the following command:

```
mkdir jrc_seeker
cd jrc_seeker
git clone https://github.com/b-kolar/jrc_seeker.git
```

## Do a test run

Now that you have JRC_seeker installed and ready to go, let's run an example to see if the pipeline is fully working. Snakemake only works if you are in the same directory as the ```Snakefile```, so ensure you are in the ```/jrc_seeker``` directory where the ```Snakefile``` is.

In the ```/sample_data``` folder you can find some sample files. Be sure to edit the configuration file ```/sample_data/test_config.json``` with the correct paths. Later in this tutorial, you can find a detailed explantion of the contents of this sample data folder.

To run the example, run the following:

```
conda activate jrc_seeker
snakemake --cores 1 --configfile sample_data/test_config.json
```

On your screen, you will notice the steps of the pipeline running one-by-one. You can find the output of this test run in ```/sample_data/output``` and can compare your output to the expected output in the ```/sample_data/expected_output``` directory.

**Note**
You will notice that there is a warning in the ```label_states``` step stating: "ERROR: UNMETHYLATED STATE NOT CLASSIFIED CORRECTLY. RUN ON MORE DATA". The sample BAM file used does not contain enough data for ChromHMM to label all four states (methylated, unmethylated, intermediately methylated, no data). With this small dataset, the unmethylated state is not found. You can also see this below, as pictured in the ```/sample_data/expected_output/chromhmm/output_files/webpage_4.html``` file, where there is no unmethylated state in the emissions matrix of the ChromHMM output:

![image](https://user-images.githubusercontent.com/52743495/177798760-6eabb423-584f-4513-833d-17bd7b237157.png)

For example, in this image you can see that State 1 is intermediately methylated and State 2 is methylated, but there is no apparent unmethylated state.

Keep an eye on these warnings when running the pipeline and check the emissions matrix of the ChromHMM output to ensure that all states are deconvoluted.

## Process Overview
A number of steps are used to complete this analysis.

![image](https://user-images.githubusercontent.com/52743495/174032479-4844df52-63d8-4b46-90cc-15454b6cf113.png)


1. **Generate methylation data** - _pileup, meth_info, format_meth_data_
2. **Binarize methylation data** - _binarize_
3. **ChromHMM genome segmentation** - _learnmodel_
4. **BinPolish** - _binpolish_assets, label_states, binpolish_
5. **Binokulars** - _epiread, process_epiread, binokulars_

### Generate methylation data

BISCUIT is used to create a pileup VCF of DNA methylation and genetic information. After generating the VCF file, with both genetic and methylation information, beta values and coverage are extracted using the ```vcf2bed```  command to study the methylation levels at sequenced CpGs. The output of this step is a BED file of the format below, where the columns represent chromosome, CpG start base, CpG end base, beta methylation value (proportion of reads methylated), and total read coverage.

``` 
chr1    10469   10470   0.625   8
chr1    10471   10472   0.444   9
chr1    10484   10485   0.889   9
chr1    10489   10490   1.000   10
chr1    10493   10494   0.875   8
``` 
### Binarize methylation data

From the above BED file, the counts of methylated and unmethylated reads at every CpG location are calculated using the ```format_methylation.py``` script. Using the ```format_methylation.py``` script CpGs counts are combined using 200-bp bins and a binarized track is formed per chromosome for all chromosomes that contain methylation data:

``` 
methylated  unmethylated
0           1
0           0
1           1
1           0
1           1
``` 

The first column represents if the bin contains methylated reads ```(0 = False, 1 = True)``` and the second column represents if the bin contains unmethylated reads ```(0 = False, 1 = True)```. Thus: ```0   1``` is an unmethylated bin, ```1   0``` is a methylated bin, and ```1   1``` is an intermediately methylated bin. This bins are binarized from count data using their methylation value. A methylation value less than 0.2 is unmethylated and above 0.8 is methylated. Between these values is intermediately methylated. These binarization thresholds can be changed in the ```config.json``` file using the ```lower_im_methylation_bound``` and the ```upper_im_methylation_bound``` variables. The bin size can also be adjusted in the ```config.json``` file using the ```bin_size``` variable. ChromHMM, the genome segmentation software used in the following step recommends a bin size of 200. Changing this bin size here will also change it for the ChromHMM step. 

A threshold is also set to set bins with low coverage to a no-data state, aka ```0   0```. If there are less that 3 methylated/unmethylated counts (not CpGs), then this bin is set to the low-coverage state. This threshold can be set by changing the ```k_threshold``` variable in the ```config.json``` file.

### ChromHMM genome segmentation

ChromHMM, a genome segmentation and annotation software tool, is leveraged to identify intermediately methylated regions within the provided dataset. Ultimately, this step summarizes the binarized bins into larger regions throughout the genome. 

Finding intermediately methylated regions is necessary for finding JRCs. ChromHMM uses a multivariate hidden Markov Model to infer states from binarized data, allowing these regions to be identified and for the rest of the dataset/genome to be removed (as this is not of interest for the purpose of finding CpGs). ChromHMM segments all chromosomes provided into state regions and annotates them with a number. Each of these numbers corresponds to one of the four states of interest to us (no-data, unmethylated, methylated, intermediately methylated).

 ![image](https://user-images.githubusercontent.com/52743495/174035708-41f2e402-666b-48dd-845a-892f3d7194d9.png)

The output of ChromHMM is a list of segments that correspond to a given state. Using the emission matrix outputted by ChromHMM, the states can be annotated and the regions corresponding to intermediately methylated regions can be identified. 

### BinPolish

After the ChromHMM segmentation, there are often a large number of intermediately methylated regions, some of which are very small (200 bp). The goal of BinPolish is to remove some of these sparse, small regions or to summarize them into larger intermediately methylated blocks, thus reducing the large number of small intermediately methylated regions into a set of larger, robust ones. Essentially, we are polishing intermediately methylated segments.

![BinPolish](https://user-images.githubusercontent.com/52743495/174038503-de253ac7-e8c3-4e08-86a5-7d6f0d60af16.png)

Intermediately methylated regions are ignored that overlap with known regions that have low mappability (source) or are documented blacklisted regions (source). The number of CpGs per region is also calculated using the original data from BISCUIT. Following this, a the main steps employed that merge or remove intermediately methylated (IM) regions are:

1. Merge IM regions that are separated by 200bp
2. If small IM regions are within two larger regions, classify as the state of the earlier large region
3. If two IM regions are separated by a region that has less than one CpG, merge these IM regions
4. If IM region is classified wrongly, re-classify
5. Merge IM regions that are separated by 200bp (again)
6. If a state <=600bp is within two IM states that have higher cpg density, merge
7. States with CpG density <= 2 are turned to no-data states
8. If small IM regions are within two larger regions, classify as the state of the earlier large region (again)
9. Remove small regions (<= 200 bp)

### Binokulars

Using the set of intermediately methylated regions from BinPolish, binokulars assigns a p-value to each region, allowing users to identify significant JRCs. Binokulars uses the list of JRCs and a Single Fragment Epiread format file (produced by BISCUIT) to complete its analysis. 

Binokulars, named after its method (Binomial likelihood function-based bootstrap hypothesis test for co-methylation within reads), can be found [here](https://github.com/BenjaminPlanterose/Binokulars).

## Quality Control

When running JRC_seeker, you will notice that the labels of ChromHMM states are inferred in the snakemake step ```label_states```. In this step, the emissions file from ChromHMM is used to assign each state as representing intermediately methylated (IM), methylated(M), unmethylated (U), or no data (ND) regions. If there is not enough data, these states cannot be inferred correctly.

See below the result of a run with low coverage:

```
state  methylated  unmethylated label
0    E3    0.999999      0.970137    IM
1    E2    0.609009      0.006016     M
2    E4    0.126744      0.106534     U
3    E1    0.000658      0.000123    ND
```

You can notice that the state ```E4``` is labelled as unmethylated, but does not seem to represent this state (the emissions matrix shows that it is almost tied between methylated and unmethylated).

While it is most important that the IM states are inferred correctly, as these are retained, incorrect state labels can affect the accuracy of later steps like BinPolish. Be cautious and double check this matrix during the pipeline run. An error message will also be outputted if the states are not classified correctly.

As an example, a correctly labelled state matrix would be as follows:

```
state  methylated  unmethylated label
0    E3    0.91          0.97        IM
1    E2    0.63          0.01         M
2    E4    0.10          0.70         U
3    E1    0.0006        0.0001       ND
```

## Getting started
Now that you have JRC_seeker and its dependencies installed, you are ready to roll! To run with using JRC_seeker, there are three main steps to follow.

1. **Prepare input files** - Get your data ready for analysis and prepare input files in required formats.
2. **Edit the configuration file** - Adjust settings in the ```config.json``` file to suit your analysis needs.
3. **Run JRC_seeker** - Run the snakemake pipeline.

### Prepare input files
To use JRC_seeker, the following files are required:
1. BAM file
2. Reference Genome
3. Chromosomes file

If you are using a reference genome other than Hg19 or Hg38, these files will also be needed. Otherwise, they are provided in the ```/assets``` directory:

5. Mappability file
6. Blacklist region file

#### BAM 
A sorted, indexed BAM with duplicate marked reads is required to run this analysis. We used BISCUIT to align our data and produced an index (.csi) for the BAM file using:

```samtools sort --write-index -o my_output.bam -O BAM -```

#### Reference Genome

A reference genome is required to run JRC_seeker and JRC_seeker has asset files supporting the use of the hg19 and hg38 reference genomes.

The reference genome you use must be indexed using the BISCUIT ```index``` command. To download and index a reference genome, you can run the following command (for hg19):

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz
gunzip hg19.fa.gz
biscuit index hg19.fa
```

Be aware that indexing a reference genome can take time, but only needs to be done once for each reference.

#### Chromosomes file

The Snakemake pipeline requires a text file that contains a list of chromosomes found in the BAM file. Depending on how you intend to use JRC_seeker, you can create this file two different ways:

###### #1 - Process entire BAM

Generally, you will run your entire BAM through JRC_seeker. In this case, do the following:

To generate the chromosomes file, ChromHMM must also be installed and you must have the location to the CHROMSIZES file corresponding to your genome version (e.g. hg19.txt). This file can be created using the following Linux shell commands:

```
samtools idxstats your_file.bam | grep 'chr' |  cut -f 1 | tee temp.txt
awk 'NR==FNR{A[$1];next}$1 in A' /ChromHMM/CHROMSIZES/hg19.txt temp.txt > chromosomes.txt
rm temp.txt
```

For this run, remember to set your the ```region``` parameter in your config file to ```none```.

###### #2 - Process specific region

If you want to process a single region, you should manually create a chromosomes.txt file. In this case, just save the chromosome number in a text file titles ```chromosomes.txt```.

For example, if your region of interest is ```chr1:3000-50000```, your chromosomes.txt file should be:

```
chr1
```

For this run, remember to set your the ```region``` parameter in your config file to the region you are investigating. In the above example, you would set the ```region``` parameter in your config file to ```chr1:3000-50000```.

Please be aware that even if you have a reference genome for a single chromosome and you have extracted reads from a single chromsome (as is with chromosome 21 in the JRC_seeker example), you must still treat this as a run that is processing as specific region.

#### Asset files

##### Mappability Files

Regions with low mappability can result in false positive JRC regions, which is why regions with low mappability are removed during BinPolish. In the ```/assets/mappability_files``` directory, two files are included (for hg19 and hg38) that contain a list of regions of the genome that are uniquely mappable by at least one k-mer (in our case, k=100).

These files are the Bismap individual k-mer files for Human hg19 and Human hg38 genomes (k100 Single-read), which were downloaded from: https://bismap.hoffmanlab.org/

##### Blacklist Regions

Regions that overlap with regions have been labelled as "blacklist" regions are also removed during BinPolish. For hg19, a blacklist file was created by combining regions from the following two blacklist region datasets:

1) DAC blacklisted regions ([DBR - ENCFF001TDO.bed](https://www.encodeproject.org/annotations/ENCSR636HFF/))

2) Duke Excluded Regions ([DER - ENCFF001THR.bed](https://www.encodeproject.org/annotations/ENCSR797MUY/))

For hg38, the hg19 file was lifted over using [UCSC LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver).

Further details on the LiftOver settings and the blacklist regions files are listed in the ```/assets/blacklist_regions/README.md```.

### Edit the configuration file

Snakemake uses a configuration file to locate external files and parameter values. Make a copy of the sample configuration file ```/sample_data/test_config.json``` and edit the file and directory paths to point towards the respective input files and directories on your machine (use absolute paths). Furthermore, be sure to edit the parameters in the configuration file. 

The configuration file is organized into three sections, for readability:
_DIRECTORIES_: Directories to add your own path to.
_FILES_: Files to add your own path to.
_PARAMETERS_: Parameters to change. The listed default values for parameters that are listed at the end of this tutorial are the suggested values per parameter. 

Here are a couple of important reminders:
- If running JRC_seeker on your entire BAM file, be sure to change the region value to none:
```"region" : "none"```
- Ensure you change the number of binokulars cores to an amount your machine can handle, such as 4:
```"binokulars_cores" : 4```

A detailed list of the parameters in the configuration file are found at the end of this tutorial.

### Run JRC_seeker

Go into the JRC_seeker directory by:
```
cd <Path to JRC_seeker>
```

Make sure to activate the conda environment by: 
```
conda activate jrc_seeker
```

To make sure there are no problems with your configuration file, do a dry run of the snakemake pipeline:

```
snakemake -n --configfile [path to config file]
```

To run the JRC_Seeker pipeline, simply run the command below, where the number of cores and the path to your edited configuration file is specified.
```
snakemake --cores [amount of cores] --configfile [path to config file]
```

*Note:* It is recommended to run with 1 core and specify in the ```config.json``` file the number of cores to run the binokulars step with. If the number of snakemake cores is above 1, the number of cores the final binokulars step uses will be the product of the snakemake cores and the value of the ```binokulars_cores``` field in the ```config.json``` file.

## A note on time 

The final step of the snakemake pipeline is the actual binokulars run. From a 98 GB BAM file of pooled blood samples of 12 individuals and a 3.1 GB reference genome (entire hg19 reference genome), ~400,000 regions are found to be intermediately methylated regions after BinPolish. It takes approximately 2 seconds for binokulars to process each region. As a result, for 4 cores we expect a binokulars run to take ~55 hours. Using 20 cores, this takes ~11 hours. 

The binokulars step is the bottleneck of this pipeline, but it is also the most parallelizable. The remainder of the pipeline takes under 10 hours for the above mentioned dataset. We recommend using ~20 cores to run the entire pipeline within a day.

## Overview of sample data

In the ```/sample_data``` folder you can find some sample files to test if JRC Seeker is running properly on your machine. Below you can find an explanation of where these files come from and how we created them, just in case you were interested:

**sample_data.bam - Sample BAM file**

This BAM file is of data from chromosome 20 of the pooled blood data from old and young men from the following study: https://www.ebi.ac.uk/ena/browser/view/PRJEB28044?show=reads.

```
samtools view -b pooled_blood_hg38.bam chr20 > chr_20.bam
```

To reduce the file size, The first million aligned reads were selected using the following command:

```
samtools view -h chr_20.bam \
    | head -n 1000000 \
    | samtools view -b -o sample_data.bam
```

**sample_data.bam.csi - BAM index**

The BAM file was indexed using the following command:

```
samtools index -c sample_data.bam
```

**chromosomes.txt - chromosomes list**

A list of chromosmes, as outlined earlier. As this comprised of only data from chromosome 20, the chromosomes.txt file is simply one line:

```
chr20
```

**test_config.json - sample config file**

The sample config file contains the settings for this run. As we are only looking at chromosome 20, the region is set to ```"chr20"```.

**Reference genome**

In the ```/reference_genome``` directory, you can find ```chr20.fa```, a FASTA file of the reference genome of hg38 for chromosomome 20, which was downloaded [from this Genome Reference Consortium source](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/) using the following command: 

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz
gunzip chr20.fa.gz
```

The FASTA file was indexed by BISCUIT using the following command, which generated the additional files:

```
biscuit index chr20.fa
```

## Config file parameters

DIRECTORIES:
```
output_folder : folder all files will be outputted in (absolute path)

path_to_scripts : path to JRC_seeker scripts folder

temp_folder : some computations might be too large for default temp folders. Put the default temp folder path here to add a new path (absolute path)

path_to_jrc_seeker : absolute path to JRC_seeker (e.g. /home/opt/jrc_seeker). This is the folder where the Snakefile is located.
```

FILES:
```
path_to_config_file : path to Snakemake config.json file (absolute path)

path_to_reference_genome : path to FASTA reference genome (absolute path)

path_to_bam : path to BAM file (absolute path)

path_to_chrom_length_file : path to chromosome length text file used by ChromHMM. Found in the CHROMSIZES folder of ChromHMM. Use the file corresponding the the genome version your are using (e.g. hg19.txt) (absolute path)

chromhmm : path to ChromHMM jar file (absolute path)

path_to_mappability_file : path to mappability file. For hg19 or hg38, use the ones in the assets/mappability_files folder. Otherwise, add the absolute path to the version for your genome.

chromosomes_file : path to chromosomes.txt file outlined in the "Input Files" step of this tutorial (absolute path)

path_to_blacklist : path to blacklist regions file. For hg19 or hg38, use the ones in the assets/blacklist_regions folder. Otherwise, add the absolute path to the version for your genome.
```

PARAMETERS:
```
sample_name : name of your sample (don't use spaces or the following characters: "/" "," "." "\")

bin_size : size of bins for binzarization and ChromHMM segmentation. Measured in base pairs. Recommended not to change. (default: 200)

binokulars_output_directory : name of directory for binokulars output (don't use spaces or the following characters: "/" "," "." "\")

region : if a specific region of the BAM file is to be investigated OR if the BAM file contains one chromosome, specify this here (e.g. "chr1" or "chr20:0-1000"). Only one region is permitted and if changing this setting, ensure that the chromosomes.txt file is manually created (see instructions above). In most circumstances, the entire BAM file is to be processed and this should thus be set to "none". (default: "none")

lower_im_methylation_bound : intermediately methylated methylation value for lower boundary for binarization. Recommended not to change (default: 0.2)

upper_im_methylation_bound : intermediately methylated methylation value for upper boundary for binarization. Recommended not to change (default: 0.8)

data_assignment_setting : if (un)methylated counts are on the bin boundary, they are added to the earlier or "left" bin. Recommended not to change (default: "left", "right" is other option)

k_threshold : threshold for number of reads needed for bin to not be set to a no-data state (see "Binarize methylation data" section). (default: 3)

n_states : number of states for ChromHMM to predict. Binarize subroutine depends on 4 states. Recommended not to change. (default: "4")

chromhmm_it : number of ChromHMM iterations. Recommend 500 to ensure convergence, but ChromHMM default for maxiterations in the LearnModel command is 200. (default: "500")

map_threshold : threshold overlap coverage (as a percent) of intermediately methylated regions overlapping with low mappability regions for them to be discarded. Recommend not to change (default: 0.95)

segment_min_sz : BinPolish discards regions equal or smaller to this threshold size. (default: 200)

permutation_iterations : number of permutations that binokulars runs per region. (default: 1000)

seed_binokulars : binokulars seed value. (default: 4)

binokulars_cores : number of cores that the binokulars subroutine uses. (default: 4)

flank_length : number of base pairs that binokulars flanks regions by. Recommend not to change (default: 500)
```

## Sources

ChromHMM: automating chromatin-state discovery and characterization, Nature methods, 2012, Jason Ernst & Manolis Kellis

Faust, G.G. and Hall, I.M., “SAMBLASTER: fast duplicate marking and structural variant read extraction,” Bioinformatics Sept. 2014; 30(17): 2503-2505.

Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z

Karimzadeh M, Ernst C, Kundaje A, Hoffman MM. 2018. Umap and Bismap: quantifying genome and methylome mappability. doi: https://doi.org/10.1093/nar/gky677 Nucleic Acids Research, Volume 46, Issue 20, 16 November 2018, Page e120.
