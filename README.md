

<img src="https://user-images.githubusercontent.com/52743495/173834089-526c540a-df4b-452f-964e-26104bf6f261.png" width="350" />

**JRC Seeker** is a tool for identifying jointly regulated CpGs (JRCs) in the human methylome. These regions can be classified into a number of genomic phenomenon, such as imprinting regions or methylation quantitative trait loci (mQTLs). Developed by the Genetic Identification Lab at Erasmus Medical Center, a BAM file of WGBS reads can be inputted and the result of this tool is a list of JRC locations and their associated p-values.

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
3. If two IM regions are separated by a 400bp region that has less than two CpGs, merge these IM regions
4. Merge IM regions that are separated by 200bp (again)
5. States with CpG density <= 2 are turned to no-data states
6. If small IM regions are within two larger regions, classify as the state of the earlier large region (again)
7. Remove small regions (<= 200 bp)

### Binokulars

Using the set of intermediately methylated regions from BinPolish, Binokulars assigns a p-value to each region, allowing users to identify significant JRCs.
Binokulars uses the list of JRCs and a Single Fragment Epiread format file (produced by BISCUIT) to complete its analysis. 

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
JRC Seeker is built using a Snakemake pipeline that combines Python scripts with Linux shell commands. To get started with using JRC Seeker, there are three main steps to follow.

1. **Prepare input files** - Get your data ready for analysis and prepare input files in required formats.
2. **Download JRC Seeker** - Download the source code of this project and install dependencies.
3. **Edit the configuation file** - Adjust settings in the ```config.json``` file to suit your analysis needs.
4. **Run JRC Seeker** - Run the snakemake pipeline.

### Input files
To use JRC Seeker, the following files are required:
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

A referenc genome is required to run JRC_seeker and JRC_seeker has asset files supporting the use of the hg19 and hg38 reference genomes.

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

These files are the Bismap individual k-mer files for Human hg19 and Human hg38 genomes (k100 Single-read), which were downloaded from:

https://bismap.hoffmanlab.org/

##### Blacklist Regions

Regions that overlap with regions have been labelled as "blacklist" regions are also removed during BinPolish. For hg19, a blacklist file was created by combining regions from the following two blacklist region datasets:

1) DAC blacklisted regions ([DBR - ENCFF001TDO.bed](https://www.encodeproject.org/annotations/ENCSR636HFF/))

2) Duke Excluded Regions ([DER - ENCFF001THR.bed](https://www.encodeproject.org/annotations/ENCSR797MUY/))

For hg38, the hg19 file was lifted over using [UCSC LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver).

Further details on the LiftOver settings and the blacklist regions files are listed in the ```/assets/blacklist_regions/README.md```.

## Downloading JRC Seeker

Clone the repository with the following command:

```
$ git clone https://github.com/b-kolar/binokulars.git
```
Enter the repository by cd binokulars.

## Install dependencies

### [binokulars](https://github.com/BenjaminPlanterose/Binokulars)

Follow the installation instructions outlined on: https://github.com/BenjaminPlanterose/Binokulars.

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

### [samtools](http://www.htslib.org/doc/samtools.html)

```
conda install -c bioconda samtools
```

### [samblaster](https://github.com/GregoryFaust/samblaster)

```
conda install -c bioconda samblaster
```

### [BISCUIT](https://huishenlab.github.io/biscuit/)

```
conda install -c bioconda biscuit
```

## Set up conda environment
- Download [Conda](https://www.anaconda.com/products/individual) with Python 3.9

```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
```

Install Conda
```
bash Anaconda3-2021.11-Linux-x86_64.sh
```

- [Snakemake](https://snakemake.readthedocs.io/) (at least v4.3.1) and [biopython-1.79](https://biopython.org/docs/1.79/api/Bio.html)

Install Mamba into your Conda-based python distribution
```
conda install -n base -c conda-forge mamba
```
Activate the Conda base environment (which now includes Mamba).
```
conda activate base
```
Create a new conda environment called ```jrc_seeker``` with snakemake in it.
```
mamba create -c conda-forge -c bioconda -n snakemake jrc_seeker
```
Activate the ```jrc_seeker``` conda environment.
```
conda activate jrc_seeker
```
Check whether Snakemake is succesfully installed by running the following command:
```
snakemake --help
```
Install the following packages: biopython=1.79 .
```
conda install -c conda-forge biopython=1.79
```

## Running JRC Seeker
Go into the JRC Seeker directory by:
```
cd <Path to JRC Seeker>
```
Make sure to activate the conda environment by: 
```
conda activate jrc_seeker
```
```
snakemake --cores [amount of cores] --configfile [path to config file]
```

### Test run
In order to perform a testrun of JRC Seeker go into the JRC Seeker directory and execute the following command:
```
snakemake --cores [amount of cores] --configfile test_data/config_test.json
```

It is recommended to run with 1 core and specify in the ```config.json``` file the number of cores to run the binokulars step with. If the number of snakemake cores is above 1, the number of cores the final binokulars step uses will be the product of the snakemake cores and the value of the ```binokulars_cores``` field in the ```config.json``` file.


## Config file

Paths to change:

```
output_folder : folder all files will be outputted in (absolute path)

path_to_config_file : path to Snakemake config.json file (absolute path)

path_to_scripts : path to Binokulars scripts folder

path_to_reference_genome : path to FASTA reference genome (absolute path)

path_to_bam : path to BAM file (absolute path)

sample_name : name of your sample (don't use spaces or the following characters: "/" "," "." "\")

temp_folder : some computations might be too large for default temp folders. Put the default temp folder path here to add a new path (absolute path)

path_to_chrom_length_file : path to chromosome length text file used by ChromHMM. Found in the CHROMSIZES folder of ChromHMM. Use the file 
corresponding the the genome version your are using (e.g. hg19.txt) (absolute path)

chromhmm : path to ChromHMM jar file (absolute path)

path_to_mappability_file : path to mappability file. For hg19 or hg38, use the ones in the assets/mappability_files folder. Otherwise, add the absolute path to the version for your genome.

chromosomes_file : path to chromosomes.txt file outlined in the "Input Files" step of this tutorial (absolute path)

binokulars_output_directory : name of directory for binokulars output (don't use spaces or the following characters: "/" "," "." "\")

path_to_blacklist : path to blacklist regions file. For hg19 or hg38, use the ones in the assets/blacklist_regions folder. Otherwise, add the absolute path to the version for your genome.
```

Model parameters:

```
bin_size : size of bins for binzarization and ChromHMM segmentation. Measured in base pairs. Recommended not to change. (default: 200)

lower_im_methylation_bound : intermediately methylated methylation value for lower boundary for binarization. Recommended not to change (default: 0.2)

upper_im_methylation_bound : intermediately methylated methylation value for upper boundary for binarization. Recommended not to change (default: 0.8)

data_assignment_setting : if (un)methylated counts are on the bin boundary, they are added to the earlier or "left" bin. Recommended not to change (default: "left", "right" is other option)

k_threshold : threshold for number of reads needed for bin to not be set to a no-data state (see "Binarize methylation data" section). (default: 3)

n_states : number of states for ChromHMM to predict. Binarize subroutine depends on 4 states. Recommended not to change. (default: "4")

map_threshold : threshold overlap coverage (as a percent) of intermediately methylated regions overlapping with low mappability regions for them to be discarded. Recommend not to change (default: 0.95)

segment_min_sz : BinPolish discards regions equal or smaller to this threshold size. (default: 200)

permutation_iterations : number of permutations that binokulars runs per region. (default: 1000)

seed_binokulars : binokulars seed value. (default: 4)

binokulars_cores : number of cores that the binokulars subroutine uses. (default: 4)

flank_length : number of base pairs that binokulars flanks regions by. Recommend not to change (default: 500)
```

## Sample Data

In the ```/sample_data``` folder you can find some sample files to get binokulars running. Below you can find an explanation of where these files come from and how we created them, just in case you were interested:

**sample.bam - Sample BAM file**

This BAM file is of data from chromosome 20 of the pooled blood data from old and young men from the following study: https://www.ebi.ac.uk/ena/browser/view/PRJEB28044?show=reads.

```
samtools view -b pooled_blood_hg38.bam chr20 > chr_20.bam
```

To reduce the file size, 10% of aligned reads were selected using the following command:

```
samtools view -bo sample.bam -s 123.1 chr_20.bam 
```

**sample.bam.csi - BAM index**

The BAM file was indexed using the following command:

```
samtools index -c sample.bam
```

**chromosomes.txt - chromosomes list**

A list of chromosmes, as outlined earlier. As this comprised of only data from chromosome 21, the chromosomes.txt file is simply one line:

```
chr21
```

**test_config.json - sample config file**

The sample config file contains the settings for this run. As we are only looking at chromosome 21, the region is set to ```"chr21"```.
**Reference genome**

In the ```/reference_genome``` directory, you can find ```chr21.fa```, a FASTA file of the reference genome of Hg38 for chromosomome 21, which was downloaded from: https://www.ncbi.nlm.nih.gov/nuccore/CM000683.2?report=fasta

The FASTA file was indexed by BISCUIT using the following command, which generated the additional files:

```
biscuit index chr21.fa
```

Reference genomes are quite heavy files. To get set up with a reference genome, you can download version Hg19 and index it with BISCUIT using the code below. Be aware that indexing a reference genome can take time, but only needs to be done once.

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz
biscuit index hg19.fa.gz
```

## Time 

The final step of the snakemake pipeline is the actual binokulars run. To process ~400,000 regions from pooled blood samples of 12 individuals, it takes approximately 2 seconds per site. As a result, for 4 cores we expect a binokulars run to take ~55 hours. Using 20 cores, this takes ~11 hours. 

The binokulars step is the bottleneck of this pipeline, but it is also the most parallelizable. The remainder of the pipeline takes under 10 hours for the above mentioned dataset. We recommend using ~20 cores to run the entire pipeline within a day.
 

## Sources

ChromHMM: automating chromatin-state discovery and characterization, Nature methods, 2012, Jason Ernst & Manolis Kellis

Faust, G.G. and Hall, I.M., “SAMBLASTER: fast duplicate marking and structural variant read extraction,” Bioinformatics Sept. 2014; 30(17): 2503-2505.

Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z

Karimzadeh M, Ernst C, Kundaje A, Hoffman MM. 2018. Umap and Bismap: quantifying genome and methylome mappability. doi: https://doi.org/10.1093/nar/gky677 Nucleic Acids Research, Volume 46, Issue 20, 16 November 2018, Page e120.

# cite binokulars!!!
