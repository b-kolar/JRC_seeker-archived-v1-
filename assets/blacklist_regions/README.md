
#############################################
Blacklist Regions Data
hg19 Downloaded on May 23 2022
hg38 Downloaded on June 7 2022
Downloaded by Bronte Kolar
#############################################


###############################################################
File 1 - hg19_blacklist_regions.bed

list of blacklist regions for hg19 genome
concatenation (using cat command) of the following two BED files below:
	1) ENCFF001TDO.bed - DAC blacklisted regions (DBR)
	2) ENCFF001THR.bed - DUKE Excluded Regions (DER)

###############################################################

File 2 - hg38_blacklist_regions.bed

Old file name: hglft_genome_37ed9_f29ef0.bed 
UCSC liftover of blacklist_regions.bed (hg19) to hg38
Liftover done on June 7 2022
Parameters:
	- Allow multiple output regions: on
	- Minimum ratio of bases that must remap: on
	- Conversion failed on 16 records
	- Sucessfully converted 2044 records
	- List of failed record included (hglft_genome_37ed9_f29ef0.err)

Further information on these two files for the hg19_blacklist_regions.bed is listed below:

################################################################

1) DAC blacklisted regions (DBR)
https://www.encodeproject.org/annotations/ENCSR636HFF/

Summary for annotation file set ENCSR636HFF
Status
    released
Accession
    ENCSR636HFF
Description
    The DAC Blacklisted Regions aim to identify a comprehensive set of regions in the human genome that have anomalous, unstructured, high signal/read counts in next gen sequencing experiments independent of cell line and type of experiment. There were 80 open chromatin tracks (DNase and FAIRE datasets) and 20 ChIP-seq input/control tracks spanning ~60 human tissue types/cell lines in total used to identify these regions with signal artifacts. These regions tend to have a very high ratio of multi-mapping to unique mapping reads and high variance in mappability. Some of these regions overlap pathological repeat elements such as satellite, centromeric and telomeric repeats. However, simple mappability based filters do not account for most of these regions. Hence, it is recommended to use this blacklist alongside mappability filters. The DAC Blacklisted Regions track was generated for the ENCODE project.
Biosample summary
    (Homo sapiens)
Organism
    human
Annotation type
    blacklist

Encyclopedia version
    Blacklists
Lab
    Ewan Birney, EBI
Award
    U01HG004695 (Ewan Birney, EBI)
External resources

        UCSC-ENCODE-hg19:wgEncodeEH001432




2) Duke Excluded Regions (DER)
https://www.encodeproject.org/annotations/ENCSR797MUY/

Summary for annotation file set ENCSR797MUY
Status
    released
Accession
    ENCSR797MUY
Description
    Human Mapability Blacklist
Biosample summary
    (Homo sapiens)
Organism
    human
Annotation type
    blacklist

Encyclopedia version
    ENCODE2
Lab
    Gregory Crawford, Duke
Award
    U54HG004563 (Gregory Crawford, Duke)
External resources

        UCSC-ENCODE-hg19:wgEncodeEH000322


