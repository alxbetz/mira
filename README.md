# MiRa
Microarray Reanalysis is an R package that enables the re-analysis of microarray datasets in the Gene Expression Omnibus(GEO) using up-to-date transcriptome annotations.

## Installation

### Requirements

* MacOS 10.13.6+ or any flavor of Linux with Kernel version 4.4+
* R version 3.6.0+
* bowtie 1.2.2+

### Dependencies

This packages requires the ungapped short sequence aligner bowtie for microarray tag re-alignment. You can either install bowtie manually from sourcefourge:

http://bowtie-bio.sourceforge.net/index.shtml

Or install it via bioconda:

https://bioconda.github.io/recipes/bowtie/README.html

### Mira package installation

```R
install.packages("devtools")
library(devtools)
install_github("alxbetz/mira")
```

## Usage

These instructions require familiarity with the GEO database structure. If you are unfamiliar, please check their overview:

https://www.ncbi.nlm.nih.gov/geo/info/overview.html

### Preparation

1. First, we need to set up the folder structure for the pipeline:
```R
pipeline_path = '/path/to/pipeline'
mira::setup_folders(pipeline_path)
```
This will set up folders for the input, output and the reference transcriptomes.

For each GEO series that we want to analyze, we need to assign the sample filenames to the corresponding experimental groups. Therefore, create a folder for each series in '/path/to/pipeline/input' and then create a tab-separated file  '/path/to/pipeline/input/GSEXXXXX/annotation.tsv' that has 2 columns:

group | file
-----|-----
control | a1.txt.gz
control | a2.txt.gz
control | a3.txt.gz
treatment1 | b1.txt.gz
treatment1 | b2.txt.gz
treatment1 | b3.txt.gz
treatment2 | c1.txt.gz
treatment2 | c2.txt.gz
treatment2 | c3.txt.gz
  
where group is the experimental group and file must correspond to the filename in the series supplement file GSEXXXX_RAW.tar.

2. Check if the microarray platform definition file in the GEO platform entry has a column for probe ids and a column for probe sequences and note the exact column names. If there are no probe sequences in the platform definition file, you need to check the microarray manufacturer website and save the file to the folder '/path/to/pipeline/platforms/'.

3. Download the reference transcriptomes 
from Ensembl.
http://www.ensembl.org/index.html
http://plants.ensembl.org/index.html
and save them into
'/path/to/pipeline/genomicDB'

4. Create the job submission file where all input files and parameters for the run are defined and save it as '/path/to/pipeline/runfile.tsv'. This is again a tab-separated file with the following columns:


field | values | description
-----|----- | -----
database| [GEO,AE] | database of dataset, either Gene Expression Omnibus or ArrayExpress
species| string | biomaRt code for the species, usually a concatenation of the first letter of the first scientific name and the last scientific name (athaliana,celegans,drerio)
eid| string | GEO series id
platformID| string | GEO platform id. '-' if the platform file needs to be specified manually, because the GEO platform file does not have a sequence column.
platformFile| string | location of platform file with probe_id and sequence columns '/path/to/pipeline/output/platforms/platformFileName.txt'
probeIdColName| string | exact column name of the microarray probe id columns
probeSeqColName| string | exact column name of the microarray sequence
channelCount| [1,2] | indicates whether it is a dual or single channel microarray
cDNAFile| string | file location of the ensembl transcriptome
contrastString| string | A set of contrast strings that are accepted by limma's makeContrasts. All contrasts need to be separated by ','. E.g. 'treatment1-control,treatment2-control'. These must correspond the the groups you assigned in the first step.
pvalDE| double ]0-1] | adjusted p-value cutoff to determine which genes are differentially expressed
lfcDE| double ]0-inf] | log2 fold change cutoff to determine which genes are differentially expressed

Two example runfiles can be found in the `example` folder.

### Job submission
```R
jobDef = read_tsv(file.path(pipeline_path,'runfile.tsv'))
#to run a single dataset
singleResult  = mira::mira(jobDef[1,],pipedir=pipeline_path)

#to run all datasets in the job definition file
multiResult = apply(jobDef,1,mira::mira,pipedir=pipeline_path)
```

# Trouble shooting

#### Manual installation of dependencies

The R package dependencies should be installed automatically. 
If there are problems with the automatic dependency installation from the package, you can install them manually:

Bioconductor packages:

* ArrayExpress
* GEOquery
* Rbowtie
* Biostrings
* limma
* biomaRt
* topGO

CRAN packages

* tibble
* dplyr
* ggplot2
* fcros
* tidyr




