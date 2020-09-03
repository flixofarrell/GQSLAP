![Alt text](GQlogo/GQlogo.png?raw=true "Title")

[![Build Status](https://travis-ci.com/flixofarrell/GQSLAP.svg?branch=master)](https://travis-ci.com/flixofarrell/GQSLAP)


# GQSLAP
GWAS/QTL SimpLe Automated Pipeline

GQSLAP is an automated pipeline designed to process UKBB plink files on HPCs. The pipeline is built ontop of CGAT-core and uses GCTA's fastGWA linear mixed model GWAS algorithm.

## Overview

GQSLAP performs the following:
   * Searches for plink files in user designated directories
   * Creates a GRM by parts to limit the amount of memory required 
   * Concatenates the GRM parts into master GRMs for each chromosome
   * Runs GCTA's lmm GWAS algorithm on the plink files and GRMs
   * Concatenates the fastGWA output into a single master summary file
   
   
Author: Felix O'Farrell
Date: June 2020


## Installation

Currently required to clone this github directory and download GCTA's latest release and place it inside the GQSLAP dir.

## Dependencies 

CGAT-core is required and can be downloaded from anaconda

```bash
conda install -c bioconda cgatcore
conda install -c bioconda/label/cf201901 cgatcore
```
fastGWA is a feature in the GCTA package available [here](https://cnsgenomics.com/software/gcta/#Download)

## Usage

Go into the GQSLAP directory and open the pipeline.yml file. Alter the yml to configure with your system

```yml
general:
    author_name: Felix O'Farrell
    project_name: GQSLAP
    licence: MIT
    version: 1.0

# Pipeline specific options:

#dir of UKBB files
files: /path/to/your/UKBB/plink/files

#title of GRM directory (made by pipeline)
GRM_dir: GRMs

#dir of GCTA64 command 
gcta_dir: path/to/gcta/gcta64

#pheno file dir
pheno: /path/to/your/UKBB/pheno/file

#how many parts to split GRM into (save memory) 
#gcta authors advise 250 
part: 3

#will add param for covariates soon...
```

and run

```python
python GQSLAP.py make full 
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
