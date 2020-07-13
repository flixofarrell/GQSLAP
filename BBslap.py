"""========================================================
	  		BBslap V1.0
===========================================================

Overview
========
The pipeline performs the following:
   * 
   * 
   * 

Author: Felix O'Farell
Date: June 2020

   
   
Usage
=====


Configuration
-------------
The pipeline requires a pipeline.yml configuration file. This is located
within the fast-mlm directory.
Input
-----

Output
-----



Code
====
"""

import cgatcore.experiment as E
from ruffus import *
from ruffus.combinatorics import *
from cgatcore import pipeline as P
import cgatcore.iotools as iotools
import os
import sys
import re
import pandas as pd
import numpy as np
#import hail as hl
#from bokeh.io import export_png


#get parameters from YML file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
    "../pipeline.yml",
    "pipeline.yml"])


c_title = []
bfiles = []
files = PARAMS["files"]


for file in os.listdir(files):
    if file.endswith(".bed"):
       c_title.append(os.path.join(files, file))
       
for i in c_title:
    bfiles.append(i[:-4])




GRMs = os.path.abspath(
       os.path.join(PARAMS["files"] + "/GRMs"))

SPs = os.path.abspath(
       os.path.join(PARAMS["files"] + "/SPs"))

GWAS = os.path.abspath(
       os.path.join(PARAMS["files"] + "/GWAS"))

##################################################
# Generate GRM from raw plink files
##################################################
@originate(bfiles) 
@mkdir(GRMs)        
def generate_GRM(title):	

	'''
	Function to generate Genetic Relation Matrix from the
	FASTGWA-MLM analysis.
	'''	
	
	gcta = PARAMS["gcta_dir"]

	part = PARAMS["part"]

	GRMs = os.path.abspath(
      	   os.path.join(PARAMS["files"] + "/GRMs"))

	chrom = re.findall(r'chr[1-9][0-9]?$|^100$', title, re.I)
	chrom = ''.join(chrom)



	statement = '''

	scripts/GRM_loop.sh %(gcta)s %(title)s %(part)s %(GRMs)s/%(chrom)s

				'''

	P.run(statement)


##################################################
# Concatenate the GRM files 
##################################################
@follows(generate_GRM)
@originate(bfiles)
@mkdir('logs')
def concat1(title):

	'''
	Function to concatenate the .grm.id/bin/N.bin 
	files for each chromosome.
	'''
	GRMs = os.path.abspath(
      	   os.path.join(PARAMS["files"] + "/GRMs"))
	
	chrom = re.findall(r'chr[1-9][0-9]?$|^100$', title, re.I)
	chrom = ''.join(chrom)

	part = PARAMS["part"]


	statement = '''

	cat %(GRMs)s/%(chrom)s.part_%(part)s_*.grm.id > %(GRMs)s/%(chrom)s.grm.id &&
	cat %(GRMs)s/%(chrom)s.part_%(part)s_*.grm.bin > %(GRMs)s/%(chrom)s.grm.bin &&
	cat %(GRMs)s/%(chrom)s.part_%(part)s_*.grm.N.bin > %(GRMs)s/%(chrom)s.grm.N.bin

				'''

	P.run(statement)


##################################################
# Create sparse GRM from merged GRMs
##################################################
@originate(bfiles)
@follows(concat1)
@mkdir(SPs)
def sparse_GRM(title):

	'''
	Function to create sparse GRM for each chromosome.
	'''

	gcta = PARAMS["gcta_dir"]

	GRMs = os.path.abspath(
      	   os.path.join(PARAMS["files"] + "/GRMs"))

	SPs = os.path.abspath(
      	   os.path.join(PARAMS["files"] + "/SPs"))

	chrom = re.findall(r'chr[1-9][0-9]?$|^100$', title, re.I)
	chrom = ''.join(chrom)


	statement = '''
	
	%(gcta)s --grm %(GRMs)s/%(chrom)s --make-bK-sparse 0.05 --out %(SPs)s/%(chrom)s_sp


				'''

	P.run(statement)


##################################################
# Run mlm-fastgwa
##################################################
@follows(sparse_GRM)
@originate(bfiles)
@mkdir(GWAS)
def mlm_gwa(title):

	'''
	Function to run fastGWA on the plink files and GRM.
	'''	
	gcta = PARAMS["gcta_dir"]

	plinks = os.path.abspath(
      	   os.path.join(PARAMS["files"]))

	GRMs = os.path.abspath(
      	   os.path.join(PARAMS["files"] + "/GRMs"))

	GWAS = os.path.abspath(
      	   os.path.join(PARAMS["files"] + "/GWAS"))	

	SPs = os.path.abspath(
      	   os.path.join(PARAMS["files"] + "/SPs"))

	pheno = PARAMS["pheno"]

	
	chrom = re.findall(r'chr[1-9][0-9]?$|^100$', title, re.I)
	chrom = ''.join(chrom)


	statement = ''' 

	%(gcta)s --bfile %(plinks)s/%(chrom)s
	--grm-sparse %(SPs)s/%(chrom)s_sp --fastGWA-mlm 
	--pheno %(pheno)s  
	--out %(GWAS)s/%(chrom)s

				'''

	P.run(statement)


##################################################
# Concat the individual fastGWA's
##################################################

@follows(mlm_gwa)
@mkdir('masterGWA') 
def concat2():

	'''
	Function to append the previous .fastGWA summary
	file to the next. This outputs a masterfastGWA 
	summary file which contains SNPs across entire
	genome.
	'''	

	GWAS = os.path.abspath(
      	   os.path.join(PARAMS["files"] + "/GWAS"))	

	statement = '''


	cat %(GWAS)s/*.fastGWA | grep -v -e CHR -e Inverse -e Egger > %(GWAS)s/master.tsv ;
    echo -e "CHR\\tSNP\\tPOS\\tA1\\tA2\\tN\\tAF1\\tBETA\\tSE\\tP" | cat - %(GWAS)s/master.tsv > %(GWAS)s/master2.tsv ;
    mv -f %(GWAS)s/master2.tsv %(GWAS)s/master.tsv ;


				'''

				
	P.run(statement)



##################################################
# Clean up
##################################################
@follows(concat2)
@mkdir('logs')
def clean_up():


	'''
	Function to move log files into a log directory 
	'''

	statement = '''


	mv *.sh *.times *.log logs


				'''
	P.run(statement)


##################################################
# Dummy function 
##################################################

@follows(clean_up)
def full():
	
    pass

##################################################
# Pass aguments 
##################################################


#Pass arguments
@follows(full)
def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))

#end
