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
within the BBSlap directory.
Input
-----

Output
-----



Code
====
"""

import time
from pyfiglet import Figlet
import cgatcore.experiment as E
from ruffus import *
from ruffus.combinatorics import *
from cgatcore import pipeline as P
import cgatcore.iotools as iotools
import os
import sys
import re



#get parameters from YML file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
    "../pipeline.yml",
    "pipeline.yml"])



##################################################
# Opening GQSLAP message
##################################################
if PARAMS["multiple_chr"] == 1:
  mc = 'Yes'

elif PARAMS["multiple_chr"] == 0:
  mc = 'No'

x = ' '
t = Figlet(font='doh',width = 1000)
print (t.renderText('GQSLAP'))
#print (10*x)
a = Figlet(font='slant')
print (a.renderText("Parameters"))
print(f"{'='*100}\nSingle files dir" + 3*x + "--->" + 3*x + PARAMS["files"])
print(f"{'='*100}\nMultiple Chr" + 3*x + "--->" + 3*x + mc)
print(f"{'='*100}\nPrefix" + 3*x + "--->" + 3*x + PARAMS["prefix"])
print(f"{'='*100}\nPhenotype file dir" + 3*x + "--->" + 3*x + PARAMS["pheno"])
print(f"{'='*100}\nCovariate file dir" + 3*x + "--->" + 3*x + PARAMS["cov"])
print(f"{'='*100}\nQcovariate file dir" + 3*x + "--->" + 3*x + PARAMS["qcov"])
print(f"{'='*100}\nRemove file dir" + 3*x + "--->" + 3*x + PARAMS["rmv"])
print(f"{'='*100}\nGRM construction parts" + 3*x + "--->" + 3*x + str(PARAMS["part"]))
print(f"{'='*100}\nOut dir" + 3*x + "--->" + 3*x + PARAMS["out_dir"]+"\n")

time.sleep(3)

##################################################
# Workflow related functions
##################################################
def multiple_chr_findr(files):
    
    '''
  Function to grab input strings for plink files where 
  each chromosome has a set of .fam .bim & .bed.
    '''    

    c_title = []
    bfiles = []
    
    for file in os.listdir(files):
        if file.endswith(".bed"):
           c_title.append(os.path.join(files, file))
           
    for i in c_title:
        bfiles.append(i[:-4])
        
    return (bfiles)

def single_file_makr(prefix):
	'''
	Function to make a list of one string to ensure
	each ruffus function iterates through once and 
	not once per chromosome. 
	'''
	bfiles = [prefix]
	return (bfiles)

##################################################
# Start
##################################################

GRMs = os.path.abspath(
      	os.path.join(PARAMS["out_dir"] + "/GRMs"))

SPs = os.path.abspath(
      	os.path.join(PARAMS["out_dir"] + "/SPs"))

GWAS = os.path.abspath(
      	os.path.join(PARAMS["out_dir"] + "/GWAS"))

plots = os.path.abspath(
        os.path.join(GWAS + "/Plots"))

prefix = PARAMS["prefix"]
files = PARAMS["files"]

if PARAMS["multiple_chr"] == 1:
  bfiles = multiple_chr_findr(files) 

elif PARAMS["multiple_chr"] == 0:
	bfiles = single_file_makr(prefix)
print (bfiles)

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

  if PARAMS["multiple_chr"] == 1:


    gcta = PARAMS["gcta_dir"]

    part = PARAMS["part"]

    rmv = PARAMS["rmv"]

    GRMs = os.path.abspath(
         os.path.join(PARAMS["out_dir"] + "/GRMs"))

    chrom = re.findall(r'chr[1-9][0-9]?$|^100$', title, re.I)
    chrom = ''.join(chrom)

    statement = '''
        ./scripts/GRM_loop.sh %(gcta)s %(title)s %(part)s %(GRMs)s/%(chrom)s %(rmv)s

              '''

  elif PARAMS["multiple_chr"] == 0:

    gcta = PARAMS["gcta_dir"]

    part = PARAMS["part"]

    rmv = PARAMS["rmv"]

    files = PARAMS["files"]

    prefix = PARAMS["prefix"]

    GRMs = os.path.abspath(
          os.path.join(PARAMS["out_dir"] + "/GRMs"))

    statement = '''

        ./scripts/GRM_loop.sh %(gcta)s %(files)s/%(prefix)s %(part)s %(GRMs)s/%(prefix)s %(rmv)s

              '''
  P.run(statement)

##################################################
# Concatenate the GRM files 
##################################################


@follows(generate_GRM)
@mkdir('logs')
@originate(bfiles) 
def concat1(title):

  '''
  Function to generate Genetic Relation Matrix from the
  FASTGWA-MLM analysis.
  '''

  if PARAMS["multiple_chr"] == 1:

    GRMs = os.path.abspath(
           os.path.join(PARAMS["out_dir"] + "/GRMs"))
  
    chrom = re.findall(r'chr[1-9][0-9]?$|^100$', title, re.I)
    chrom = ''.join(chrom)

    part = PARAMS["part"]


    statement = '''
        
      cat %(GRMs)s/%(chrom)s.part_%(part)s_*.grm.id > %(GRMs)s/%(chrom)s.grm.id &&
      cat %(GRMs)s/%(chrom)s.part_%(part)s_*.grm.bin > %(GRMs)s/%(chrom)s.grm.bin &&
      cat %(GRMs)s/%(chrom)s.part_%(part)s_*.grm.N.bin > %(GRMs)s/%(chrom)s.grm.N.bin

              '''

  elif PARAMS["multiple_chr"] == 0:


    GRMs = os.path.abspath(
             os.path.join(PARAMS["out_dir"] + "/GRMs"))
    
    part = PARAMS["part"]

    prefix = PARAMS["prefix"]


    statement = '''

    cat %(GRMs)s/%(prefix)s.part_%(part)s_*.grm.id > %(GRMs)s/%(prefix)s.grm.id &&
    cat %(GRMs)s/%(prefix)s.part_%(part)s_*.grm.bin > %(GRMs)s/%(prefix)s.grm.bin &&
    cat %(GRMs)s/%(prefix)s.part_%(part)s_*.grm.N.bin > %(GRMs)s/%(prefix)s.grm.N.bin      

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

  if PARAMS["multiple_chr"] == 1:

    gcta = PARAMS["gcta_dir"]

    GRMs = os.path.abspath(
             os.path.join(PARAMS["out_dir"] + "/GRMs"))

    SPs = os.path.abspath(
             os.path.join(PARAMS["out_dir"] + "/SPs"))

    chrom = re.findall(r'chr[1-9][0-9]?$|^100$', title, re.I)
    chrom = ''.join(chrom)

    statement = '''
        
    %(gcta)s --grm %(GRMs)s/%(chrom)s --make-bK-sparse 0.05 --out %(SPs)s/%(chrom)s_sp


              '''

  elif PARAMS["multiple_chr"] == 0:

    gcta = PARAMS["gcta_dir"]

    GRMs = os.path.abspath(
             os.path.join(PARAMS["out_dir"] + "/GRMs"))

    SPs = os.path.abspath(
             os.path.join(PARAMS["out_dir"] + "/SPs"))


    statement = '''

    %(gcta)s --grm %(GRMs)s/%(prefix)s --make-bK-sparse 0.05 --out %(SPs)s/%(prefix)s_sp
 

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

  if PARAMS["multiple_chr"] == 1:

    gcta = PARAMS["gcta_dir"]

    plinks = os.path.abspath(
             os.path.join(PARAMS["files"]))

    GRMs = os.path.abspath(
             os.path.join(PARAMS["out_dir"] + "/GRMs"))

    GWAS = os.path.abspath(
             os.path.join(PARAMS["out_dir"] + "/GWAS")) 

    SPs = os.path.abspath(
             os.path.join(PARAMS["out_dir"] + "/SPs"))

    pheno = PARAMS["pheno"]

    cov = PARAMS["cov"]

    qcov = PARAMS["qcov"]

    
    chrom = re.findall(r'chr[1-9][0-9]?$|^100$', title, re.I)
    chrom = ''.join(chrom)


    statement = '''
    
    %(gcta)s --bfile %(plinks)s/%(chrom)s
    --grm-sparse %(SPs)s/%(chrom)s_sp --fastGWA-mlm 
    --pheno %(pheno)s
    --out %(GWAS)s/%(chrom)s

              '''

  
  elif PARAMS["multiple_chr"] == 0:

    gcta = PARAMS["gcta_dir"]

    GRMs = os.path.abspath(
             os.path.join(PARAMS["out_dir"] + "/GRMs"))

    GWAS = os.path.abspath(
             os.path.join(PARAMS["out_dir"] + "/GWAS")) 

    SPs = os.path.abspath(
             os.path.join(PARAMS["out_dir"] + "/SPs"))

    pheno = PARAMS["pheno"]

    cov = PARAMS["cov"]

    qcov = PARAMS["qcov"]

    statement = '''
    %(gcta)s --bfile %(files)s/%(prefix)s
    --grm-sparse %(SPs)s/%(prefix)s_sp --fastGWA-mlm 
    --pheno %(pheno)s 
    --out %(GWAS)s/%(prefix)s
 

              '''

  P.run(statement)


##################################################
# Concat the individual fastGWA's
##################################################
@active_if(PARAMS["multiple_chr"] == 1)
@follows(mlm_gwa)
def concat2():

  '''
  Function to append the previous .fastGWA summary
  file to the next. This outputs a masterfastGWA 
  summary file which contains SNPs across entire
  genome.
  ''' 

  GWAS = os.path.abspath(
           os.path.join(PARAMS["out_dir"] + "/GWAS")) 

  statement = '''


  cat %(GWAS)s/*.fastGWA | grep -v -e CHR -e Inverse -e Egger > %(GWAS)s/master.tsv ;
    echo -e "CHR\\tSNP\\tPOS\\tA1\\tA2\\tN\\tAF1\\tBETA\\tSE\\tP" | cat - %(GWAS)s/master.tsv > %(GWAS)s/master2.tsv ;
    mv -f %(GWAS)s/master2.tsv %(GWAS)s/%(prefix)s.fastGWA ;


        '''
       
  P.run(statement)


##################################################
# Wrangle fastGWA 
##################################################
@follows(concat2)
def wrnglfastGWA():

  '''
  Function to wrangle the fastGWA output into 
  qqman format. 
  '''
    
  GWAS = os.path.abspath(
           os.path.join(PARAMS["out_dir"] + "/GWAS"))


  statement = '''

    sed -i -e '1s/POS/BP/' %(GWAS)s/%(prefix)s.fastGWA

              '''
  P.run(statement)

##################################################
# Get Man + QQ plots
##################################################
@follows(wrnglfastGWA)
@mkdir(plots) 
def plots():

  '''
  Function to create QQ and Manhattan plots
  in ggplot2.
  '''
    
  GWAS = os.path.abspath(
           os.path.join(PARAMS["out_dir"] + "/GWAS"))

  plots = os.path.abspath(
           os.path.join(GWAS + "/Plots"))


  statement = '''

   Rscript ./scripts/GQ_plots.R %(GWAS)s/%(prefix)s.fastGWA %(plots)s/ %(prefix)s

              '''
  P.run(statement)

##################################################
# Clean up
##################################################
@follows(plots)
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


##################################################
# End
##################################################












