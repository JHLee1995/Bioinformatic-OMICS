#! /usr/bin/env python3

import subprocess
import os

def run_orthofinder():
    result = subprocess.run(" orthofinder2 -f /home/unhTW/share/mcbs913_2020/omics/files/protein_files/ -S diamond -M msa -A mafft", shell=True, stdout=subprocess.PIPE)
    output = result.stdout.decode('ascii')

def get_vars(file_path_fa):

    wd = os.getcwd()
    path = wd+"/vcf"

    try:
        #a new dir is created to hold the vcf out put files 
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s " % path)

    try:
        #for every .fa files created in orthofinder the prgram snp-sites is run to call variatns 
        for file in os.listdir(file_path_fa):
            result = subprocess.run("snp-sites -v -o {0}/{1}.vcf {2}{1}".format(path, file, file_path_fa), shell=True, stdout=subprocess.PIPE)
            #output = result.stdout.decode('ascii')

    except IOError:
        print("IO Error with Orthogroups.csv")

run_orthofinder()
get_vars("/home/unhTW/share/mcbs913_2020/omics/files/protein_files/Results_Apr15/Orthologues_Apr15/Alignments/")

