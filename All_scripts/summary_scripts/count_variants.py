#! /usr/bin/env python3

import os
import operator
import csv

def count_var():

    names={}
    dir = "/home/unhTW/share/mcbs913_2020/omics/scripts/vcf/"

    for file in os.listdir(dir):
        ortho = file.split(".")
        with open(dir+file, "r") as vcf:
            names[ortho[0]]=sum(1 for line in vcf if line.strip() and not line.startswith('#'))


    #with open("/home/unhTW/share/mcbs913_2020/omics/files/protein_files/Results_Apr15/Orthogroups.csv", "r") as og:
        #for line in og:
            #column = line.split("\t")
            #if column[0] in names:
                #names[column[0]].append(column[1:])

    #print(names)

    var_file = csv.writer(open("var.txt", "w"))
    for key, value in names.items():
        #listToStr = ' '.join(map(str, value)) 
        var_file.writerow([key, value])


    var_file.close()
    
    key_max = max(names.keys(), key=(lambda k: names[k]))
    key_min = min(names.keys(), key=(lambda k: names[k]))

    print('Maximum Value: ',names[key_max], )
    print('Minimum Value: ',names[key_min])

count_var()
