#! /usr/bin/env python3

import Bio.SeqIO
import csv

def count_seq(input_file):

    protein = {}

    try:
        for record in Bio.SeqIO.parse(input_file , "fasta"):
            protein[record.id] = len(record.seq)

    except IOError:
        print("IO error")

    #print(protein)

    seq_file = csv.writer(open("proteins_PEMA.csv", "w"))
    for key, value in protein.items():
        #listToStr = ' '.join(map(str, value)) 
        seq_file.writerow([key, value])


    key_max = max(protein.keys(), key=(lambda k: protein[k]))
    key_min = min(protein.keys(), key=(lambda k: protein[k]))

    print('Maximum Value: ',protein[key_max], )
    print('Minimum Value: ',protein[key_min])

count_seq("/home/unhTW/share/mcbs913_2020/omics/files/protein_files/Peromyscus_maniculatus__GCA_003704035.1_HU_Pman_2.1_genomic.fna_v3.functional.proteins.fasta")
