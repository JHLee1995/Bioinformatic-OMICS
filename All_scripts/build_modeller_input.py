#! /usr/bin/env python3

import Bio.SeqIO
import os
import argparse

def get_protein(input_file):
    protein = {}

    try:
        #parse the fasta file
        for record in Bio.SeqIO.parse(input_file, "fasta"):
            protein[record.id] = record.seq

    except IOError:
        print("IO Error with genome file")

    print("protein made")
    return protein

def get_seq(protein_dict):

    path = os.getcwd()
    output = ""

    #loop though and write the proteins in the correct modeller input
    for key, value in protein_dict.items():
        seq = protein_dict[key] + "*"
        out =">P1;{1}{0}{2}{3}{4}{0}{5}{0}".format("\n", key, "sequence:", key, ":::::::0.00: 0.00", seq)
        output += out

        #write the data to a file
        try:
            with open("{0}/{1}_modller_input.ali".format(path, species), "w") as fasta:
                reads = output
                fasta.write(reads)

        except IOError:
            print("IO Error with making new fasta")

parser = argparse.ArgumentParser()
parser.add_argument("--species", help = "name of species")
parser.add_argument("--fasta", help = "path to fasta file")
args = parser.parse_args()

species = args.species
get_seq(get_protein(args.fasta))
