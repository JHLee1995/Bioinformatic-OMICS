#! /usr/bin/env python3

def count_id():

    species_id = set()
    kbid = set()
    pdb = set()

    try:
        with open("/home/unhTW/share/mcbs913_2020/omics/scripts/PEER_modller_input.fasta", "r") as species:
            for line in species:
                if line.startswith(">"):
                    ids = line.split("|")
                    species_id.add(ids[1])

    except IOError:
        print("IO error with fasta file")

    try:
        with open("/home/unhTW/share/mcbs913_2020/omics/scripts/idmapping_subset2-1.txt", "r") as idmap:
            for line in idmap:
                column = line.split("\t")
                if column[1] == "PDB":
                    pdb.add(column[0])
                else:
                    kbid.add(column[0])

    except IOError:
        print("IO error with idmapping file")


    print("number of PDBs in the PEMA gff:  ", len(species_id.intersection(pdb)))
    print("number of UniProtKB-IDin the PEMA gff: ", len(species_id.intersection(kbid)))

count_id()
