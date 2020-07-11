import os
import math
import sys
import Bio
from scipy.spatial import distance
import Spanning_Tree
from Bio.PDB import *

'''modeller_pdb_directory is the path for the output pdb files from modeller
variants_directory is the path for the files containing the variants for their corresponding PDB files
The input is the pdb file that is outputted from modeller and their corresponding variants file'''

modeller_pdb_directory = "/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/PDB_files"
variants_directory = "/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller"
if len(sys.argv) == 1:
    print("Please enter the threshold value as command line argument")
    sys.exit()
threshold = float(int(sys.argv[1])/100)
variants = []
atom_coordinates = []
variant_coordinates = []
for filename in os.listdir(modeller_pdb_directory):  # loops through all files in the modeller_pdb_directory
    if filename.endswith(".pdb"):  # checks if the file is a PDB file
        for file in os.listdir(variants_directory):  # loops through all the files in the variants_directory
            if file.__eq__(filename.rstrip(".pdb")):  # checks if the file is variants file
                file1 = open(os.path.join(variants_directory, file), "r") # opens the variants file
                for line in file1:
                    variants.append(line.rstrip("\n")) # reads the variants and stores it in the variants list
                break
        f = open(os.path.join(modeller_pdb_directory, filename), "r")  # opens the PDB file
        line = f.readline().split(" ")
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        data = parser.get_structure(line[35].rstrip("\n"), os.path.join(modeller_pdb_directory, filename))
        model = data.get_models()
        models = list(model)
        a = 0
        for m in models:
            chains = list(models[a].get_chains())
            a += 1
            b = 0
            for n in chains:
                residue = list(chains[b].get_residues())
                b += 1
                c = 0
                for o in residue:
                    atoms = list(residue[c].get_atoms())
                    c += 1
                    i = 0
                    for atom in atoms:
                        atom_coordinates.append(atoms[i].get_coord())  # gets the coordinates of all the atoms in the
                        # PDB file and stores it in atom_coordinates
                        i += 1

        for line in variants:  # loops through all the variants, converts to floating point numbers and stores the x,
            # y,z coordinates in variant_coordinates
            my_variants = line.split("\t")
            p = (float(my_variants[1]), float(my_variants[2]), float(my_variants[3]))
            variant_coordinates.append(p)

        flagg = True
        xmin = 0
        ymin = 0
        zmin = 0
        xmax = 0
        ymax = 0
        zmax = 0

        for i in atom_coordinates:  # loops through all the atom coordinates and finds the minimum x,y,z coordinate
            if flagg:
                xmin = i[0]
                ymin = i[1]
                zmin = i[2]
                flagg = False
            if i[0] < xmin:
                xmin = i[0]
            if i[1] < ymin:
                ymin = i[1]
            if i[2] < zmin:
                zmin = i[2]

        flagg = True

        for i in atom_coordinates:  # loops through all the atom coordinates and finds the maximum x,y,z coordinate
            if flagg:
                xmax = i[0]
                ymax = i[1]
                zmax = i[2]
                flagg = False
            if i[0] > xmax:
                xmax = i[0]
            if i[1] > ymax:
                ymax = i[1]
            if i[2] > zmax:
                zmax = i[2]

        abs_x = abs(xmax - xmin)
        abs_y = abs(ymax - ymin)
        abs_z = abs(zmax - zmin)

        pd = math.sqrt((abs_x * abs_x) + (abs_y * abs_y) + (abs_z * abs_z)) # calculates the dimensions for the
        # protein sequence using pythagorean theorem

        Spanning_Tree.s_tree(pd, variant_coordinates, filename, threshold) # calls the s_tree method from Spanning_Tree and creates clusters
        print("Atoms in current file clustered.")
        print()
        print()
        continue
    else:
        continue
