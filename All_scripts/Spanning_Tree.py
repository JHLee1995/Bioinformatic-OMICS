from scipy.spatial import distance
from collections import defaultdict
from networkx import *


def s_tree(pd, variant_coordinates, filename, threshold):
    list_variants = []
    distance_matrix = []
    dict = {}

    for i in variant_coordinates:  # loops through all the variants and stores it in the list
        list_variants.append(i)

    for i in range(list_variants.__len__()):  # Creates a dictionary of all variant coordinates
        dict.update({str(i): list_variants[i]})

    '''The following loop takes each pair of (x,y,z) coordinates and calculates the euclidean distance between them.
    After calculating the euclidean distance it creates a distance matrix and stores it in distance_matrix list in the 
    following format [source, destination, euclidean distance]'''
    count = 0
    for i in dict:
        countj = 0
        for j in dict:
            if dict[i] is not dict[j]:
                if (countj > count):
                    d = distance.euclidean(dict[i], dict[j])
                    distance_matrix.append([i] + [j] + [d])
            countj += 1
        count += 1

    while dict.__len__() > 1:  # loops over the dictionary of variant coordinates
        lowest = 0
        source = 0
        dest = 0
        flag = True
        for i in distance_matrix:  # loops over the distance matrix and gets the coordinate pairs with the smallest
            # euclidean distance
            if flag:
                flag = False
                lowest = i[2]
                source = i[0]
                dest = i[1]
                continue

            if i[2] < lowest:
                lowest = i[2]
                source = i[0]
                dest = i[1]

        '''The following if condition normalizes the euclidean distance and checks if it is greater than the 
        dimension threshold, if the normalized distance is greater than the threshold, it breaks the loop and outputs 
        the final number of clusters '''
        if (lowest / pd) > threshold:
            break

        x_s = dict[source][0]  # gets the source point's x-coordinate
        y_s = dict[source][1]  # gets the source point's y-coordinate
        z_s = dict[source][2]  # gets the source point's z-coordinate

        x_d = dict[dest][0]  # gets the destination point's x-coordinate
        y_d = dict[dest][1]  # gets the destination point's y-coordinate
        z_d = dict[dest][2]  # gets the destination point's z-coordinate

        x = (x_s + x_d) / 2  # calculates the centroid for x-coordinate
        y = (y_s + y_d) / 2  # calculates the centroid for y-coordinate
        z = (z_s + z_d) / 2  # calculates the centroid for z-coordinate

        dict.pop(source)  # removes the source data point from the dictionary
        dict.pop(dest)  # removes the destination data point from the dictionary

        dict.update({str(source) + str(dest): (x, y, z)})  # adds the centroid to the dictionary

        distance_matrix.clear()  # clears the entire distance matrix

        '''The following loop again takes each pair of coordinates from the dictionary, calculates the euclidean 
        distance between them and creates a distance matrix '''
        count = 0
        for i in dict:
            countj = 0
            for j in dict:
                if dict[i] is not dict[j]:
                    if (countj > count):
                        # print(str(dict[i]) + "," + str(dict[j]))
                        d = distance.euclidean(dict[i], dict[j])
                        num = d / pd
                        # print("Euclidean distance: ", d)
                        distance_matrix.append([i] + [j] + [d])
                countj += 1
            count += 1

    file = open("Cluster_Output.txt", "a")  # Creates and stores the output to the Cluster_Output.txt file
    file.write("Clustered Output for PDB File: " + filename + "\n")
    print("Clustered Output for PDB file: " + filename)
    counter = 1
    for i in dict:
        print("Cluster " + str(counter) + ": " + str(dict[i]))
        file.write("Cluster " + str(counter) + ": " + str(dict[i]) + "\n")
        counter += 1
    file.write("\n")
    file.close()
