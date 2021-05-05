#!/usr/bin/env python

#!pip install regex
import re
#!pip install scipy
from scipy.spatial import distance
#!pip install numpy
import numpy as np
#!pip install pandas
import pandas as pd
#!pip install matplotlib
import matplotlib.pyplot as plt



def distanceMap(pdbFile):
    #Extract PDB code
    patternPDBCode = re.compile("HEADER.*([A-Z0-9]{4}) *$")
    #Using regex to extract coordinates of alfa carbons of the pdb file
    patternATOM = re.compile("ATOM.*CA +([A-Z]{3}) +[A]{1}")
    patternCOOR_CA = re.compile("^ATOM.*CA +([A-Z]{3}) +[A]{1} *([0-9]+) +(\-?[0-9\.]+)+ +(\-?[0-9\.]+) +(\-?[0-9\.]+) +(\-?[0-9\.]+)")

    coor = []

    for i, line in enumerate(open(pdbFile)):
        pdbCodeSearch = re.search(patternPDBCode, line)
        if pdbCodeSearch:
            pdbCode = pdbCodeSearch.group(1)
            print(pdbCode)
        for match in re.finditer(patternATOM, line):
            #print(line)
            coordinateSearch = re.search(patternCOOR_CA,line)
            if coordinateSearch:
                aa= coordinateSearch.group(1)
                n = coordinateSearch.group(2)
                x = float(coordinateSearch.group(3))
                y = float(coordinateSearch.group(4))
                z = float(coordinateSearch.group(5))

                coor.append([aa+' '+n , x , y , z])

    pdbCoor = pd.DataFrame( coor , columns = ['AA', 'x', 'y', 'z' ] )
    #Calculate Distances
    xyz= pdbCoor[['x','y','z']]
    distances = []
    naa = len(pdbCoor)
    for i in range(naa):
        distances.append([])
        for j in range(naa):
            dist = distance.pdist([xyz.iloc[i].values, xyz.iloc[j].values])
            distances[i].append(dist[0])
    #Square matrix
    for i in reversed(range(naa)):
        for j in reversed(range(1, naa)):
            if len(distances[i])< naa:
                distances[naa-i][naa-j]= distances[i][j]
    #PlotDistances
    plt.imshow(distances)

    cbar = plt.colorbar();
    cbar.set_label('Distance (Å)') 
    cbar.set_ticks([10,20,30])

    #Name of the file from the PDB Code
    nameDistanceMap = pdbCode+"distMap.png"

    plt.savefig( nameDistanceMap )

distanceMap('example_file_2.pdb')