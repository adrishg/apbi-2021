#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!pip install regex
import re
#!pip install scipy
from scipy.spatial import distance
#!pip install numpy
import numpy as np
#!pip install pandas
import pandas as pd


# In[2]:


#Using regex to extract coordinates of alfa carbons of the pdb file
def extractAlfaC(pdbFile):
    
    patternATOM = re.compile("ATOM.*CA +([A-Z]{3}) +[A]{1}")
    patternCOOR_CA = re.compile("^ATOM.*CA +([A-Z]{3}) +[A]{1} *([0-9]+) +(\-?[0-9\.]+)+ +(\-?[0-9\.]+) +(\-?[0-9\.]+) +(\-?[0-9\.]+)")

    coor = []

    for i, line in enumerate(open(pdbFile)):
        #print(line)
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
    return pdbCoor


# In[16]:


#Calculate 3D distances all vs all c_αs
def calculateDistances(pdbCoor):
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
    return distances


# In[4]:


# Mirror already calculated distances to create a square matrix
def squareMatrix(distances, naa):
    for i in reversed(range(naa)):
        for j in reversed(range(1, naa)):
            if len(distances[i])< naa:
                distances[naa-i][naa-j]= distances[i][j]
    return distances


# In[5]:


#!pip install matplotlib
import matplotlib.pyplot as plt

def distanceMapCreation(distances):
    plt.imshow(distances)

    cbar = plt.colorbar();
    cbar.set_label('Distance (Å)') 
    cbar.set_ticks([10,20,30])

    plt.show()


# In[17]:


def distanceMap(pdbFile):
    pdbCoor = extractAlfaC(pdbFile)
    distances = calculateDistances(pdbCoor)
    distanceMap = distanceMapCreation(distances)
    


# In[18]:


distanceMapExample2 = distanceMap('example_file_2.pdb')

