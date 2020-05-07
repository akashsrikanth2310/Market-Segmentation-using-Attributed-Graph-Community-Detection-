# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 18:15:29 2020

@author: Vidhya
"""
import pandas as pd
import numpy as np
import sys
import math
from igraph import Graph, VertexClustering
from scipy.spatial import distance



#return the attributes as a list associated with a vertex in the graph

def findAttribute(vertex):
    listreturner = []
    for iter in vertex.attributes().values(): #retunrs attribute values of each vertex [Reference link : https://igraph.org/python/doc/igraph.GraphBase-class.html]
        listreturner.append(iter)
    #print(listreturner)    
    return listreturner


# function to calculate modularity of graph 
def calculateModularity(graph,listofcommunity,iterator,secondVal,flag,storing):
    if(flag == 0): #flag checks for old or new modullarity calculation
        listofcommunity[iterator] = storing
        return graph.modularity(listofcommunity, weights = 'weight') #[Reference link : https://igraph.org/python/doc/igraph.Graph-class.html]
    else:
        listofcommunity[iterator] = secondVal #changes the list at the community value and calculate new modularity value
        return graph.modularity(listofcommunity, weights = 'weight')



#function to calculate DeltaQNewMan algorithm 

def newman(graph,listofcommunity,iterator,secondVal):
    flag = 0
    flag1 = 1
    storing = listofcommunity[iterator]
    #subtracts new and old modularity
    valueCalculate = calculateModularity(graph,listofcommunity,iterator,secondVal,flag1,storing) - calculateModularity(graph,listofcommunity,iterator,secondVal,flag,storing)
    listofcommunity[iterator] = storing 
    return valueCalculate
    

#function to calculate the deltaQAttribute algorithm

    
def attti(alphavalue,graph,listofcommunity,iterator,secondval,simmat):
    summer = 0
    counter = 0    
    counter2 = 0
    for community in listofcommunity:
        summer = summer + simmat[iterator][counter] if community == secondval else summer #checks if community mateches and based on it adds values to summer
        counter2 = counter2 + 1 if community == secondval else counter2
        counter += 1
    finalVal = counter2 * len(set(listofcommunity))
    return summer / finalVal





# executing phase 1 algorithm 

def executingPhase(graph,alphavalue,simmat):
    #creating a list for community that will hold each community for itself to start with
    listofcommunity = []
    iterator = 0 
    length = len(graph.vs)
    while(iterator < length): #since number of vertices is 324
        listofcommunity.append(iterator) #appending each node for itself
        iterator = iterator + 1
    iterator = 0
    #print(listofcommunity)
    #print(listofcommunity)
    lengthofvertex = len(graph.vs)
    convergancecheck = False   #boolean to check if we have converged
    iterationscheck = 1      #we can only run till 15 iterations
    while(iterationscheck < 16): # while loop takes care of the convergance checking
        #print('executing first iteration')
        iterator = 0
        if not convergancecheck:     #iterations is taken care with the if loop
            while(iterator < lengthofvertex):  #since number of vertices is 324
                community_picked = listofcommunity[iterator]      #picks first community to compare
                max_dq_value = -sys.maxsize -1                    #sets the most minimum value possible
                max_community = None                              #initially sets maximum community to None
                iterator2 = 0
                while(iterator2 < len(listofcommunity)):                           #since number of vertices is 324
                    if(community_picked != listofcommunity[iterator2]):                    #skips if same vertex
                        dq = alphavalue * newman(graph,listofcommunity,iterator,listofcommunity[iterator2])                         #calculates dq
                        dq += (1-alphavalue) * attti(alphavalue,graph,listofcommunity,iterator,listofcommunity[iterator2],simmat)
                        Valmin = min(max_dq_value,dq)             #finds maxval and determines which to keep
                        if(Valmin == max_dq_value and max_dq_value != dq):
                            max_dq_value = dq
                            max_community = listofcommunity[iterator2]       #changes to max_coommunity if found
                        else:
                            max_dq_value = max_dq_value
                            max_community = max_community
                    iterator2 = iterator2 + 1    
                #print(max_community)
                convergancecheck = False if max_dq_value > 0 else True              #does not change convergance if max val found > 0
                listofcommunity[iterator] = max_community if max_dq_value > 0 else listofcommunity[iterator]   #does not change convergance if max val found <= 0
                    
                iterator = iterator + 1
        else:
            break
        iterationscheck = iterationscheck + 1
    
    return listofcommunity
   
    
# function that helps us rebase the cluster values in graph 
    
def cluster_restore(community,clustermapper,listenhance,mapped):
    count = 0
    iterator = 0
    while(iterator < 324):
        createdKeys = mapped.keys()
        check = mapped.get(community[iterator], 'sure')  #to get dictionary values [reference link: https://www.tutorialspoint.com/python/dictionary_get.htm]
        listenhance.append(mapped[community[iterator]]) if(check != 'sure') else listenhance.append(count)  #if not present the append
        if(check != 'sure'):
            vals = mapped[community[iterator]]
            clustermapper[vals].append(iterator)
            iterator += 1
        else:
            mapped[community[iterator]] = count
            listcreated = []
            listcreated.append(iterator)
            clustermapper[count] = listcreated
            count = count + 1
            iterator += 1
            
    return listenhance,clustermapper

#the contract graph based on vertices and simply functions 

def contractgraph(graph,communitulist):
    graph.contract_vertices(communitulist,combine_attrs="mean")  #reference link [https://igraph.org/python/doc/igraph.Graph-class.html]
    graph.simplify(combine_edges=sum)  #reference link [https://igraph.org/python/doc/igraph.Graph-class.html]
    return graph


#function to execute phase 2 of algorithm

def executingPhase2(graph,communitulist,alphavalue,simmat):
    clustermapper = {} #creates a new map to calculate the mapped clusters
    listenhance = [] #new list that holds all the added communities
    newmap = {} #new map to perfrom same operations
    communitylist, clustermapper1 = cluster_restore(communitulist,clustermapper,listenhance,newmap) #calls the cluster restore method
    #print(len(communitylist))
    graph = contractgraph(graph,communitylist) # call contractgraph method
    
    #Calculating similarity matrix for the new vertices created
    length = len(graph.vs) #holds length of graph vertices
    simmat1 = [[None for i in range(length)] for i in range(length)] #creates similarity matrix with the none values again
    iterator1 = 0
    while(iterator1 < len(graph.vs)): 
        iterator2 = 0
        attrlist1 = findAttribute(graph.vs[iterator1]) #returns the values of attributes associated with each vertex
    #print(attrlist1)
        while(iterator2 < len(graph.vs)):
            if simmat1[iterator1][iterator2] is None:
                attrlist2 = findAttribute(graph.vs[iterator2])
                #print(attrlist2)
                cosineDistance = cosineDistance = distancefind(attrlist1, attrlist2,len(attrlist1))  #calculates the distance find parameters
                #print(cosineDistance)
                neededSimilarityVal = cosineDistance
                simmat1[iterator1][iterator2] = neededSimilarityVal
                simmat1[iterator2][iterator1] = neededSimilarityVal   #since similarity matrix is commutative
                iterator2 = iterator2 + 1
            else:
                iterator2 = iterator2 + 1       #since commuataive we skip certain rows
        iterator1 = iterator1 + 1

    
    #print('similarity recalculated')
    #Exdecuting phase 1 algorithm again
    community_list = executingPhase(graph,alphavalue,simmat1)
    #community_list = phase1(graph,alphavalue,simmat)
    #print('phase 1 reloaded')
    return community_list, clustermapper1 #returns the chagnes community list and the cluster mapped


#function to aid calculation of similarity matrix using distance

def distancefind(a1,a2,length):
    s1 = 0
    s2 = 0
    s3 = 0
    iterator = 0  
    while(iterator < length):
        s1 = s1 + a1[iterator]**2 #[using the ** operator - [Reference link : https://treyhunner.com/2018/10/asterisks-in-python-what-they-are-and-how-to-use-them/]]
        s2 = s2 + a2[iterator]**2
        val = a1[iterator]*a2[iterator]
        s3 = s3 + val
        iterator = iterator + 1
    valcalc = s1*s2
    valdiv = np.sqrt(valcalc)
    finalval = s3/valdiv
    return finalval
    


#program executing from here


length_of_arguments = len(sys.argv) #checking length of arguments in sys
if(length_of_arguments < 2):
    print('please enter the alpha value')
    sys.exit(1)

alphavalue = sys.argv[1] #picking the first argument value i.e the alpha value
alphavalue = float(alphavalue) #converting alpha value to float
#print(alphavalue)

# reads the graph as a edgelist from the given edge list file -> convertes as vertex -> edge values
graphobtained = Graph.Read_Edgelist('data/fb_caltech_small_edgelist.txt', directed=False) #[Reference Link : #https://igraph.org/c/doc/igraph-Foreign.html}

#reading a dataframe of the attributes for each node.
attributes = pd.read_csv('data/fb_caltech_small_attrlist.csv') #[Reference Link : https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_csv.html]



#to set all the vertex attributes

iterator = 0
while(iterator < 65): #65 since we have 65 vertices
   graphobtained.vs[attributes.keys()[iterator]] =  attributes[attributes.keys()[iterator]] #using vertex sqeuence to set all the attributes for the vertex #[Reference link : https://igraph.org/python/doc/igraph.VertexSeq-class.html]
   iterator = iterator + 1



#Setting all the 6005 edges weights to 1
length_of_edge_sequence = len(graphobtained.es) #selects the graph as a sequence and obtains length
iterator = 0
while iterator < length_of_edge_sequence:
    graphobtained.es[iterator]['weight'] = 1 #sets all edges to weight of one [Reference link : https://igraph.org/python/doc/igraph.EdgeSeq-class.html]
    iterator = iterator + 1




#creating a similarity matrix and assign all values to None initially
simmat = [[None for i in range(324)] for i in range(324)]
#simmat = [[None]*(len(graphobtained.vs)]*len(graphobtained.vs)
#running two loops to find the similarity martix between every pair of vertices
iterator1 = 0
while(iterator1 < 324): #324 becasue we have 0-323 vertices
    iterator2 = 0
    attrlist1 = findAttribute(graphobtained.vs[iterator1]) #returns the values of attributes associated with each vertex [Reference link : https://igraph.org/python/doc/igraph.VertexSeq-class.html]
    #print(attrlist1)
    while(iterator2 < 324):
        if simmat[iterator1][iterator2] is None:
            attrlist2 = findAttribute(graphobtained.vs[iterator2]) #[Reference link : https://igraph.org/python/doc/igraph.VertexSeq-class.html]
            #print(attrlist2)
            #cosineDistance = distance.cosine(attrlist1, attrlist2) #caluculates the cosine similarity of 1d arrays [Reference Link : https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cosine.html]
            cosineDistance = distancefind(attrlist1, attrlist2,len(attrlist1))
            #print(cosineDistance)
            #neededSimilarityVal = 1 - cosineDistance
            neededSimilarityVal = cosineDistance
            simmat[iterator1][iterator2] = neededSimilarityVal
            simmat[iterator2][iterator1] = neededSimilarityVal   #since similarity matrix is commutative
            iterator2 = iterator2 + 1
        else:
            iterator2 = iterator2 + 1       #since commuataive we skip certain rows
    iterator1 = iterator1 + 1


#print(simmat)
#running phase1 algorithm
community = executingPhase(graphobtained,alphavalue,simmat)
#running phase 2
print('phase 1 executed')
listofcommunity, clustermapping = executingPhase2(graphobtained,community,alphavalue,simmat)
print('phase 2 executed')

# clusterting vertices together and printing output
output = ''
for clust in VertexClustering(graphobtained,listofcommunity):  #{Reference link : https://igraph.org/python/doc/igraph.clustering.VertexClustering-class.html}
    if not clust:
        continue
    else:
        orig = []
        iterator =0 
        while(iterator < len(clust)): #runs till end of cluster
            extendedVal =  clustermapping[clust[iterator]] 
            orig.extend(extendedVal) #[Reference link : https://www.programiz.com/python-programming/methods/list/extend]
            iterator = iterator  + 1
        output += ','.join(map(str,orig)) #[Reference link : https://www.geeksforgeeks.org/python-map-function/]
        output += '\n'
output = output[:-2]

#Opening file and writing the required output
#print(alphavalue)
if(alphavalue == 0.5):
    alphavalue = 5
else:
    alphavalue = int(alphavalue)

file = open("communities_"+str(alphavalue)+".txt", 'w+')
file.write(output)
file.close()
