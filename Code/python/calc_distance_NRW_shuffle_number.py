#----------------------------------------------------------------
# calculate network distance with different numbers of shuffled data points 
# for Nordrhein-Westfalen (NRW)
# author: shanshan wang
# email: shanshan.wang@uni-due.de
# Feb. 29, 2024
#----------------------------------------------------------------
import pandas as pd
import numpy as np
import osmnx as ox
import networkx as nx
import geopandas as gpd
import csv
import time

start_time=time.time()
name1='NRW'
name2='nordrhein-westfalen'
dir1='data/'
dir2='results/'

#--change the number of shuffled data points---
num=200
#----------------------------------------------

G=ox.graph.graph_from_xml(dir1+name2+'_filter_by_motorway.osm',  bidirectional=True, simplify=False, retain_all=True)

sectab=pd.read_csv(dir1+'NRW_number/'+name1+'_mw_less_shuffle_number_'+str(num)+'.csv')
n=sectab.shape[0]
nodes=[]
l=[]
length=[]
for s in range(0,n):
    nodes.append(ox.distance.nearest_nodes(G,sectab.iloc[s,:].lon, sectab.iloc[s,:].lat))
for s in range(0,n):
    length = nx.shortest_path_length(G, source=nodes[s],  weight='length',method='dijkstra') # meters by default
    l.append(length)
b=[]
for s in range(0,n):
    a=[]
    for t in range(0,n):
        try:
            a.append(round(l[s][nodes[t]]/1000,4))
        except KeyError:
            a.append("NaN")
    b.append(a)

del sectab
del nodes
del l

# save data
with open(dir2+name1+"_networkdistance_motorway_shuffle_number_"+str(num)+".csv", "w", newline="") as f:
    writer = csv.writer(f, delimiter =',')
    writer.writerows(b)    
print("---%s second---" % (time.time()-start_time))
