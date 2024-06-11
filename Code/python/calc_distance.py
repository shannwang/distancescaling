#----------------------------------------------------------------
# calculate network distance with the coordinates and motorway networks 
# for different countries and regions
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
#--------change names,see name1_name2_for_cal_dist.txt-----------
name1='Germany'
name2='germany'
#----------------------------------------------------------------
dir1='data/'
dir2='results/'
sectab=pd.read_csv(dir1+name1+'_mw_less.csv')
G=ox.graph.graph_from_xml(dir1+name2+'_filter_by_motorway.osm',  bidirectional=True, simplify=False, retain_all=True)

n=sectab.shape[0]
nodes=[]
l=[]
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
            a.append(round(l[s][nodes[t]]/1000,4)) # change distance unit from meter to kilometer
        except KeyError:
            a.append("NaN")
    b.append(a)



# save data
with open(dir2+name1+"_networkdistance_motorway.csv", "w", newline="") as f:
    writer = csv.writer(f, delimiter =',')
    writer.writerows(b)
    
print("---%s second---" % (time.time()-start_time))
