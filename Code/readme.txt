Author: Shanshan Wang (shanshan.wang@uni-due.de)
Reference: 
Shanshan Wang, Henrik M. Bette, Michael Schreckenberg, and Thomas Guhr. How much longer do you have to drive than the crow has to fly? npj Complex. 1,22 (2024), DOI:[10.1038/s44260-024-00023-x](https://doi.org/10.1038/s44260-024-00023-x), Preprint: arXiv:2406.06490 (2024).

Program list

To extract coordinates (latitude and longitude) from osm file,
run python/coordinates.ipynb

To calculate network distance for different countries and regions,
run python/calc_distance.py

To calculate network distance with 2000 shuffled data points for North Rhine-Westphalia (Nordrhein-Westfalen, or NRW),
run python/calc_distance_NRW_shuffle.py

To calculate network distance with different numbers of shuffled data points for North Rhine-Westphalia (Nordrhein-Westfalen, or NRW),
run python/calc_distance_NRW_shuffle_number.py

To find scaling factors and draw distributions of distances (Fig. 2 in the paper) for 9 regions (including China, US, Germany, France, Spain, Great Britain, California, Ontario, North Rhine-Westphalia),
run matlab/distance_countries.m

To find scaling factors and draw distributions of distances for 16 states in Germany (Fig. S2 in Supplementary Information),
run matlab/distance_Germany.m

To check the influence of different choices of data points for the NRW motorway network (Fig. S3 in Supplementary Information) and to find scaling factors and draw distributions of distances for the NRW motorway network (Figs. 3bc in the paper),
run matlab/distance_NRW.m

To transform NRW motorway network to NRW region network half manually and to obtain NRW region information,
run python/NRW_motorway_network_to_NRW_region_network.ipynb

To draw NRW motorway network (Fig.3a in the paper), NRW region network (Fig.3b in the paper), and fully connected region network (Fig. S1b in the Supplementary Information),
Run python/DrawMaps.ipynb

To find scaling factors and draw distributions of distances for the NRW region network (Figs. 3ef in the paper),
run matlab/distance_NRW_region_net.m

To construct fully random network with random nodes and random connections (Figs. 4a-c in the paper, and Figs. S4 and S5 in Supplementary Information)
run matlab/fully_rand_net.m

To construct random grid networks (Figs. 4d-f in the paper, and Figs. S6 and S7 in Supplementary Information)
run matlab/fully_rand_grid_net.m

To construct random geometric networks (Figs. 4g-i in the paper, and Figs. S8 and S9 in Supplementary Information)
run matlab/rand_geometric_graph.m

To construct random k-neighbour networks (Figs. S10 and S11 in Supplementary Information)
run matlab/k_neighbor_graph.m

To construct partially and fully random motorway networks in NRW with different connection fractions (Fig.5 in the paper and Figs. S12-S15 in the Supplementary Information),
run matlab/part_rand_net_with_fraction.m

To draw the distributions of geodetic distances for locations uniformly distributed in a disk, and the empirical distribution of geodetic distances for China, France, Germany and North Rhine-Westphalia, fitted with analytic equation (S4) in the Supplementary Information (Figs. S16-S17),
run matlab/understand_distance_distribution.m



