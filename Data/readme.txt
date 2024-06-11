
Contents

* Folder filtered_motorway_coordinates contains the original coordinates on each motorway network  obtained from the corresponding osm file, where the osm file is downloaded from Geofabrik (https://download.geofabrik.de). The file also contains reduced coordinates for the calculation of distances and for the map in Fig. 1 in the paper. 

* Folder natworkdistance contains the network distances between sections on each motorway network. 

* Folder NRW_number contains the different numbers of coordinates on the NRW motorway network.

* Folder NRW_shuffle contains 2000 coordinates shuffled on the NRW motorway network.

* Folder NRW_region_info contains region information obtained by program NRW_motorway_network_to_NRW_region_network.ipynb

* The source data on population density is from Statistische Ämter des Bundes und der Länder, Germany (https://regionalatlas.statistikportal.de)


Raw data processing
1. Download data *.osm.bz2 from Geofabrik (https://download.geofabrik.de) 
2. Uncompress *.osm.bz2 file as *.osm file
3. Extract the motorway network from *.osm file and save it as *_filter_by_motorway.osm file with osmosis command line in Terminal, e.g.,
osmosis  --read-xml region-latest.osm --tf accept-ways highway=motorway,motorway_link --tf reject-relations --used-node --write-xml region_filter_by_motorway.osm

4. Use python program coordinates.ipynb to extract coordinates (latitude and longitude) and save them to *_motorway.csv


