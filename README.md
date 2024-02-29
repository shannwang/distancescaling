# distance scaling

## Data

### Procedure of raw data processing
1. Download data *.osm.bz2 from Geofabrik (https://download.geofabrik.de) 
2. Uncompress *.osm.bz2 file as *.osm file
3. Filter the motorway network from *.osm file and save it as *_filter_by_motorway.osm file with osmosis command line in Terminal, e.g.,\
`
osmosis  --read-xml region-latest.osm --tf accept-ways highway=motorway,motorway_link --tf reject-relations --used-node --write-xml region_filter_by_motorway.osm
`
4. Use python program coordinates.py to filter the data of coordinates (latitude and longitude) and save it as *_motorway.csv
