function dkm=geodistkm(latlon1,latlon2)
% d1km: distance in km based on Haversine formula
% (Haversine: http://en.wikipedia.org/wiki/Haversine_formula)
%
% --Inputs:
%   latlon1: latlon of origin point [lat lon]
%   latlon2: latlon of destination point [lat lon]
%
% --Outputs:
%   dkm: distance calculated by Haversine formula
%   
% --Example:
%   latlon1=[-43 172];
%   latlon2=[-44  171];
%   dkm=geodistkm(latlon1,latlon2)
%   dkm =
%           137.365669065197 (km)
%--------------------------------------------------------------------------
radius=6371;
lat1=latlon1(1)*pi/180;
lat2=latlon2(1)*pi/180;
lon1=latlon1(2)*pi/180;
lon2=latlon2(2)*pi/180;
deltaLat=lat2-lat1;
deltaLon=lon2-lon1;
a=sin((deltaLat)/2)^2 + cos(lat1)*cos(lat2) * sin(deltaLon/2)^2;
c=2*atan2(sqrt(a),sqrt(1-a));
dkm=radius*c;    %Haversine distance
end