function vec = LatLonDist2ECEF(Azimuth,elevation,Distance)
    %This function takes a list of azimut, elevation angle and a list of distance between obs and the object.
    %pos0 is a vector define as [Lat_obs,long_obs,h_obs]
    wgs84 = wgs84Ellipsoid('kilometers');
    lat0 = 45;
    lon0 = 0;
    h0 = 0;
    [x, y, z] = aer2ecef(Azimuth, elevation, Distance, lat0, lon0, h0, wgs84);
    vec = [x(:)'; y(:)'; z(:)'];  
end
