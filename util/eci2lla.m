% Compute the lat, long, and alt of a
% point in space given a position vector in ECI
function lat_long_alt = eci2lla(eci_pos_vect, time)
    earth_radius = 6378;  % km
    
     % lattitude in degrees
    lat= atand(eci_pos_vect(3)/sqrt(eci_pos_vect(1)^2 + eci_pos_vect(2)^2))/pi*180;
    
    % negative means west of the vernal equinox and positive means east
    % angle from vernal equinox in radians
    ang_from_vern = atan2(eci_pos_vect(2), eci_pos_vect(1));
    
    % Get the greenwich hour angle
    gha = get_gha(time); % degrees
    longitude = ang_from_vern - gha; % degrees
    
    % ensure that longitude is between -180 and 180 degrees
    if longitude < -180
        longitude = longitude + 360;
    elseif longitude > 180
        longitude = longitude - 360 ;   
    end 
    
    altitude = norm(eci_pos_vect) - earth_radius;
    
    lat_long_alt = [lat, longitude, altitude];
end