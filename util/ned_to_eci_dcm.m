% Compute DCM relating a North-East-Down 
% frame to the eci frame. Inputs are:
% lattitude, longitude, greenwich hour angle (gha)
function dcm = ned_to_eci_dcm(lat, long, gha)
    ang = long + gha; % angle between vernal equinox and point
    % If >0 the point is east of the vernal equinox
    dcm = [-sind(lat)*cosd(ang) -sind(lat)*sind(ang) cosd(lat);
            -sind(ang) cosd(ang) 0;
            -cosd(lat)*cos(ang) -cosd(lat)*sind(ang) -sind(lat)];
end