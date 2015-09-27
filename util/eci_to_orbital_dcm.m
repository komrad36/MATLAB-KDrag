% Compute DCM from ECI such that x axis of orbital frame is aligned with
% the velocity vector, z azis is aligned with the nadir vector, and the
% y axis completes the right handed coordinate frame
function dcm = eci_to_orbital_dcm(pos_vect, vel_vect)
    % unit velocity vector in the orbital frame 
    vel_unit_orb = [1 0 0]';
    nadir_vect_orb = [0 0 1]';
    vel_unit_eci = vel_vect/norm(vel_vect);
    nadir_vect_eci = -pos_vect/norm(pos_vect);
    
    % Calculate the dcm relating the eci frame to the orbital frame using
    % the triad algorithm
    dcm = getAttitude(nadir_vect_eci, nadir_vect_orb, vel_unit_eci, vel_unit_orb);
end