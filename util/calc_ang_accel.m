% This function calculates the angular acceleration vector of a body based
% on the body's moment of inertia matrix, angular momentum vector,
% angular momentum derivative, and angular velocity vector, and a 
% vector of applied torques
function ang_accel = calc_ang_accel(torque, ang_vel, ang_momentum, moi_matrix)
    % Calculations based on equation dH/dt = H x w + T
    ang_mom_der = cross(ang_momentum, ang_vel) + torque;
    ang_accel = (inv(moi_matrix))*ang_mom_der;
end