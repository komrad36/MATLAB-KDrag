% This function calculates the torques required to acheive a given satellite 
% angular acceleration. These 'required torques' can be acieved by applying a
% torque directly to the body or by accelerating a reaction wheel such that
% the wheel accelerates opposite the direction of the desired body torque
% vector
function req_torques = calc_req_torques(des_ang_accel, sat_ang_vel, ang_momentum, moi_matrix)
    req_torques = moi_matrix*des_ang_accel - cross(ang_momentum, sat_ang_vel);
end