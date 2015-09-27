% Quaternion derivative based on a position 
% quaternion and an angular velocity vector
function quat_der = get_quat_der(quat_ang_vel)
    quat_der = 0.5*quat_mult([0; quat_ang_vel(5:7)], quat_ang_vel(1:4));
end