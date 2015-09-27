% This function rotates a vector by a quaternion
function new_vect = v_rot_q(vect, quat )
    vect_quat = [0; vect];
    new_vect_quat = quat_mult(quat_mult(quat, vect_quat), quat_inv(quat));
    new_vect = new_vect_quat(2:4);
end

