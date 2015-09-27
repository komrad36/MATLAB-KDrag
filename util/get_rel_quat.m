% Given two quaternions defined relative to the same reference frame, get
% quaternion relating cur_quat to des_quat. Can be used
% to calculate the error quaternion based on the
% current quaternion and desired quaternion
function rel_quat = get_rel_quat(cur_quat, desQuat)
    rel_quat = quat_mult(desQuat, quat_inv(cur_quat));
end