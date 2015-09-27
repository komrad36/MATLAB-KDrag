function qinv = quat_inv( quat )
% This function calculates the quaternion inverse
% The quaternion inverse magnitude should always be one.
    qinv = [quat(1); -quat(2:4)];
end

