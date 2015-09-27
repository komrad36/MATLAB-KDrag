% Calculates a quaternion based on a rotation vector and rotation angle in
% degrees
function ret = makequat(direction, angle)
    halfangle = angle*pi/360;
    ret = [cos(halfangle); direction./norm(direction).*sin(halfangle)];
end %function