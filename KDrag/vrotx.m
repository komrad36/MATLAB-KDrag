function ret = vrotx(v, angle_rad)
    Rx = [ 1,       0,          0 ;
           0, cos(angle_rad), -sin(angle_rad)
           0, sin(angle_rad),  cos(angle_rad) ];
%     if ([1, 3] == size(v))
        ret = (Rx * v')';
%     else
%         ret = (Rx * v)';
%     end %if
end %function