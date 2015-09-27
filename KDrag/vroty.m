function ret = vroty(v, angle_rad)
    Rx = [ cos(angle_rad),       0,          sin(angle_rad) ;
           0,                1,                 0   ;
           -sin(angle_rad),      0,         cos(angle_rad) ];
%     if ([1, 3] == size(v))
        ret = (Rx * v')';
%     else
%         ret = (Rx * v)';
%     end %if
end %function