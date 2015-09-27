function ret = vrotq(vec, quat)
    % RIGHT-HANDED ROTATION!
    % MATLAB uses left-handed in the aero toolbox
    % BE CAREFUL!
    ret = qmul(qmul(quat, [0, vec]), qconj(quat));
    ret = ret(2:4);
    % btw there is a faster algo but this is clearer for now
end %function