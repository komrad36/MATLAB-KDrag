function ret = vrotqi(vec, quat)
    % LEFT-HANDED ROTATION!
    % (like matlab toolbox)
    % BE CAREFUL!
    ret = qmul(qmul(qconj(quat), [0, vec]), quat);
    ret = ret(2:4);
    % btw there is a faster algo but this is clearer for now
end %function