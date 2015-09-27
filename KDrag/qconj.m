function ret = qconj(q)
% unit quats only
    ret = [q(1), -q(2), -q(3), -q(4)];
end %function