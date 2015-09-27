function [quat] = dcm2Quat(dirCos)
    quat = zeros(4,1);
    traceDirCos = trace(dirCos); % sum of diagonal terms
    quat(1) = sqrt((traceDirCos+1)/4); % cos(theta/2)
    % correspond to sin(theta/2)*(rotation axis)
    quat(2) = (dirCos(3,2) - dirCos(2,3))/(4*quat(1));
    quat(3) = (dirCos(1,3) - dirCos(3,1))/(4*quat(1));
    quat(4) = (dirCos(2,1) - dirCos(1,2))/(4*quat(1));
end