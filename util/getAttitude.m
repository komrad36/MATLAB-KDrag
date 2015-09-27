% Returns a direction cosine matrix and a quaternion relating the body
% frame to the reference frame. The dirCos matrix times a ref frame vector
% will give that vector in the body frame
function [dirCos quat] = getAttitude(A, B, C, D)
    vect1Ref = A;
    vect1Body = B;
    vect2Ref = C;
    vect2Body = D;

    vect2Ref = vect2Ref/norm(vect2Ref);
    vect1Ref = vect1Ref/norm(vect1Ref);
    vect2Body = vect2Body/norm(vect2Body);
    vect1Body = vect1Body/norm(vect1Body);

    crossMagRef = cross(vect1Ref, vect2Ref); % cross product of reference frame vectors
    crossMagRef = crossMagRef/norm(crossMagRef); % make into a unit vector
    crossMagBody = cross(vect1Body, vect2Body);
    crossMagBody = crossMagBody/norm(crossMagBody);

    vect1CrossRef = cross(vect1Ref,crossMagRef); % cross vect1 with cross product of vect1 and vect2
    vect1CrossRef = vect1CrossRef/norm(vect1CrossRef); % make unit vector
    vect1CrossBody = cross(vect1Body,crossMagBody);
    vect1CrossBody = vect1CrossBody/norm(vect1CrossBody);

    refVects = [vect1Ref crossMagRef vect1CrossRef];
    bodyVects = [vect1Body crossMagBody vect1CrossBody];
    %bodyVects = A*RefVects
    dirCos = bodyVects*refVects';

    % calculate quaternion
    quat = dcm2Quat(dirCos);
    