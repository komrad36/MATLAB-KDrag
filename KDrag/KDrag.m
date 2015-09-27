% *******************************************************************
% *   KDrag.m
% *   KDrag
% *   https://github.com/komrad36
% *
% *	9/26/2015
% *   This program is entirely my own work.
% *******************************************************************
%
% KDrag is my MATLAB implementation of my orbital aerodynamics simulation.
% MATLAB is a bit slow for this; there also exist native C++ and
% CUDA versions.
% 
% KDrag operates on a satellite or body defined as a series of polygons
% in 3-D space. When supplied with an orbital to body frame quaternion,
% density, and velocity, it approximates the drag force (if enabled)
% and torque on the body by simulating collisions and accumulating
% impulses and angular impulses per unit time. Note that force is
% returned in the orbital frame and torque is returned in the body frame.
% 
% This module can be tested standalone (see the 'caller.m' function),
% but it is intended for use by a 6 DoF orbital/attitude propagator which
% calls this module from its integrator to obtain forces and/or torques.
%
% Inputs:
%   alt   - altitude in km
%   q     - orbital to body quaternion
%   pa    - dart panel angle in degrees
%   v_sat - scalar velocity of satellite in km/s
%
% Outputs:
%   sum_F - net force in N, already converted to orbital frame
%   sum_T - net torque in N
function [sum_F, sum_T] = KDrag(alt, q, pa, v_sat)

% needed in m/s
v_sat = v_sat * 1000.0;

% rays travel in -x

PAD = 0.000001;
PITCH = 0.01;

% folding panel length
H = 0.3;

cm = [0.0; 0.0; 0.0];
F_scale_factor = 15075.514*getDensity(alt)*PITCH*PITCH;

% points of 3D panel polygons for dart
Ps =        [   [0.05 -0.05 -0.15] ;
                [0.05 0.05 -0.15] ;
                [0.05 0.05 0.15] ;
                [0.05 -0.05 0.15] ];
            
Ps(:,:,2) = [   [-0.05 -0.05 -0.15] ;
                [-0.05 0.05 -0.15] ;
                [-0.05 0.05 0.15] ;
                [-0.05 -0.05 0.15] ];
            
Ps(:,:,3) = [   [-0.05 0.05 -0.15] ;
                [0.05 0.05 -0.15] ;
                [0.05 0.05 0.15] ;
                [-0.05 0.05 0.15] ];
            
Ps(:,:,4) = [   [-0.05 -0.05 -0.15] ;
                [0.05 -0.05 -0.15] ;
                [0.05 -0.05 0.15] ;
                [-0.05 -0.05 0.15] ];
            
Ps(:,:,5) = [   [-0.05 -0.05 0.15] ;
                [0.05 -0.05 0.15] ;
                [0.05 0.05 0.15] ;
                [-0.05 0.05 0.15] ];
            
Ps(:,:,6) = [   [-0.05 -0.05 -0.15] ;
                [0.05 -0.05 -0.15] ;
                [0.05 0.05 -0.15] ;
                [-0.05 0.05 -0.15] ];
            
Ps(:,:,7) = [   [0.05 0.05 -0.15] ;
                [0.05 -0.05 -0.15] ;
                [0.05+H*sind(pa) -0.05 -0.15+H*cosd(pa)] ;
                [0.05+H*sind(pa) 0.05 -0.15+H*cosd(pa)] ];
            
Ps(:,:,8) = [   [-0.05 0.05 -0.15] ;
                [-0.05 -0.05 -0.15] ;
                [-0.05-H*sind(pa) -0.05 -0.15+H*cosd(pa)] ;
                [-0.05-H*sind(pa) 0.05 -0.15+H*cosd(pa)] ];
            
Ps(:,:,9) = [   [0.05 0.05 -0.15] ;
                [-0.05 0.05 -0.15] ;
                [-0.05 0.05+H*sind(pa) -0.15+H*cosd(pa)] ;
                [0.05 0.05+H*sind(pa) -0.15+H*cosd(pa)] ];
            
Ps(:,:,10) = [   [0.05 -0.05 -0.15] ;
                [-0.05 -0.05 -0.15] ;
                [-0.05 -0.05-H*sind(pa) -0.15+H*cosd(pa)] ;
                [0.05 -0.05-H*sind(pa) -0.15+H*cosd(pa)] ];
% ray
% RV = [-1 0 0]; % invariant

%% setup
num_poly = size(Ps, 3);
N = zeros(num_poly, 3);
precomp = zeros(1, num_poly);
min_y = zeros(1, num_poly);
max_y = zeros(1, num_poly);
min_z = zeros(1, num_poly);
max_z = zeros(1, num_poly);

%% for each time step
% precompute polygon stuff 
for i = 1:num_poly
    % rotate polygons into orbital frame
    for j = 1:size(Ps(:,:,i),1)
        Ps(j,:,i) = vrotqi(Ps(j,:,i), q);
    end %for
    
    P = Ps(:,:,i);
    N(i, :) = cross(P(2,:) - P(1,:), P(2,:) - P(3,:));
    N(i, :) = N(i, :) / norm(N(i, :));
    precomp(i) = N(i, 1)*P(1,1) + N(i, 2)*P(1,2) + N(i, 3)*P(1,3);
    min_y(i) = min(P(:,2));
    max_y(i) = max(P(:,2));
    min_z(i) = min(P(:,3));
    max_z(i) = max(P(:,3));
end %for

% get rectangular extent of sat
total_min_y = min(min_y) + PAD;
total_max_y = max(max_y) - PAD;
total_min_z = min(min_z) + PAD;
total_max_z = max(max_z) - PAD;

sum_F = [0; 0; 0];
sum_F_compensation = [0; 0; 0];
sum_T = [0; 0; 0];
sum_T_compensation = [0; 0; 0];

% for each ray
for y = [0:-PITCH:total_min_y, PITCH:PITCH:total_max_y]
    for z = [0:-PITCH:total_min_z, PITCH:PITCH:total_max_z]
        best_x = -99999999;
        for p = 1:num_poly
            % bail early if outside bounding box
            if y < min_y(p) || y > max_y(p) || z < min_z(p) || z > max_z(p)
                continue
            end %if
            
            P = Ps(:,:,p);
            num_vtx = size(P, 1);
%           otherwise, perform point-in-polygon anaylsis.
%           Looks scary but isn't too bad. Pretend coordinate system
%           is such that test point is at origin. Send the point to the right (+x)
%           and flip a flag each time it crosses a segment of the polygon. Flag will
%           end up flipped if inside polygon because it will cross an odd number of
%           segments. This works even for concave polygons.
%           For each segment of the polygon:
% 				if both ends are above or below x axis (i.e. share same sign), not a crossing
% 					otherwise, if both ends are left of the y axis, not a crossing
% 						otherwise, we do the slow bit (but usually don't have to due to the above):
% 						see if the segment intersects the x axis right of 0. if it does, crossing!
            j = num_vtx;
            odd_nodes = false;
            for i = 1:num_vtx
                odd_nodes = xor(odd_nodes, (((P(i, 3) < z && P(j, 3) >= z) || (P(j, 3) < z...
                    && P(i, 3) >= z)) && (P(i, 2) <= y || P(j, 2) <= y)) && ((P(i, 2)...
                    + (z - P(i, 3)) / (P(j, 3) - P(i, 3))*(P(j, 2) - P(i, 2)) < y)));
                j = i;
            end %for
            
            % if in polygon, compute collision location
            if odd_nodes && N(p, 1)
                x = (precomp(p) - (N(p, 2)*y + N(p, 3)*z))/N(p, 1);
                % and if it's the first thing this ray hits,
                % update best_x
                if x > best_x
                    best_x = x;
                    best_N = N(p,:);
                end %if
            end %if
        end %for
        if best_x ~= -99999999
%             fprintf('Cull-valid collision at (%g, %g, %g)\n', best_x, y, z);
            % got a valid bonk! do the math on it
            % m_sat*delta_v=2*m_particle*v_particle
            % but we just want force and will be off by a constant anyway
            % so assume 1 sec for impulse, allowing direct force convert
            % F = constant*v_norm
            % where v_particle is component NORMAL TO PANEL
            % best_N appears twice so it corrects any sign error
            v_norm = (best_N*best_N(1)*-v_sat);
            force = F_scale_factor*v_norm' - sum_F_compensation;
            
            interm_F = sum_F + force;
            sum_F_compensation = (interm_F - sum_F) - force;
            sum_F = interm_F;
            
            r = [best_x, y, z] - vrotqi(cm', q);
%             torque = cross(r, F_scale_factor*v_norm)' - sum_T_compensation;
            torque = vrotq((cross(r, force)' - sum_T_compensation)', q)';
            
            interm_T = sum_T + torque;
            sum_T_compensation = (interm_T - sum_T) - torque;
            sum_T = interm_T;
        end %for
    end %for
end %for