% *******************************************************************
% *   aer.m
% *   KDrag
% *   https://github.com/komrad36
% *
% *	9/26/2015
% *   This program is entirely my own work.
% *******************************************************************
%
% aer is the analytical portion of KDrag. Rather than simulating
% particle collisions, it simply estimates drag using the
% conventional drag formula off of each panel of the
% satellite, assuming the entire panel is exposed. Although
% it culls backward-facing panels, this means it does NOT
% cull obscured panels, potentially leading to some error,
% depending on the geometry of the satellite. Thus this module
% is for reference, simplicity, and speed, not necessarily
% realism, although for shapes that do not have partial
% obscuration (i.e. rectangular prism) accuracy should be
% quite good.

function [force, torque] = aer(alt, q, panel_angle_deg, v)

v = 1000.0 * v;

rho = getDensity(alt);
c_d = 4.0;
k = -1/2 * rho * c_d;
panel_angle_rad = panel_angle_deg * pi / 180.0;

% Wind Direction in Orbital Frame:
w_o = [-1 0 0];

% Surface normals for table, in Body Frame:
                            % body faces:
normals_b = [   1  0  0 ;   % +x
               -1  0  0 ;   % -x
                0  1  0 ;   % +y
                0 -1  0 ;   % -y
                0  0  1 ;   % +z
                0  0 -1 ;   % -z
                
                            % panel faces:
         vroty([1  0  0], panel_angle_rad) ;   % +x of panel hinged along +x, -z
         vroty([-1 0  0], panel_angle_rad) ;   % -x of panel hinged along +x, -z
         vroty([-1 0  0], -panel_angle_rad) ;   % +x of panel hinged along -x, -z
         vroty([1  0  0], -panel_angle_rad) ;   % -x of panel hinged along -x, -z
         vrotx([0  1  0], -panel_angle_rad) ;   % +y of panel hinged along +y, -z
         vrotx([0 -1  0], -panel_angle_rad) ;   % -y of panel hinged along +y, -z
         vrotx([0 -1  0], panel_angle_rad) ;   % +y of panel hinged along -y, -z
         vrotx([0  1  0], panel_angle_rad) ];  % -y of panel hinged along -y, -z
num_surf = 14;

% dart

%center of mass offsets, in m, from center
cm = [0, 0, 0];

% in m^2
a = [.03, .03, .03, .03, .01, .01,... %body
    .03, .03, .03, .03, .03, .03, .03, .03]; %panels

r = [.05-cm(1),  -cm(2), -cm(3) ;
     -.05-cm(1),  -cm(2), -cm(3) ;
     cm(1),  .05-cm(2), -cm(3) ;
     cm(1),  -.05-cm(2), -cm(3) ;
     cm(1),  -cm(2), .15-cm(3) ;
     cm(1),  -cm(2), -.15-cm(3) ;
   
   .05-cm(1)+.15*sin(panel_angle_rad), -cm(2), -.15-cm(3)+.15*cos(panel_angle_rad) ;
   .05-cm(1)+.15*sin(panel_angle_rad), -cm(2), -.15-cm(3)+.15*cos(panel_angle_rad) ;
   -.05+cm(1)-.15*sin(panel_angle_rad), -cm(2), -.15-cm(3)+.15*cos(panel_angle_rad) ;
   -.05+cm(1)-.15*sin(panel_angle_rad), -cm(2), -.15-cm(3)+.15*cos(panel_angle_rad) ;
   
   -cm(1), .05-cm(2)+.15*sin(panel_angle_rad), -.15-cm(3)+.15*cos(panel_angle_rad) ;
   -cm(1), .05-cm(2)+.15*sin(panel_angle_rad), -.15-cm(3)+.15*cos(panel_angle_rad) ;
   -cm(1), -.05+cm(2)-.15*sin(panel_angle_rad), -.15-cm(3)+.15*cos(panel_angle_rad) ; 
   -cm(1), -.05+cm(2)-.15*sin(panel_angle_rad), -.15-cm(3)+.15*cos(panel_angle_rad) ]; 
   
% faster to only do a single conversion, of the wind
% instead of the normals:
w_b = vrotq(w_o, q);

% now project wind (w_b) onto each surface
% normals are unit so:
% proj of a onto b is (a dot b)*b

f = zeros(num_surf, 3);
t = zeros(num_surf, 3);
for i = 1:num_surf
    dprod = dot(-w_b, normals_b(i,:));
    if dprod > 0
        f(i,:) = k * a(i)*(v*dprod)^2*normals_b(i,:);
        t(i,:) = cross(r(i,:), f(i, :));
    end %if
end %for

force = vrotqi(sum(f), q)';
torque = sum(t)';

end %function