% *******************************************************************
% *   caller.m
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
% This module is a simple test dispatcher for the MATLAB version
% of the KDrag simulation, as well as the 'aer' fast approximation
% analytical module.

alt = 500;
% q = [0 1 0 0];
% q =          [sqrt(2)/2
%                          0
%         -sqrt(2)/2
%                          0];
% q=         [0.737277336810124
%                          0
%          -0.67559020761566
%                          0];
% q = [      0.76604
%             0
%      -0.64279
%             0];
%   q=[    0.64279
%             0
%      -0.76604
%             0];
% q = [1 0 0 0];

q = [sqrt(2)/2 0 sqrt(2)/2 0];
% q = makequat([0; 0; -1], -45);

v = sqrt(398600.44189/(6371+alt));
panel_angle = 135;

clc

[faer, taer] = aer(alt, q, panel_angle, v);
[fkdrag, tkdrag] = KDrag(alt, q, panel_angle, v);
faer
taer
fkdrag
tkdrag