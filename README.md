# KDrag
Orbital and attitude propagator with B-dot and *dynamic* aerodynamic drag simulation, including torque computation for aero-stabilized bodies.

This project is in three parts:

1) MATLAB orbital and attitude propagator with rotating orbital frame, B-dot,
and calls to an aerodynamics engine.

2) The full MATLAB aerodynamics engine, KDrag.

3) A fast, low-precision reference aero engine called 'aer'.

KDrag is my MATLAB implementation of my orbital aerodynamics simulation.
MATLAB is a bit slow for this; there also exist native C++ and
CUDA versions on my github.

KDrag operates on a satellite or body defined as a series of polygons
in 3-D space. When supplied with an orbital to body frame quaternion,
density, and velocity, it approximates the drag force (if enabled)
and torque on the body by simulating collisions and accumulating
impulses and angular impulses per unit time. Note that force is
returned in the orbital frame and torque is returned in the body frame.

This module can be tested standalone (see the 'caller.m' function),
but it is intended for use by a 6 DoF orbital/attitude propagator which
calls this module from its integrator to obtain forces and/or torques.

Such an integrator is available as O-A-prop.m in this repo.

Please note that the IGRF data used for magnetorquer B-dot are not
provided by me; rather, they are from Drew Compston of MathWorks:
http://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field--igrf--model

Please see the sample output 'Demo.png'!
