function [coordsPast, coordsPres, newVelocity] = evolveSystem(Accelerations,coordsPres,coordsPast,dt,L)
dts = dt^2;

New = 2*coordsPres - coordsPast + Accelerations*dts;
New = mod(New, L);

newVelocity = (New - coordsPres)/dt;

coordsPast = coordsPres;
coordsPres = New;