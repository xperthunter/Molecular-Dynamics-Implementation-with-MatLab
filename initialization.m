function [coords, distances, closestVectors, velocity, L] = initialization(N,density,tau,T)
% Setting the simulation %
%
% === INPUTS === %
% N = number of particles
% density = density we want for the simulation
%   This has to be less than the maximum density. 
%   Program, densities.m, calculates the max density
%   for a vector of N values and a specific tau.
%   It also spits out candidate densities if you 
%   don't know what to chose. 
% tau = radius of sphere of excluded volume around each
%   particle. This is just for initialization. This is not
%   to say that any two particles, through the course of the
%   simulation will never be within tau of each other.
%   It is just important to make sure that initially, the particles
%   are not too close. 
% T = temperature of the system. 
%
% === OUTPUTS === %
% coords = [x,y,z] coordinates for each particle. coords is a (N,1,3)
%   stack. Basically a stack of three (N,1) vectors.
% distances = (N,N) matrix containing the distances between each pair of
%   particles
% closestVectors = [x,y,z] componets of the vector representing the closest
%   distance between every pair of points. This is a matrix that is
%   (N,N,3).
% velocity = [x,y,z] components of initial velocity vectors for each 
%   particle. Again, velocity is a (N,1,3) stack. 
%
% === PROGRAM === %

% From the input temperatue, lets compute the total kinetic energy. 
% kT = 2*KE/(3N-3). Its -3 because momentum is set to zero. 

initialKE = (3*N-3)*0.5*T;

% Get initial positions, distances, the vectors that are closest
%   between pairs of points, and length of box. 

[coords, distances, closestVectors, L] = initializationPositions(N,density,tau);

% Velocity initialization
% Pick uniform number between [-0.5,0,5] for each component of each
%   particles velocity. 

velocity = rand(N,1,3) - 0.5;

% Set momentum to zero. 
% Average velocity
vsum = sum(velocity)/N;
% Shift velocities so that the momentum center of mass is at the center. 
velocity = velocity - vsum;

% Function to get magnitude of velocity squared. 
magVelocity = @(vx, vy, vz) vx.^2 + vy.^2 + vz.^2;

% Sum of velocity squares
vSquare = sum(magVelocity(velocity(:,1,1),velocity(:,1,2),velocity(:,1,3)));

% Initial kinetic energy
KEinit = 0.5*vSquare;

% Re-scale factor to get the desired temperature
vRescale = sqrt(initialKE/KEinit);
% Scale the velocities by re-scale
velocity = velocity*vRescale; 

% If you want to check and make sure this is working, uncomment these lines
%vCheck = sum(magVelocity(velocity(:,1,1),velocity(:,1,2),velocity(:,1,3)));
%vCheck = vCheck/(3*N-3);
%disp(vCheck);

