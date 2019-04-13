function [points] = polymerInitialization(N,l)
% Position Initialization Program
% === PURPOSE === %
% This program produces the initial configuration of the 
% polymer. 
% =============== %
%
% === INPUTS === %
% N = number of segments along the polymer
% l = length of each segment
% =============== %
%
% === OUTPUTS === %
% Positionsxyz = A N x 3 matrix where coordinates for each 
%   particle will be stored
% =============== %
%

% === PROGRAM === %

% Set the first point, initialP
initialP = [0,0,0];

% Initialize a Nx1x3 stack of coordinates for the N+1 points
points = zeros(N+1,1,3);
% Set the first point coordinates
points(1,1,:) = initialP;

disp(points)
% Make Nx1 vectors of random Theta and Phi angles
Theta = rand(N,1)*2*pi;
Phi   = rand(N,1)*pi;

% Take the cos and sin of these angles, these vectors will be used 
% in the for loop
cosTheta = cos(Theta);
sinTheta = sin(Theta);
cosPhi = cos(Phi);
sinPhi = sin(Phi);

% The length of each rod, l
l=1;

for i=1:N
    points(i+1,1,1) = points(i,1,1) + l*cosTheta(i)*sinPhi(i);
    points(i+1,1,2) = points(i,1,2) + l*sinTheta(i)*sinPhi(i);
    points(i+1,1,3) = points(i,1,3) + l*cosPhi(i); 
end

distfunc = @(x1,x2) sqrt((x2(1)-x1(1))^2+(x2(2)-x1(2))^2+(x2(3)-x2(3))^2);

dist = distfunc(points(1,1,:),points(N+1,1,:));

L = floor(2*dist);
disp(L)
close all
grid on
ax=gca;
ax.GridColor = [0 0.5 0.5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
set(ax,'Xlim',[-L L]);
set(ax, 'Xtick',(-L:0.5:L));
set(ax,'Ylim',[-L L]);
set(ax, 'Ytick',(-L:0.5:L));
set(ax,'Zlim',[-L L]);
set(ax, 'Ztick',(-L:0.5:L));
ax.Layer = 'bottom';
hold(ax);
view(3);

plot3(points(:,1,1),points(:,1,2),points(:,1,3),'k-','LineWidth',3)