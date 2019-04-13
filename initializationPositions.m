function [Positionsxyz, distances, disVectors, L] = initializationPositions(N,density,tau)
% Position Initialization Program
% === PURPOSE === %
% This program produces the initial (x,y,z) coordinates
%   for each N particles.
% =============== %
%
% === INPUTS === %
% N = number of particles
% density = user defined density
% tau = user defined radius of excluded volume around each particle
% =============== %
%
% === OUTPUTS === %
% Positionsxyz = A N x 3 matrix where coordinates for each 
%   particle will be stored
% =============== %
%
% === IMPORTANT NOTES === %
% When initializing a system, it is important to not put
%   any two particles too close to each other. We want to
%   initially place each particle such that each particle 
%   lives in its own sphere of radius tau. This condition
%   places a restriction on the possible densities you can
%   simulate. The maximum density that can be achieved when
%   the number of particles is a perfect cube is 1/(tau)^3.
%   However, we are not in a perfect cube; we have periodic 
%   boundary conditions. This places a new restriction on the
%   the density. Say we have a total N particles, and n is its
%   perfect cube (n^3 = N). Then the density will be 
%   d = N /(((n+1)*tau)^3. So in the limit as n is large, the 
%   density approaches 1/tau^3. But small N, N that we can simulate,
%   there is an extra factor of tau to the volume. 
%
%   What happens when N is not a perfect cube? So the question is:
%   What is the smallest cube that you can pack N balls of volume
%   ~(tau^3). It terms out that this is kind of a hard problem, from
%   some reading I was doing. Just consider two particles. Well the
%   smallest cube would be the cube with 2tau on the diagonal and the
%   two particles diagonally from each other. But what about 3? Or 4?
%   We know 8. But it was not clear to me how to construct the minimal
%   volume per N. However, I am pretty sure I can get close. 
%
%   For simplicity, find the smallest perfect cube, n, such that
%   n^3 > N, your total number of particles. For N = 2, n = 2, n^3 = 8. 
%   So the cube would just be a 2x2x2 tau cube. You can always place
%   the two particles such that they are >= tau away from each other. This 
%   is not the smallest cube, but this approach produces a clear cut way to
%   to construct the box for any N. 
%   
%   When we get the max density from this construction, we can check to
%   make sure that any user defined density is less than less max density.
%   Densities that satisfy this will have linear dimensions such that each
%   particle will have have its own excluded volume. 
% 
% ======================= %

% === PROGRAM === %

% Check if input density is less than max density %
% Have to first calculate max density from the N and tau. 

% Find the perfect cube

nCube = 2;
while(nCube^3 < N)
    nCube = nCube + 1;
end

% Calculate max density
% Minus because we are in a periodic boundary conditions, the edges
% see each other. 

maxDensity = N/(((nCube - 1)*tau)^3);

% Check if we can keep going 
if(density > maxDensity)
    fprintf('Error! : Density is larger than max density.\n');
    return;
end

% Initialize Positionsxyz matrix with zeros %
Positionsxyz=zeros(N, 1, 3);

% Calculate the side length of a cube that corresponds %
%   to the input N and density.                        %

L = (N/density)^(1/3)+tau;
%disp(L);

% step = is the length of the tick marks along each axis. %
%   We place each particle into their own cube of length  %
%   step, and this will satisfy the tau condition only if %
%   density is less than max density.                     %
stepsize = L/(nCube);
%disp(stepsize);
%disp(stepsize/2);

% We still want a random configuration given the constraint %
%   from the density. There is some "wiggle" room between   %
%   the tau we specify and the step to pick a random        %
%   position. randtau is that "wiggle" room.                %

randtau = (stepsize - tau)/2;
%disp(randtau);

% Now we need to start placing particles down. %
% Because are particles now live in a cube, 
%   all we gotta is just place them in their 
%   respective little cubes. 
%  
% If you number number each small cube, snaking all the
% way through the cube, you can generate the coordinates 
% of the middle of that small cube with just knowing what
% number it is in the snake. So all we gotta do is just
% select N random numbers 1 through nCube^3, and that gives 
% use the cube each N particle lives in. Now we have to transform
% back into [x,y,z] coordinates. randperm(N,K) does that, selects
% K random integers between 1:N without replacement. 

ww = randperm(nCube^3,N);

Positionsxyz(:,1,1) = (mod(ww-1,nCube) + 0.5)*stepsize;
Positionsxyz(:,1,2) = (mod(floor((ww-1)./nCube),nCube) + 0.5)*stepsize;
Positionsxyz(:,1,3) = (floor((ww-1)./nCube^2) + 0.5)*stepsize;

%disp(Positionsxyz);
%return;

% Randomize the new position. Pick a random displacement, R, and angles
%   theta, and Phi, using spherical coordinates to translated back to
%   x, y, z coordinates. 
    
R      = rand(N,1)*randtau;
Theta  = rand(N,1)*2*pi;
Phi    = rand(N,1)*pi;
Positionsxyz(:,1,1) = Positionsxyz(:,1,1) + R.*cos(Theta).*sin(Phi);
Positionsxyz(:,1,2) = Positionsxyz(:,1,2) + R.*sin(Theta).*sin(Phi);
Positionsxyz(:,1,3) = Positionsxyz(:,1,3) + R.*cos(Phi);
    
DiffFunction = @(x1,x2)sign(x2-x1).*min(max(x2,x1)-min(x2,x1),min(L+x1,L+x2)-max(x1,x2));
positionsPermuted = permute(Positionsxyz,[2 1 3]);
diffs = bsxfun(DiffFunction,Positionsxyz,positionsPermuted);
%disp(diffs);
disVectors = diffs;
diffSquared = diffs.^2;
distances = sum(diffSquared, 3);
distances = sqrt(distances);

shifted = abs(tau*eye(N)+distances)-tau;
%disp(shifted);
if(~isempty(find(shifted < 0)))
    fprintf("Didn't do it right\n");
    return;
end




