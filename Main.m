function [allPositions,L] = Main(N, density, T, dt, TotalTime, fid, fid2)

tau = 2^(1/6);
% Parameter for Andersen Thermostat
%nu = 0.5;
%nudt = nu*dt;
initialKE = (3*N-3)*0.5*T;

[coordsPres, distances, closestVectors, velocity, L] = initialization(N,density,tau,T);


% coordsPres is a Nx1x3 stack containing the x, y, z coordinates of each
%   particle. 
% distances is a NxN matrix of distances. 
% closestVectors are the components of the vectors of the closest distance
%   between each particle. NxNx3 
% L is the length of the box.
% velocity is a Nx1x3 stack containing the components of the velocity
%   vectors.

coordsPast = coordsPres - dt.*velocity;
allPositions = zeros(N,1,3,TotalTime+1);
allPositions(:,:,:,1) = coordsPres;

%allVelocities = zeros(N,1,3,TotalTime+1);
%allVelocities(:,:,:,1) = velocity;


for t=1:TotalTime
    [Accelerations] = getAccelerations(distances,closestVectors,N, velocity);
%    disp(Accelerations);
%    return;
    % Accelerations is Nx1x3
    [coordsPast, coordsPres, newVelocity] = evolveSystem(Accelerations,coordsPres,coordsPast,dt,L);
    allPositions(:,:,:,t+1) = coordsPres;
%    allVelocities(:,:,:,t+1) = newVelocity;
    velocity = newVelocity;
    [distances,closestVectors] = getDistances(coordsPres,L);
    
%    if(rand() < nudt)
%        vsum = sum(velocity)/N;
        % Shift velocities so that the momentum center of mass is at the center. 
%        velocity = velocity - vsum;

        % Function to get magnitude of velocity squared. 
%        magVelocity = @(vx, vy, vz) vx.^2 + vy.^2 + vz.^2;

        % Sum of velocity squares
%        vSquare = sum(magVelocity(velocity(:,1,1),velocity(:,1,2),velocity(:,1,3)));

        % Initial kinetic energy
%        KEinit = 0.5*vSquare;

        % Re-scale factor to get the desired temperature
%        vRescale = sqrt(initialKE/KEinit);
        % Scale the velocities by re-scale
%        velocity = velocity*vRescale;
%    end
    [Esystem, magV] = computeEnergy(N, distances, newVelocity);
    
%    fprintf('closest\n');
%    disp(closestVectors);
%    fprintf('distances\n');
%    disp(distances);
%    fprintf('acc\n');
%    disp(Accelerations);
%    disp(newVelocity);
%    return;
    Tsys = 2*Esystem/(3*N-3);
    if(t*dt > 0.4608 & mod(t,10) == 0)
        fprintf(fid,'%.4f \t %.4f \t %.4f\n', t*dt, Esystem, Tsys);
    end
    if(mod(t/TotalTime,0.25) == 0)
        fprintf(fid2, '%.4f\n', magV);
    end
end