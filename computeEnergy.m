function [E,magV] = computeEnergy(N, distances, velocity)

distances = distances + eye(N);
dummy = (1./(distances.^12)) - (1./(distances.^6));
dummy = 4*dummy;
ELJ = sum(sum(dummy));

ELJ = 0.5*ELJ;

magVelocity = @(vx, vy, vz) vx.^2 + vy.^2 + vz.^2;
magV = magVelocity(velocity(:,1,1),velocity(:,1,2),velocity(:,1,3));
magV = sqrt(magV);

% Sum of velocity squares
vSquare = sum(magV.^2);
Ekinetic = 0.5*vSquare;

E = ELJ + Ekinetic;
