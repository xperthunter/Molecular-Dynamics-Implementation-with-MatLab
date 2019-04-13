function [Accelerations] = getAccelerations(distances, closestVectors, N, velocity)

%Kdrag = 1e-2;
q = 0.5*eye(N);
velocity = 2*velocity;
distances = distances + eye(N);
aa = distances < 0.8;
bb = distances > 0.8;
distances = aa.*61 + distances.*bb;
Forces = (1./(distances.^12)) - 0.5*(1./(distances.^6));
Forces = Forces.*(1./(distances.^2));
Forces = Forces - q;
Forces = 48*Forces;
%fprintf('Forces\n');
%disp(Forces);
Accelerations = zeros(N,N,3);
Accelerations = Forces.*closestVectors;
%fprintf('before sum\n');
%disp(Accelerations);
Accelerations = sum(Accelerations,1);
Accelerations = permute(Accelerations, [2 1 3]);
Accelerations = Accelerations; % - Kdrag.*(velocity.^2);

% Accelerations are Nx1x3;

