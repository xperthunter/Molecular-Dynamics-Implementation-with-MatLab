function [distances,closestVectors] = getDistances(coordsPres, L)

DiffFunction = @(x1,x2)sign(x2-x1).*min(max(x2,x1)-min(x2,x1),min(L+x1,L+x2)-max(x1,x2));

positionsPermuted = permute(coordsPres,[2 1 3]);
closestVectors = bsxfun(DiffFunction,coordsPres,positionsPermuted);
distances = closestVectors.^2;
distances = sum(distances, 3);
distances = sqrt(distances);