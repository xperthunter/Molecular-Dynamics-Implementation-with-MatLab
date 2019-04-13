function D = densities(N,tau)

dilations = [0.2];

D = zeros(length(dilations),1);
%L = zeros(length(N),9);
for ii=1:length(N)
    nCube = 2;
    while(nCube^3 < N(ii))
        nCube = nCube + 1;
    end
    maxD = N(ii)/(((nCube - 1)*tau)^3);
%    D(ii,1) = maxD;
%    L(ii,1) = (N(ii)/maxD)^(1/3)+tau;
    for d=1:length(dilations)
        D(d,ii) = maxD*dilations(d);
%        L(ii,d+1) = (N(ii)/D(ii,d+1))^(1/3)+tau;
    end
end