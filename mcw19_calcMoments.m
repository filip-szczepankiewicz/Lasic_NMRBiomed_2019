function m = mcw19_calcMoments(g, tmat, dt)
% function m = mcw19_calcMoments(g, tmat, dt)

for i = 0:5    
    m(i+1) = cumsum(g.*tmat.^i)*dt;
end