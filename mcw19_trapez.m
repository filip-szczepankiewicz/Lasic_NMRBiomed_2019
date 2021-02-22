function g_lobe = mcw19_trapez(x,t1,t2,ramp)


% FILIP REWRITE
% This code is more robust to arbitary timing/slew/moment limits, and 
% handles triangular waveform components.
% Previous implementation could render g and t vectors of different size.
% To do, remove unnecessary input.

dt = x(2)-x(1);

n_rmp = ceil(ramp/dt);
n_max = floor(numel(x)/2);

n_use = min([n_rmp n_max]);

rmp = linspace(0, 1, ceil(ramp/dt));
rmp = rmp(1:n_use);


g_lobe = nan(size(x));
g_lobe(1:(length(rmp)-1)) = rmp(2:end); % here we remove the final zero, so that stacking does not cause flat g(t).
g_lobe = flip(g_lobe);
g_lobe(1:length(rmp)) = rmp;
g_lobe(isnan(g_lobe)) = max(rmp(:));



% OLD CODE
% N1 = find(t1<=x & x<=t1+ramp);
% N2 = find(t1+ramp<x & x<t2-ramp);
% N3 = find(t2-ramp<=x & x<=t2);
% 
% rampUP = (x-t1)/ramp;
% rampDown = 1-(x-t2+ramp)/ramp;
% g_lobe = [rampUP(N1) ones(1,numel(N2)) rampDown(N3)];
