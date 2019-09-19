function g_lobe = mcw19_trapez(x,t1,t2,ramp)
N1 = find(t1<=x & x<=t1+ramp);
N2 = find(t1+ramp<x & x<t2-ramp);
N3 = find(t2-ramp<=x & x<=t2);

rampUP = (x-t1)/ramp;
rampDown = 1-(x-t2+ramp)/ramp;
g_lobe = [rampUP(N1) ones(1,numel(N2)) rampDown(N3)];
end
