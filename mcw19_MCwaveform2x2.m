function g = mcw19_MCwaveform2x2(tvec, ramp, delta1)
% symmetric LTE
%(Welsh CL, DiBella EV and Hsu EW. Higher-Order Motion-Compensation for In Vivo Cardiac Diffusion Tensor Imaging in Rats. IEEE Trans Med Imaging. 2015;34:1843-53)
% (+g1)-(-g2)-d-(+g2)-(-g1)
% (delta1+2*ramp)-(delta2+2*ramp)-d-(delta2+2*ramp)-(delta1+2*ramp)
% M0, M1 and M2 compensated but requires fliped repetition after 180

T = max(tvec);
ramp = T*ramp;
delta1 = T*delta1;
Delta1 = delta1 + 2*ramp;

%correct
%delta2 = Delta1 - 2*ramp +(2*T - 4*Delta1 + ramp)/2 - (8*Delta1^2 - 16*Delta1*T + 4*T^2 + 4*T*ramp + ramp^2)^(1/2)/2

%also correct
d = 2*(T-Delta1-delta2-2*ramp);
Delta = Delta1+delta2+2*ramp+d;
delta2 = -(2*Delta*ramp - 3*Delta1*ramp + 2*ramp^2 - Delta*Delta1)/(Delta - 2*Delta1 + ramp)

Delta2 = (Delta1*(Delta - ramp))/(Delta - 2*Delta1 + ramp);
delta2 = Delta2-2*ramp;

% ﻿Welsh CL, Dibella EVR, Member S, Hsu EW, Cardiac Diffusion Tensor Imaging in Rats, Trans. Med. Imaging, 2015, 34, (9), 1843–1853.
%     Delta2 = Delta1*(Delta-ramp)/(Delta-2*Delta1+ramp);
%     delta2 = Delta2 -2*ramp;

g = [];
for n = 1:length(tvec)
    t = tvec(n);
    if t>=0 & t<= ramp
        g = [g t/ramp];
    end
    
    t1 = ramp;
    if t>t1 & t< t1 + delta1
        g = [g 1];
    end
    
    t1 = ramp+delta1;
    if t>=t1 & t< t1 + ramp
        g = [g 1-(t-t1)/ramp];
    end
    
    t1 = 2*ramp + delta1;
    if t>=t1 & t< t1 + ramp
        g = [g -(t-t1)/ramp];
    end
    
    t1 = 3*ramp + delta1;
    if t>=t1 & t< t1 + delta2
        g = [g -1];
    end
    
    t1 = 3*ramp + delta1 + delta2;
    if t>=t1 & t< t1 + ramp
        g = [g -1+(t-t1)/ramp];
    end
    
    t1 = 4*ramp + delta1 + delta2;
    if t>=t1
        g = [g 0];
    end
end
end


