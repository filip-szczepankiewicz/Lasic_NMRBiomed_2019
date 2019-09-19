function res = mcw19_vectorEulerRot(v, alpha, beta, gamma);
% rotates a vector arround the Euler angles alpha, beta, gamma (ZYZ rotations)
% angles or vector components can be vectors but not both
% size(v) = 3 x N
% ZYZ*v
if size(v,1) > size(v,2)
    v = v';
end

v1 = v(1,:);
v2 = v(2,:);
v3 = v(3,:);

N = length(alpha);
sA = sin(alpha);
cA = cos(alpha);
sB = sin(beta);
cB = cos(beta);
sG = sin(gamma);
cG = cos(gamma);

res(:,1) = conj(v2).*(cA.*sG + cB.*cG.*sA) - conj(v1).*(sA.*sG - cA.*cB.*cG) - cG.*sB.*conj(v3);
res(:,2) = conj(v2).*(cA.*cG - cB.*sA.*sG) - conj(v1).*(cG.*sA + cA.*cB.*sG) + sB.*sG.*conj(v3);
res(:,3) = cB.*conj(v3) + cA.*sB.*conj(v1) + sA.*sB.*conj(v2);
end

