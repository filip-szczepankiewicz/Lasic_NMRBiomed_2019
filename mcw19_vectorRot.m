function vout = mcw19_vectorRot(v, Ax, Ay, Az)
%rotate vectors v(n,3) around ZYX

v1 = v(:,1);
v2 = v(:,2);
v3 = v(:,3);

sAx = sin(Ax);
cAx = cos(Ax);
sAy = sin(Ay);
cAy = cos(Ay);
sAz = sin(Az);
cAz = cos(Az);

vout(:,1) = conj(v3).*(sAx.*sAz + cAx.*cAz.*sAy) - conj(v2).*(cAx.*sAz - cAz.*sAx.*sAy) + cAy.*cAz.*conj(v1);
vout(:,2) = conj(v2).*(cAx.*cAz + sAx.*sAy.*sAz) - conj(v3).*(cAz.*sAx - cAx.*sAy.*sAz) + cAy.*sAz.*conj(v1);
vout(:,3) = cAy.*sAx.*conj(v2) - sAy.*conj(v1) + cAx.*cAy.*conj(v3);