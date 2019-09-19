function R = mcw19_rotXYZ(phiX,phiY,phiZ)
% Rotation matrix around XYZ

Rx = [1 0 0
    0 cos(phiX) -sin(phiX)
    0 sin(phiX) cos(phiX)];

Ry = [
    cos(phiY) 0 sin(phiY)
    0 1 0
    -sin(phiY) 0 cos(phiY)];

Rz = [
    cos(phiZ) -sin(phiZ) 0
    sin(phiZ) cos(phiZ) 0
    0 0 1];

R = Rz*Ry*Rx;

