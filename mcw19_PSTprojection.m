function projection = mcw19_PSTprojection(PS,n)
% calculate projection of power spectrum tensor along (n.x,n.y,n.z)
% PS(frq,m,n): frq channels, m,n = 1,2,3

projection = PS(:,1,1)'.*n.x.^2 + PS(:,2,2)'.*n.y.^2 + PS(:,3,3)'.*n.z.^2 ...
    + PS(:,1,2)'.*n.x.*n.y + PS(:,1,3)'.*n.x.*n.z + PS(:,2,3)'.*n.y.*n.z ...
    + PS(:,2,1)'.*n.x.*n.y + PS(:,3,1)'.*n.x.*n.z + PS(:,3,2)'.*n.y.*n.z;
