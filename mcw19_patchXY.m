function [px,py] = mcw19_patchXY(x,y1,y2,n)
l = length(x);
j = 1;
for i = 1:n:(l-n)
    
    px(1,j) = x(i);
    py(1,j) = y1(i);
    
    px(2,j) = x(i);
    py(2,j) = y2(i);
    
    px(3,j) = x(i+n);
    py(3,j) = y2(i+n);
    
    px(4,j) = x(i+n);
    py(4,j) = y1(i+n);
    
    j = j + 1;
end