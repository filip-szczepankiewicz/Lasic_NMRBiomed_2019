function gij = mcw19_fouter2(g)
tim = size(g,1);
dir = size(g,2);

gij = zeros(tim,dir,dir);

for i = 1:dir
    for j = 1:dir
        gij(:,i,j) = g(:,i).*g(:,j);
    end
end

end