function gout = mcw19_scaleG(g, b,dt)
%scale and rotate g(:,1:3) so that we get desired eigenvalues b(1:3)
gmr = 26.75e7;
q = gmr*cumsum(g)*dt;

for i = 1:3
    for j = 1:3
        B(i,j) = sum(q(:,i).*q(:,j))*dt;
    end
end

[V L] = eig(B);
l = [L(1,1) L(2,2) L(3,3)];
%V*L*V'-B

gr = 0*g;

for i = 1:size(g,1)
    gr(i,:) = V'*g(i,:)';
end

ind = find(gr == 0);

gr = gr./sqrt(l).*sqrt(b);
gr(ind) = 0;


for i = 1:size(gr,1)
    gout(i,:) = V*gr(i,:)';
end

end


