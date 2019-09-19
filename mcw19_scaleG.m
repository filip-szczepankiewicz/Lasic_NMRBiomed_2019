function gout = mcw19_scaleG(g, b, dt, gmr)
%scale and rotate g(:,1:3) so that we get desired eigenvalues b(1:3)

q = gmr*cumsum(g)*dt;

B = q'*q*dt;

[V, L] = eig(B);
l = diag(L)';

gr = 0*g;

for i = 1:size(g,1)
    gr(i,:) = V'*g(i,:)';
end

gr = gr./sqrt(l).*sqrt(b);
gr(gr==0) = 0;


for i = 1:size(gr,1)
    gout(i,:) = V*gr(i,:)';
end
