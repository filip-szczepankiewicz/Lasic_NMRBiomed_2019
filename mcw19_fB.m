function [V, L] = fB(q, dt)
for i = 1:3
    for j = 1:3
        B(i,j) = dt*sum(q(:,i).*q(:,j));
    end
end

[V L] = eig(B);
end