function A = mcw19_findRotation1(g,h,N)

Ax = 2*pi*linspace(0,1,N);
Ay = Ax;
Az = Ax;

n = size(g);

dmin = 1e12;
A = [0 0 0]';

for i = 1:N
    for j = 1:N
        for k = 1:N
            gtmp = fVectorRot(g, Ax(i), Ay(j), Az(k));
            
            % the actual gradient - not the "effective" one
            %glab = gtmp.*frepvec(h,n);
            
            
            K(1) = sum(h.*gtmp(:,1).*gtmp(:,1));
            K(2) = sum(h.*gtmp(:,2).*gtmp(:,2));
            K(3) = sum(h.*gtmp(:,3).*gtmp(:,3));
            K(4) = sum(h.*gtmp(:,1).*gtmp(:,3));
            K(5) = sum(h.*gtmp(:,2).*gtmp(:,3));
            
            %d = sum(K.^2);
            d = sum(abs(K));
            if d<dmin
                dmin = d;
                A(1) = Ax(i);
                A(2) = Ay(j);
                A(3) = Az(k);
            end
            
        end
    end
end

