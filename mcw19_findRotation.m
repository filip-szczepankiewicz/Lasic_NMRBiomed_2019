function A = mcw19_findRotation(g, h, Ax0, Ay0, Az0, dAx, dAy, dAz, N, Nlevels)
n = size(g);

while Nlevels > 0
    Nlevels = Nlevels-1;
    Ax = Ax0 + dAx/2*linspace(-1,1,N);
    Ay = Ay0 + dAy/2*linspace(-1,1,N);
    Az = Az0 + dAz/2*linspace(-1,1,N);
    
    dAx = Ax(2)-Ax(1);
    dAy = Ay(2)-Ay(1);
    dAz = Az(2)-Az(1);
    
    dmin = 1e12;
    A = [0 0 0]';
    
    for i = 1:N
        for j = 1:N
            for k = 1:N
                gtmp = mcw19_vectorRot(g, Ax(i), Ay(j), Az(k));
                
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
                    Ax0 = Ax(i);
                    Ay0 = Ay(j);
                    Az0 = Az(k);
                end
            end
        end
    end
end
A(1) = Ax0;
A(2) = Ay0;
A(3) = Az0;

