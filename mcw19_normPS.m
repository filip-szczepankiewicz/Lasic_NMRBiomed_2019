function [f, PS] = mcw19_normPS(F,dt,mode)
% Power Spectrum of F

N = size(F,2);
NT = size(F,1);

f = 1/dt*(0:(NT/2))/NT;

if (exist('mode', 'var'))
    if strcmp(mode,'positive')
        mode = 1;
        PS = zeros(length(f),N,N);
    end
else
    mode = 0;
end

if ~mode
    f = [-fliplr(f(2:end)) f(1:end)];
    PS = zeros(NT+1,N,N);
end


S = zeros(NT,N);


for n = 1:N % spectra
    S(:,n) = fft(F(:,n));%/NT;
    ps = S(:,n).*conj(S(:,n));
    
    if mode
        PStmp = ps(1:length(f));
    else
        PStmp = fftshift(ps);
        PStmp = [PStmp; PStmp(1)];
        % PStmp = PStmp(1:NT/2+1);
        % PStmp(2:end-1) = 2*PStmp(2:end-1);
    end
    
    Norm(n) = sum(PStmp);
    if Norm(n) ~= 0
        PStmp = PStmp/Norm(n);
    end
    
    PS(:,n,n) = PStmp;
    
end


for m = 1:N
    for n = 1:N
        if m ~= n
            ps = S(:,m).*conj(S(:,n));
            
            if mode
                PStmp = ps(1:length(f));
            else
                PStmp = fftshift(ps);
                PStmp = [PStmp; PStmp(1)];
            end
            
            if Norm(m)*Norm(n) == 0 ; PStmp = PStmp*0; else
                PStmp = PStmp/sqrt(Norm(m)*Norm(n));
            end
            
            PS(:,m,n) = PStmp;
            
        end
    end
end





