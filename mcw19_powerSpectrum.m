function [f PS] = mcw19_powerSpectrum(F,dt,N0,thresh,fmax)
% obtain the normalized power spectrum of F at positive frequencies
% truncate at cumulative power < thresh
% if fmax is supplied, truncate at fmax and ignore thresh value

% zero padding
F = [F; zeros(N0,1)];
NT = length(F);
f = 1/dt*(0:(NT/2))'/NT;

SF = fft(F)/NT;
%PStmp = SF.*conj(SF);
%PStmp = PStmp(1:NT/2+1);
%PStmp(2:end-1) = 2*PStmp(2:end-1);
PS = fftshift(SF.*conj(SF));
PS = [PS; PS(1)];
PS = PS(round(end/2+1):end);

% normalize
norm = sum(PS);
if norm ~= 0 PS = PS/norm; end

if nargin==4
    ind = cumsum(PS) <= thresh;
end

if nargin==5
    ind = f <= fmax;
end

if nargin>3
    PS = PS(ind);
    f = f(ind);
end
