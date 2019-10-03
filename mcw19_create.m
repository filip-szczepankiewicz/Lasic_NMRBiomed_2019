function [gwf, rf, dt] = mcw19_create(M, Gmax, Smax, ttot, tp, dt)

if nargin<1
    M = 0;
    
    Gmax = 80e-3;
    Smax = 100;
    
    dt = 1e-4;
    ttot = 60e-3;
    tp = 10e-3;
end

gmr = 26.75e7; % gamma (gyromagnetic ratio)
B0  = 3;

waveForm = 2; %1-simusoidal 2-trapezoidal

scale_gradient = [1 1 1]; % any 3D shape

order = [2 3 1]; % to permute axes

set_b = 1;
b = sort(scale_gradient)/sum(scale_gradient);

autoFindRotation = 0; % auto rotate

if prod(scale_gradient) == 0
    autoFindRotation = 0;
end

NrotGrid = 10;
Niter = 10;


%% Create waveform
ramp = Gmax/Smax; %for trapezoidal only

N = round(ttot/dt);
Nhalf = round(0.5*N)+1; % up to 180
N0 = 2*round(tp/dt/2); % d is fraction of T
N1 = Nhalf-N0/2; % in the waveform
t1 = 0:dt:(N1-1)*dt; % in the waveform
T1 = t1(end); % duration of waveform within an encoding block
t = 0:dt:(Nhalf-1)*dt; % up to 180
T180 = t(end);
t = [t T180 + dt + t]';

switch M
    case 0
        intervals = [];
    case 1
        intervals = [1]/2;
    case 2
        intervals = [1 3]/4; % regular - even (selected by Johan)
end

if waveForm == 1
    g_tmp1 = mcw19_MCWaveform(t1, M-1, intervals);
else
    g_tmp1 = mcw19_MCWaveform(t1, M-1, intervals, ramp/T1);
end

g(:,1) = scale_gradient(1)*[g_tmp1 zeros(1,N0) -g_tmp1];


switch M
    case 0
        intervals = [1]/2;
        
    case 1
        intervals = [1 3]/4;
        
    case 2
        intervals = [1 3 5]/6;
        
end

if waveForm == 1
    g_tmp1 = mcw19_MCWaveform(t1, M, intervals);
else
    g_tmp1 = mcw19_MCWaveform(t1, M, intervals, ramp/T1);
end
g(:,2) = scale_gradient(2)*[g_tmp1 zeros(1,N0) g_tmp1];
g(:,3) = scale_gradient(3)*[g_tmp1 zeros(1,N0) -g_tmp1];


h = sign(t-T180); % dephasing direction (flips after 180)

% change order of axis
g = g(:,order);

%adjust g for b (needs to be before final rotation!)
if set_b %prod(sum(g)) ~= 0
    g = mcw19_scaleG(g, b, dt, gmr);
end

% rotate
if autoFindRotation
    A = mcw19_findRotation(g, h, 0, 0, 0, 2*pi, 2*pi, 2*pi, NrotGrid, Niter);
    Autox = A(1);
    Autoy = A(2);
    Autoz = A(3);
    R = mcw19_rotXYZ(Autox,Autoy,Autoz);
    g = mcw19_vectorRot(g, Autox, Autoy, Autoz);
    
    disp('-------------------------------------------')
    disp(['Ax = ' num2str(Autox) ' Ay = ' num2str(Autoy) ' Az = ' num2str(Autoz)])
    
end


%% Compile output that is compatible with the MD-dMRI framework
gwf = g.*h / max(abs(g(:))) * Gmax;
rf  = h;


% gwf_plot_all(gwf, rf, dt)
% mcw19_plotWF(t, gwf)



