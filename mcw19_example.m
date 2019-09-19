clear

%% Define premise
gmr = 26.75e7; % gamma (gyromagnetic ratio)
B0  = 3;

% M = 0; N = 380;
% M = 1; N = 542;
M = 2; N = 684;

waveForm = 2; %1-simusoidal 2-trapezoidal
ramp = 1.6e-3; %for trapezoidal only

dt = 1e-4; % sampling interval
d =  6.505e-3/(N*dt); %separation between encoding blocks (fraction of total encoding time)

check_LTE = 0;

% set shape - scale gradients in the initialization step
% for 3D they can be [1 1 1], but for 2D with M-nulling they need to be [1 1 0]
% to get required b values, post-initialization scalling is needed
% however, it seems we get a bit lower gmax and slew if we set shape in the initialization step

scale_gradient = [1 1 1]; % any 3D shape
%scale_gradient = [0.25 0.25 0.5]; % prolate
%scale_gradient = [0.4 0.4 0.2]; % oblate
%scale_gradient = [1 1 0]; % 2D

% this is needed to properly scale b after initialization
set_b = 1;
bval = 1e9; % 564*1e6
b = sort(scale_gradient)/sum(scale_gradient)*bval;
%b = [1 1 0]; b = sort(b)/sum(b)*bval;

order = [1 2 3]; % to permute order

autoFindRotation = 1; % auto rotate
if prod(scale_gradient) == 0
    autoFindRotation = 0;
end
NrotGrid = 10;
Niter = 10;

% rotate - Euler
alpha = 0*pi/4; %Z
beta = 0*pi/4;  %Y
gamma = 0*pi/4; %Z

% rotate - ZYX
Ax = 0*pi/4; % for Wacky 5 (0.404, not 0.392119*), for Wacky 6 (0.0028),  for Wacky 7 (0.214)
Ay = 0*pi/4;
Az = 0*pi/4;



%% Create waveform

Nhalf = round(0.5*N)+1; % up to 180
N0 = 2*round(d/(1-d)*Nhalf); % d is fraction of T
N1 = Nhalf-N0/2; % in the waveform
t1 = 0:dt:(N1-1)*dt; % in the waveform
T1 = t1(end); % duration of waveform within an encoding block
t = 0:dt:(Nhalf-1)*dt; % up to 180
T180 = t(end);
t = [t T180 + dt + t]';

if (check_LTE)
    
    g_tmp1 = fMCwaveform2x2(t1,0.02,0.24);
    
    g(:,1) = [g_tmp1 zeros(N0,1); -fliplr(g_tmp1)];
    
    g(:,2) = 0*g(:,1);
    g(:,3) = 0*g(:,1);
    
else
    
    switch M
        case 0
            intervals = [];
        case 1
            intervals = [1]/2;
        case 2
            intervals = [1 3]/4; % regular - even (selected by Johan)
            %intervals = [0.1 0.9]; %[1 3]/4 Wacky 1
            %intervals = [0.4 0.6]; %[1 3]/4 Wacky 2
            %intervals = [0.4 0.6]; %[1 3]/4 Wacky 3
            %intervals = [0.5 3.5]/4; %[1 3]/4 Wacky 4
            %intervals = [1 2]/4; %[1 3]/4 Wacky 5 asymmetric in z (selected by Johan)
            %intervals = [1 3]/4; %[1 3]/4 Wacky 6 asymmetric in xy
            %intervals = [1 2]/8; %[1 3]/4 Wacky 7 asymmetric in z and xy
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
            intervals = [1 3 5]/6; % regular - even (selected by Johan)
            %intervals = [0.4 0.5 0.6]; %[1 3 5]/6 Wacky 1
            %intervals = [0.1 0.5 0.9]; %[1 3 5]/6 Wacky 2
            %intervals = [0.4 0.5 0.6]; %[1 3 5]/6 Wacky 3
            %intervals = [1 3 5]/6; %[1 3 5]/6 Wacky 4
            %intervals = [1 3 5]/6; %[1 3 5]/6 Wacky 5 (selected by Johan)
            %intervals = [1 2 5]/6; %[1 2 5]/6 Wacky 6
            %intervals = [1 2 5]/6; %[1 2 5]/6 Wacky 7
            
    end
    
    if waveForm == 1
        g_tmp1 = mcw19_MCWaveform(t1, M, intervals);
    else
        g_tmp1 = mcw19_MCWaveform(t1, M, intervals, ramp/T1);
    end
    g(:,2) = scale_gradient(2)*[g_tmp1 zeros(1,N0) g_tmp1];
    g(:,3) = scale_gradient(3)*[g_tmp1 zeros(1,N0) -g_tmp1];
    
    
end


h = sign(t-T180); % dephasing direction (flips after 180)

% change order of axis
g = g(:,order);

%adjust g for b (needs to be before final rotation!)
if set_b%prod(sum(g)) ~= 0
    g = mcw19_scaleG(g, b, dt);
end


% rotate
if autoFindRotation
    A = mcw19_findRotation(g, h, 0, 0, 0, 2*pi, 2*pi, 2*pi, NrotGrid, Niter);
    Autox = A(1);
    Autoy = A(2);
    Autoz = A(3);
    R = mcw19_rotXYZ(Autox,Autoy,Autoz);
    g = mcw19_vectorRot(g, Autox, Autoy, Autoz);
    
    disp(['-------------------------------------------'])
    disp(['Ax = ' num2str(Autox) ' Ay = ' num2str(Autoy) ' Az = ' num2str(Autoz)])
    
end

g = mcw19_vectorRot(g, Ax, Ay, Az);

% the actual gradient - not the "effective" one
glab = g.*mcw19_repVec(h,size(g));


%% Plot result

try
    % This plot requires the MD-dMRI framework
    % https://github.com/markus-nilsson/md-dmri
    gwf_plot_all(g, h, dt)
    
catch
    mcw19_plotWF(t, g)
end
