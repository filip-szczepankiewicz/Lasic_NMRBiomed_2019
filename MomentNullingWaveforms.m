clear all

gmr = 26.75e7; % gamma (gyromagnetic ratio)
format long

saveWaveform = 0;
saveImage = 0;
showLabels = 1;
GreyScale = 0;
ImageDir = 'png/';

ImageName = 'Sin_Knull_Planar_rotX45'; %'TrapCalc_Knull';
ImageDir = 'manuscript_figs/';

%N = %380; 542; 684
% M = 0; N = 380;
% M = 1; N = 542;
M = 2; N = 684;

loadWave = 0;
filename_root= '/Volumes/Data/Dropbox (RWI)/Documents/Projects/Irvin_Teh_FlowAccelerationComp3D/Matlab/waveforms_trap_July2019/';

filename = [filename_root 'M0DW_STE_wf; b = 522.mat']; %
%filename = [filename_root 'M0DW_STE_wf; b = 564.mat']; % too large g - don't use)
%filename = [filename_root 'M1DW_STE_wf; b = 525.mat'];
%filename = [filename_root 'M2DW_STE_wf; b = 521.mat'];
% filename = [filename_root 'Welsh_LTE_wf; b = 1086.mat'];
% filename = [filename_root 'BP_LTE_wf; b = 732.mat'];
% filename = [filename_root 'PGSE_LTE_wf; b = 663.mat'];

waveForm = 1; %1-simusoidal 2-trapezoidal
ramp = 1.6e-3; %for trapezoidal only

dt = 1e-4; % sampling interval
d =  6.505e-3/(N*dt); %separation between encoding blocks (fraction of total encoding time)

check_LTE = 0;

% set shape - scale gradients in the initialization step
% for 3D they can be [1 1 1], but for 2D with M-nulling they need to be [1 1 0]
% to get required b values, post-initialization scalling is needed
% however, it seems we get a bit lower gmax and slew if we set shape in the initialization step

scale_gradient = [1 1 1]; % any 3D shape
% scale_gradient = [0.25 0.25 0.5]; % prolate
%scale_gradient = [0.4 0.4 0.2]; % oblate
%scale_gradient = [1 1 0]; % 2D

% this is needed to properly scale b after initialization
set_b = 1;
bval = 1e9; % 564*1e6
b = sort(scale_gradient)/sum(scale_gradient)*bval;
%b = [1 1 0]; b = sort(b)/sum(b)*bval;

order = [1 2 3]; % to permute order
%order = [3 2 1]; % to permute order

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


showSpectrum = 1;
showCrossTerms = 0; % look at cross-products of gradients

insetK = 0;
insetK_scale = 1.6; % for iso 1.2

% show q-trajectroy step by step
showSteps = 0;
step = 10;

% ploting parameters
showPatch = 1;
az_q = 65;
el_q = 40;
az_PS = 35;
el_PS = 40;

fs = 24;
lw = 3;
lw_curve = 3;
xmax = 1.1;

%wfm_plot_norms = [0.164733146713225  0.000641768217016 0.000034691736896 0.000002245160510]; % max
wfm_plot_norms = []; % if empty then find norm automatically




%%%%%%%%%%%%%%%%% waveform  %%%%%%%%%%

if loadWave
    
    if (1) % from mat file
        
        load(filename)
        g(:,1) = wf.all(1,:)';
        g(:,2) = wf.all(2,:)';
        g(:,3) = wf.all(3,:)';
        
        disp(['Ng = ' num2str(length(g))])
        disp(['maxG = ' num2str(max(abs(g(:))))])
        
        
    else % from txt file
        
        %filename = '/Users/Samo/Dropbox (RWI)/Documents/Projects/Irvin_Teh_FlowAccelerationComp3D/Matlab/Irvin_waveforms/FWF_CUSTOM001_A.txt';
        filename = '/Users/Samo/Dropbox (RWI)/Documents/Projects/Irvin_Teh_FlowAccelerationComp3D/Matlab/Irvin_waveforms/FWF_CUSTOM001_B.txt';
        
        delimiterIn = ' ';
        headerlinesIn = 1;
        data = importdata(filename,delimiterIn,headerlinesIn);
        data = data.data;
        N1 = size(data,1);
        N0 = round(2*d/(1-d)*N1); % d is fraction of T
        if (0)
            g(:,1) = data(:,3);
            g(:,2) = data(:,1);
            g(:,3) = data(:,2);
        else
            g(:,1) = [-data(:,1); zeros(N0,1); -data(:,1)];
            g(:,2) = [data(:,2); zeros(N0,1); -data(:,2)];
            g(:,3) = [data(:,3); zeros(N0,1); -data(:,3)];
        end
        
    end
    
    %         figure(1),clf
    %         hold on
    %         plot(g(:,1),'r.')
    %         plot(g(:,2),'go')
    %         plot(g(:,3),'bx')
    %         return
    
    N = size(g,1);
    disp(['N = ' num2str(N)]);
    T = (N-1)*dt;
    t = [0:dt:T]';
    
    T180 = T/2;
    
    
else % generate waveform
    
    Nhalf = round(0.5*N)+1; % up to 180
    N0 = 2*round(d/(1-d)*Nhalf); % d is fraction of T
    N1 = Nhalf-N0/2; % in the waveform
    t1 = 0:dt:(N1-1)*dt; % in the waveform
    T1 = t1(end); % duration of waveform within an encoding block
    t = 0:dt:(Nhalf-1)*dt; % up to 180
    T180 = t(end);
    t = [t T180 + dt + t]';
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LTE
    
    if (check_LTE)
        
        %g_tmp1 = fMCwaveform2x2(t1,0.00001,1/(2+sqrt(2)));
        g_tmp1 = fMCwaveform2x2(t1,0.02,0.24);
        
        %         figure(1),clf
        %         plot(g_tmp1)
        %         return
        g(:,1) = [g_tmp1 zeros(N0,1); -fliplr(g_tmp1)];
        %g(:,1) = [g_tmp1 -g_tmp1];
        
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
            g_tmp1 = fMCWaveform(t1, M-1, intervals);
        else
            g_tmp1 = fMCWaveform(t1, M-1, intervals, ramp/T1);
        end
        
        %g(:,1) = scale_gradient(1)*[g_tmp1 -g_tmp2];
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
            g_tmp1 = fMCWaveform(t1, M, intervals);
        else
            g_tmp1 = fMCWaveform(t1, M, intervals, ramp/T1);
        end
        %             g(:,2) = scale_gradient(2)*[g_tmp1  -g_tmp2];
        %             g(:,3) = scale_gradient(3)*[g_tmp1 g_tmp2];
        g(:,2) = scale_gradient(2)*[g_tmp1 zeros(1,N0) g_tmp1];
        g(:,3) = scale_gradient(3)*[g_tmp1 zeros(1,N0) -g_tmp1];
        
        
    end
end

h = sign(t-T180); % dephasing direction (flips after 180)

% change order of axis
g = g(:,order);

%adjust g for b (needs to be before final rotation!)
if set_b%prod(sum(g)) ~= 0
    g = fscaleG(g, b, dt);
end

disp(['-------------------------------------------'])
% rotate
if autoFindRotation
    A = findRotation(g, h, 0, 0, 0, 2*pi, 2*pi, 2*pi, NrotGrid, Niter);
    %A = findRotation1(g,h,NrotGrid);
    Autox = A(1);
    Autoy = A(2);
    Autoz = A(3);
    R = fRotXYZ(Autox,Autoy,Autoz);
    %round(R,2)
    g = fVectorRot(g, Autox, Autoy, Autoz);
    disp(['Ax = ' num2str(Autox) ' Ay = ' num2str(Autoy) ' Az = ' num2str(Autoz)])
    
end

g = fVectorRot(g, Ax, Ay, Az);

m0 = cumsum(g)*dt;
% figure(1),clf
% plot(t,m0), return
tmat = frepvec(t,size(g));
m1 = cumsum(g.*tmat)*dt;
m2 = cumsum(g.*tmat.^2)*dt;


% the actual gradient - not the "effective" one
glab = g.*frepvec(h,size(g));


gij = fouter2(glab);
h = frepvec(h,size(gij));
Maxwell_ij = squeeze(cumsum(gij.*h))*dt; % proportional to Filip's M

B0 = 3;
for i = 1:size(Maxwell_ij,1)
    K(i,:,:) = gmr / 2 / pi / (4 * B0)*[Maxwell_ij(i,3,3) 0 -2*Maxwell_ij(i,1,3);...
        0 Maxwell_ij(i,3,3) -2*Maxwell_ij(i,2,3);...
        -2*Maxwell_ij(i,1,3) -2*Maxwell_ij(i,2,3) 4*Maxwell_ij(i,1,1)+4*Maxwell_ij(i,2,2)]; % proportional to Filip's K
    
    Mt = squeeze(Maxwell_ij(i,:,:));
    MM(i,:,:) = Mt*Mt';
end


Maxwell = sqrt(MM(:,1,1) + MM(:,2,2) + MM(:,3,3)); % proportional to Filip's m
%disp(['m = ' num2str(Maxwell(end))]) % proportional to Filip's m


[V L] = fB(gmr*m0,dt);
L = L*1e-9;
disp(['eig(b)/1e9 = ' num2str(L(1,1)) ', ' num2str(L(2,2)) ', ' num2str(L(3,3)) ])


slew = diff(g)/dt;
max_slew = max(abs(slew(:)));
max_g = max(abs(g(:)));

disp(['max_g = ' num2str(max_g) ', max_slew = ' num2str(max_slew) ])
gnorm = g/max_g;
qnorm = gmr*cumsum(gnorm)*dt;
bnorm = dt*[sum(qnorm(:,1).*qnorm(:,1)) sum(qnorm(:,1).*qnorm(:,2)) sum(qnorm(:,1).*qnorm(:,3))
    sum(qnorm(:,1).*qnorm(:,2)) sum(qnorm(:,2).*qnorm(:,2)) sum(qnorm(:,2).*qnorm(:,3))
    sum(qnorm(:,1).*qnorm(:,3)) sum(qnorm(:,2).*qnorm(:,3)) sum(qnorm(:,3).*qnorm(:,3)) ];

disp(['----- per unit gradient ------'])
[V L] = eig(bnorm);
L = L*1e-9;
disp(['eig(b)/1e9 = ' num2str(L(1,1)) ', ' num2str(L(2,2)) ', ' num2str(L(3,3)) ])


slewnorm = diff(gnorm)/dt;
max_slewnorm = max(abs(slewnorm(:)));
max_gnorm = max(abs(gnorm(:)));
disp(['max_slewnorm = ' num2str(max_slewnorm) ])


if (1)
    % K
    figure(202),clf
    
    X = t*1000; %t/max(t);
    
    if GreyScale
        ph = plot(X,squeeze(K(:,1,1)),'k-',...
            X,squeeze(K(:,3,3)),'k--',...
            X,squeeze(K(:,1,3)),'k:',...
            X,squeeze(K(:,2,3)),'k-.','linewidth',lw_curve);
        
        lg = legend('K(1,1)','K(3,3)','K(1,3)','K(2,3)');
        legend boxoff
        lg.Position(1) = lg.Position(1)+0.05;
    else
        ph = plot(X,squeeze(K(:,1,1)),'r-',...
            X,squeeze(K(:,3,3)),'b--',...
            X,squeeze(K(:,1,3)),'m:',...
            X,squeeze(K(:,2,3)),'c-.','linewidth',lw_curve);
        hold on
        plot(X*xmax,0*X,'-k','linewidth',lw_curve)
        
        if showLabels
            lg = legend('K(1,1)','K(3,3)','K(1,3)','K(2,3)');
        else
            lg = legend('','','','');
        end
        legend boxoff
        lg.Location = 'best';
        %lg.Position(1) = lg.Position(1)+0.15;
    end
    
    
    maxY = max(K(:));
    minY = min(K(:));
    
    xlim([0 max(X)*xmax])
    ylim(1.1*[minY maxY])
    %ylabel('K_{ij}(t)')
    %xlabel('\itt \rm/ \tau')
    set(gca,'YTickLabel',[],'FontSize',fs, 'LineWidth',lw,'TickDir','out',...
        'TickLength',[.02 .02])%,'XAxisLocation','origin')
    box off
    
    % inset
    if (insetK)
        fh = axes(gcf,'Position',[0.3 0.57 0.4 0.23]);
        ph = plot(X,squeeze(K(:,1,1)),'r-',...
            X,squeeze(K(:,3,3)),'b--',...
            X,squeeze(K(:,1,3)),'m:',...
            X,squeeze(K(:,2,3)),'c-.','linewidth',lw_curve);
        
        xlim([0 max(X)*xmax])
        if maxY ~= 0
            ylim(insetK_scale*maxY*[-1 1])
        end
        
        %ylabel('K_{ij}(t)')
        %xlabel('t/\tau')
        set(gca,'XTickLabel',[], 'YTickLabel',[],'FontSize',fs, 'LineWidth',lw*0.75,'TickDir','out',...
            'TickLength',[.01 .01],'XAxisLocation','origin')
        box off
        %axis off
        
    end
    
    
    
    if saveImage
        s = sprintf('K_M%d_x0%0.2f', M, x180*10);
        s = [ImageDir ImageName '_' s   '.png'];
        eval(['print ' s ' -loose -dpng -r300'])
        pause(1);
    end
    
end


if GreyScale
    ls = {'k-','k-','k-'};
    lsm = {'k--','k-.','k:'};
else
    %     ls = {'r-','g-','b-'};
    %     lsm = {'c--','m-.','y:'};
    
    ls = {'k-','k-','k-'};
    lsm = {'r--','g-.','b:'};
    
end

labels = {'X','Y','Z'};

if isempty(wfm_plot_norms)
    norm_g = 1.05*max(abs(g(:)));
    norm_m0 = 1.05*max(abs(m0(:)));
    norm_m1 = 1.05*max(abs(m1(:)));
    norm_m2 = 1.05*max(abs(m2(:)));
else
    norm_g = wfm_plot_norms(1);
    norm_m0 = wfm_plot_norms(2);
    norm_m1 = wfm_plot_norms(3);
    norm_m2 = wfm_plot_norms(4);
end


for n = 1:3
    
    figure(n),clf
    hold on
    
    Np = 1; % avery Np-th point is used for pathing (sometimes can speed up things)
    
    X = t*1000; % t/max(t);
    
    %%%%%%   patches between 0 and Y
    %Y = [g(:,n)/norm_g m0(:,n)/norm_m0 m1(:,n)/norm_m1 m2(:,n)/norm_m2];
    Y = [g(:,n)/norm_g];
    
    col = 0.9*[1 1 1]';
    col(:,2)=col(:,1);
    col(:,3)=col(:,1);
    col(:,4)=col(:,1);
    
    for m = 1:size(Y,2)
        if (showPatch)
            [px,py] = patchXY(X,0*Y(:,m),Y(:,m),Np);
            p = patch(px,py,'k');
            p.FaceColor = col(:,m);
            p.FaceAlpha = 1;
            p.EdgeColor = 'none';
            p.LineStyle = '--';
            
        end
        
    end
    
    
    ph = plot(X,Y,ls{n},...
        X,m0(:,n)/norm_m0,lsm{1},...
        X,m1(:,n)/norm_m1,lsm{2},...
        X,m2(:,n)/norm_m2,lsm{3},'linewidth',lw_curve);
    
    plot(X*xmax,0*X,'-k','linewidth',lw_curve);
    
    
    xlim([0 max(X)*xmax])
    ylim([-1 1])
    if showLabels
        %xlabel('\itt \rm/ \tau')
        xlabel('\itt \rm[ms]')
        ylabel(labels{n})
    end
    
    set(gca,'YTickLabel',[],'FontSize',fs, 'LineWidth',lw,...
        'TickDir','out','TickLength',[.02 .02])%,'XAxisLocation','origin')
    %     set(gca,'YTickLabel',[],'FontSize',fs, 'LineWidth',lw,...
    %         'TickDir','out','TickLength',[.02 .02])
    box off
    
    
    
    if showLabels
        lg = legend(ph,'g','m0','m1','m2');
    else
        lg = legend(ph, '','','','');
    end
    legend boxoff
    lg.Position(1) = lg.Position(1)+0.1;
    
    
    if saveImage
        s = sprintf('wfm_%d_M%d_x0%0.2f', n, M, x180*10);
        s = [ImageDir ImageName '_' s '.png'];
        eval(['print ' s ' -loose -dpng -r150'])
        pause(1);
    end
    
    
    
end


if showCrossTerms
    
    figure(10),clf
    
    y1 = g(:,1).*g(:,2);
    y2 = g(:,1).*g(:,3);
    y3 = g(:,2).*g(:,3);
    
    
    if (0)
        y1 = y1/max(abs(y1(:)));
        y2 = y2/max(abs(y2(:)));
        y3 = y3/max(abs(y3(:)));
    else
        norm = max([y1; y2; y3]);
        y1 = y1/norm;
        y2 = y2/norm;
        y3 = y3/norm;
    end
    
    if (1)
        plot(t/max(t),cumsum(y1.*h(:,1,1)),'-',...
            t/max(t),cumsum(y2.*h(:,1,1)),'--',...
            t/max(t),cumsum(y3.*h(:,1,1)),':','linewidth',2)
    else
        plot(t/max(t),y1,'-',...
            t/max(t),y2,'--',...
            t/max(t),y3,':','linewidth',2)
    end
    
    
    xlim([0 xmax])
    %ylim([-1 1])
    xlabel('\itt \rm/ \tau')
    ylabel('cross products')
    
    set(gca,'YTickLabel',[],'FontSize',fs, 'LineWidth',lw,...
        'TickDir','out','TickLength',[.02 .02],'XAxisLocation','origin')
    %     set(gca,'YTickLabel',[],'FontSize',fs, 'LineWidth',lw,...
    %         'TickDir','out','TickLength',[.02 .02])
    box off
    
    lg = legend('xy','xz','yz');
    legend boxoff
    lg.Position(1) = lg.Position(1)+0.05;
    
    
    if saveImage
        s = sprintf('wfm_corss_M%d_x0%0.2f', M, x180*10);
        s = [ImageDir ImageName '_' s '.png'];
        eval(['print ' s ' -loose -dpng -r150'])
        pause(1);
    end
    
    
    
end



if (0)
    figure(4),clf
    %plot(t/max(t),b(:,1),'-r',t/max(t),b(:,2),'g--',t/max(t),b(:,3),'b:',t/max(t),Maxwell*scale_m,'-k','linewidth',lw)
    plot(t/max(t),Maxwell*scale_m,'-k','linewidth',lw)
    
    ylim([-1 1]*1.1)
    xlim([0 xmax])
    xlabel('t/\tau')
    ylabel('b(t), m(t)')
    set(gca,'FontSize',fs, 'LineWidth',lw,'TickDir','out','TickLength',[.02 .02])
    box off
    
    if saveImage
        s = sprintf('bANDm_M%d_x0%0.2f', M, x180*10);
        s = [ImageDir ImageName '_' s '.png'];
        eval(['print ' s ' -loose -dpng -r150'])
        pause(1);
    end
    
end

if (1)
    %lwQ = 6;% 4;
    L = max(sqrt(m0(:,1).^2+m0(:,2).^2+m0(:,3).^2));
    
    figure(100),clf
    hold on
    
    plot3(m0(:,1),m0(:,2),m0(:,3),'-k','linewidth',lw), xlabel('X'), ylabel('Y'), zlabel('Z')
    
    if GreyScale
        line(L*[-1 1],[0 0],[0 0],'Color','k','LineStyle','--','linewidth',2)
        line([0 0],L*[-1 1],[0 0],'Color','k','LineStyle','-.','linewidth',2)
        line([0 0],[0 0],L*[0 1],'Color','k','LineStyle',':','linewidth',2)
    else
        line(L*[-1 1],[0 0],[0 0],'Color','red','LineStyle','--','linewidth',2)
        line([0 0],L*[-1 1],[0 0],'Color','green','LineStyle','--','linewidth',2)
        line([0 0],[0 0],L*[0 1],'Color','blue','LineStyle','--','linewidth',2)
    end
   
    for n = 1:length(t)
        
        
        line([0 m0(n,1)],[0 m0(n,2)],[0 m0(n,3)],'Color',0.7*[1 1 1],'LineStyle','-')
        %         r = 1;
        %         colormix = [linspace(0,1,length(t)/2) linspace(1,0,length(t)/2)] ;
        %         rampDawn = linspace(1,0,length(t));
        %         rampUp = linspace(0,1,length(t));
        
        %         s = 'png/HeartArt/Heart_Art_BloodRed_3';
        %         line(r*[0 m0(n,1)],r*[0 m0(n,2)],r*[0 m0(n,3)],'Color',0.8*[1 0 0],'LineStyle','-','linewidth',1)
        %          s = 'png/HeartArt/Heart_Blue_Pink_3';
        %         line(r*[0 m0(n,1)],r*[0 m0(n,2)],r*[0 m0(n,3)],'Color',1*[rampUp(n) 0 1],'LineStyle','-','linewidth',1)
        %          s = 'png/HeartArt/Heart_Blue_Green_3';
        %         line(r*[0 m0(n,1)],r*[0 m0(n,2)],r*[0 m0(n,3)],'Color',1*[0 rampUp(n) rampDawn(n)],'LineStyle','-','linewidth',1)
        %          s = 'png/HeartArt/Heart_Red_Yellow_3';
        %         line(r*[0 m0(n,1)],r*[0 m0(n,2)],r*[0 m0(n,3)],'Color',1*[1 rampUp(n) 0],'LineStyle','-','linewidth',1)
        %          s = 'png/HeartArt/Heart_Green_Yellow_3';
        %         line(r*[0 m0(n,1)],r*[0 m0(n,2)],r*[0 m0(n,3)],'Color',1*[rampUp(n) 1 0 ],'LineStyle','-','linewidth',1)
        
        
    end
    %plot3(m0(:,1),m0(:,2),m0(:,3),'-k','linewidth',lwQ), xlabel('X'), ylabel('Y'), zlabel('Z')
    
    
    
    view(az_q,el_q)
    %     view(75,40)
    %     view(90,-60)
    
    
    %ax_lim = xlim;
    
    set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],...
        'FontSize',fs, 'LineWidth',lw,'TickDir','out','TickLength',[.02 .02])
    axis equal
    box off
    if ~showLabels
        axis off
    end
    %     eval(['print ' s ' -loose -dpng -r600'])
    %     eval(['print ' s ' -loose -dpdf -r600']), return
    
    
    if saveImage
        s = sprintf('q_M%d_x0%0.2f', M, x180*10);
        s = [ImageDir ImageName '_' s   '.png'];
        eval(['print ' s ' -loose -dpng -r300'])
        pause(1);
    end
    
end




if (0)
    % Mij = cumsum(gi*gj*h)
    figure(201),clf
    
    plot(t/max(t),squeeze(Maxwell_ij(:,1,1)),'r-',...
        t/max(t),squeeze(Maxwell_ij(:,2,2)),'g-',...
        t/max(t),squeeze(Maxwell_ij(:,3,3)),'b-',...
        t/max(t),squeeze(Maxwell_ij(:,1,2)),'c-',...
        t/max(t),squeeze(Maxwell_ij(:,1,3)),'m-',...
        t/max(t),squeeze(Maxwell_ij(:,2,3)),'y-','linewidth',lw)
    xlim([0 xmax])
    ylim([-1 1]*maxY_Maxwell)
    ylabel('M_{ij}(t)')
    xlabel('t/\tau')
    set(gca,'YTickLabel',[],'FontSize',fs, 'LineWidth',lw,'TickDir','out',...
        'TickLength',[.02 .02],'XAxisLocation','origin')
    box off
    
    if saveImage
        s = sprintf('Maxwell_M%d_x0%0.2f', M, x180*10);
        s = [ImageDir ImageName '_' s   '.png'];
        eval(['print ' s ' -loose -dpng -r300'])
        pause(1);
    end
    
end




if (0) % gi*gj
    figure(203),clf
    
    plot(t/max(t),squeeze(gij(:,1,1)),'r-',...
        t/max(t),squeeze(gij(:,2,2)),'g-',...
        t/max(t),squeeze(gij(:,3,3)),'b-',...
        t/max(t),squeeze(gij(:,1,2)),'y--',...
        t/max(t),squeeze(gij(:,1,3)),'m--',...
        t/max(t),squeeze(gij(:,3,2)),'c--','linewidth',lw)
    xlim([0 xmax])
    ylim([-1 1]*maxY_Maxwell/T)
    set(gca,'YTickLabel',[],'FontSize',fs, 'LineWidth',lw,...
        'TickDir','out','TickLength',[.02 .02],'XAxisLocation','origin')
    box off
    
end



if (showSteps)
    
    %[x, y, z] = fVectorEulerRot_abg(x,y,z,alpha,beta,gamma);
    scale = max(abs(m0(:)));
    figure(300),clf
    hold on
    hAL = animatedline('linewidth',4);
    
    %     hAL.LineStyle='-';
    %     hAL.Color='black';
    
    % defaults for the marker
    hAx=gca;      % get the axis handle
    sz=120;        % size of marker
    clr='red';      % color
    hS=scatter(hAx,nan,nan,sz,clr,'filled');  % initial point won't show but creates handle
    
    xlabel('X'),ylabel('Y'),zlabel('Z')
    set(gca,'XLim',scale*[-1 1],'YLim',scale*[-1 1],'ZLim',scale*[-1 1])
    
    view(az_q,el_q)
    %xlim = ax_lim;
    
    
    
    if step == 0
        addpoints(hAL,m0(:,1),m0(:,2),m0(:,3))
        drawnow
    else
        
        for k = 1:step:size(m0,1)
            addpoints(hAL,m0(k,1),m0(k,2),m0(k,3))
            set(hS,'xdata',m0(k,1),'ydata',m0(k,2),'zdata',m0(k,3))  % update the marker position
            drawnow
        end
    end
    
    set(gca,'FontSize',fs, 'LineWidth',lw,'TickDir','out','TickLength',[.02 .02])
    box off
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% power spectra
if (showSpectrum)
    %     addpath(genpath('/Users/Samo/Dropbox (RWI)/Documents/Matlab_Samo/MatlabUtilities'));
    %     load UDSRTriN1000
    
    
    
    ODF = fUDSR1000();
    u.x = ODF.x; u.y = ODF.y; u.z = ODF.z;
    M0shape = ones(numel(ODF.x),1);
    
    Nzeros = size(m0,1)*100;
    m0z(:,1) =[ m0(:,1); zeros(Nzeros,1)];
    m0z(:,2) =[ m0(:,2); zeros(Nzeros,1)];
    m0z(:,3) =[ m0(:,3); zeros(Nzeros,1)];
   
    [f PS] = fnormPS(m0z,dt,'positive');
 
    PSu = real(fPSTprojection(PS,u));
    CumAveragePower = cumsum(sum(PSu,1));
    CumAveragePower = CumAveragePower/CumAveragePower(end); % normalize
    ind = CumAveragePower <= 0.99;
    CumAveragePower = CumAveragePower(ind);
    PSu = PSu(:,ind); 
    f = f(ind);
    maxFrq = max(f);
    
       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RGB power spectra
   
    Ncolors = 300;%100, 64; % number of colors
    PSuColor = interp1(PSu',linspace(1,size(PSu,2),Ncolors))';
    CumAveragePowerColor = cumsum(sum(PSuColor,1));
    PSuColor = PSuColor/CumAveragePowerColor(end);
    CumAveragePowerColor = CumAveragePowerColor/CumAveragePowerColor(end);
 

%     figure(1),clf
%     plot(CumAveragePowerColor), return

%      PSuColor = PSuColor ./ repmat(sum(PSuColor,1),size(u.x));
%      PSuColor = PSuColor ./ repmat(sum(PSuColor,2),1,Ncolors);
    

%      figure(1),clf
%      plot(PSuColor), return

    
    % color frequency ranges
    Nred = max(find(CumAveragePowerColor<1/3));
    Ngreen = max(find(CumAveragePowerColor<2/3))-Nred;
    Nblue = Ncolors-Nred-Ngreen;
    % red-green-blue
    colFrq(:,1) = [1*ones(1,Nred) 0*ones(1,Ngreen) 0*ones(1,Nblue)];
    colFrq(:,2) = [0*ones(1,Nred) 1*ones(1,Ngreen) 0*ones(1,Nblue)];
    colFrq(:,3) = [0*ones(1,Nred) 0*ones(1,Ngreen) 1*ones(1,Nblue)];
    
    
    
    h = figure(1001);
    clf
    alpha = 1;
    factor = 0.8;
    patch('Faces',[1:4],'Vertices',...
        [0 0; max(find(CumAveragePowerColor<1/3)) 0; ...
        max(find(CumAveragePowerColor<1/3)) 1;...
        0 1],'EdgeColor','none','FaceColor',factor*[1 0 0],'FaceAlpha',alpha);
    patch('Faces',[1:4],'Vertices',...
        [max(find(CumAveragePowerColor<1/3)) 0; max(find(CumAveragePowerColor<2/3)) 0; ...
        max(find(CumAveragePowerColor<2/3)) 1;...
        max(find(CumAveragePowerColor<1/3)) 1],'EdgeColor','none','FaceColor',factor*[0 1 0],'FaceAlpha',alpha);
    patch('Faces',[1:4],'Vertices',...
        [max(find(CumAveragePowerColor<2/3)) 0; max(find(CumAveragePowerColor<1)) 0; ...
        max(find(CumAveragePowerColor<1)) 1;...
        max(find(CumAveragePowerColor<2/3)) 1],'EdgeColor','none','FaceColor',factor*[0 0 1],'FaceAlpha',alpha);
    
    hold on
    plot(1:length(colFrq),CumAveragePowerColor,'k-','Linewidth',6)
    xlim([0 length(colFrq)])
    %xlim([0 frq(max(find(CumAveragePowerColor<0.99)))])
    box off
    pbaspect([5 1 1])
    axis off
    
    if (saveImage)
        s = sprintf('PS_M%d_x0%0.2f', M, x180*10);
        s = [ImageDir ImageName '_' s  '.png'];
        eval(['print ' s ' -loose -dpng -r300'])
        pause(1);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for index = 1:size(PSuColor,1)
        C(index,:) = sum(length(u.x)*[PSuColor(index,:); PSuColor(index,:); PSuColor(index,:)]'.*colFrq);
    end
    
    
    % black spots
    if (0)
        d = C - mean(C);
        for index = 1:size(C,1)
            dd(index) = sqrt(sum(d(index,:).^2));
        end
        index = find(dd<8*min(dd));
        C(index,:) = 0;
    end
    
    alpha = 1;
    h = figure(1002);
    clf
    
    if (0) %%%%% make q trajectory
        hold on
        factor = 1.2;
        tmp = [F.x F.y F.z];
        tmp = tmp/max(tmp(:));
        plot3(factor*tmp(:,1),factor*tmp(:,2),factor*tmp(:,3),'-k','Linewidth',lw)
        alpha1 = 0.95;
    else
        alpha1 = alpha;
    end
    
    ODF.P = M0shape; %ones(1000,1);% (ODF.x + ODF.y + ODF.z).^2;
    ODF.verts = repmat(ODF.P,[1 3]).*[sin(ODF.theta).*cos(ODF.phi)...
        sin(ODF.theta).*sin(ODF.phi) ...
        cos(ODF.theta)];
    
    pa = patch('Faces',ODF.tri,'Vertices',ODF.verts/max(ODF.P));
    axis tight off, axis square, axis equal
    
    
    
    if (1)
        hold on
        L = 1.5;
        plot3(L*[-1:1],0*[-1:1],0*[-1:1],'-r','Linewidth',lw)
        plot3(0*[-1:1],L*[-1:1],0*[-1:1],'-g','Linewidth',lw)
        plot3(0*[-1:1],0*[-1:1],L*[-1:1],'-b','Linewidth',lw)
        
    end
    
    set(pa,'FaceColor','interp','FaceVertexCData',C,...
        'EdgeColor',.8*[1 1 1],'LineWidth',.5,'EdgeAlpha',.5,'FaceAlpha',alpha1)
    %title(['shape M_0' ', color spectra'],'FontSize',14)
    set(gca,'XTick',[],'YTick',[],'ZTick',[],'FontSize',14)
    view(az_PS,el_PS)
    
    if (saveImage)
        s = sprintf('RGB_Power_M%d_x0%0.2f', M, x180*10);
        s = [ImageDir ImageName '_' s  '.png'];
        print(h,s, '-dpng', '-r300');
    end
    
end

if saveWaveform
    s = sprintf('g_M%d_x0%0.2f.mat', M, x180*10);
    s = ['save ' ImageDir ImageName '_' s];
    eval(s)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% funtions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = fMCwaveform2x2(tvec, ramp, delta1)
% symmetric LTE
%(Welsh CL, DiBella EV and Hsu EW. Higher-Order Motion-Compensation for In Vivo Cardiac Diffusion Tensor Imaging in Rats. IEEE Trans Med Imaging. 2015;34:1843-53)
% (+g1)-(-g2)-d-(+g2)-(-g1)
% (delta1+2*ramp)-(delta2+2*ramp)-d-(delta2+2*ramp)-(delta1+2*ramp)
% M0, M1 and M2 compensated but requires fliped repetition after 180

T = max(tvec);
ramp = T*ramp;
delta1 = T*delta1;
Delta1 = delta1 + 2*ramp;

%correct
%delta2 = Delta1 - 2*ramp +(2*T - 4*Delta1 + ramp)/2 - (8*Delta1^2 - 16*Delta1*T + 4*T^2 + 4*T*ramp + ramp^2)^(1/2)/2

%also correct
d = 2*(T-Delta1-delta2-2*ramp);
Delta = Delta1+delta2+2*ramp+d;
delta2 = -(2*Delta*ramp - 3*Delta1*ramp + 2*ramp^2 - Delta*Delta1)/(Delta - 2*Delta1 + ramp)

Delta2 = (Delta1*(Delta - ramp))/(Delta - 2*Delta1 + ramp);
delta2 = Delta2-2*ramp;

% ﻿Welsh CL, Dibella EVR, Member S, Hsu EW, Cardiac Diffusion Tensor Imaging in Rats, Trans. Med. Imaging, 2015, 34, (9), 1843–1853.
%     Delta2 = Delta1*(Delta-ramp)/(Delta-2*Delta1+ramp);
%     delta2 = Delta2 -2*ramp;

g = [];
for n = 1:length(tvec)
    t = tvec(n);
    if t>=0 & t<= ramp
        g = [g t/ramp];
    end
    
    t1 = ramp;
    if t>t1 & t< t1 + delta1
        g = [g 1];
    end
    
    t1 = ramp+delta1;
    if t>=t1 & t< t1 + ramp
        g = [g 1-(t-t1)/ramp];
    end
    
    t1 = 2*ramp + delta1;
    if t>=t1 & t< t1 + ramp
        g = [g -(t-t1)/ramp];
    end
    
    t1 = 3*ramp + delta1;
    if t>=t1 & t< t1 + delta2
        g = [g -1];
    end
    
    t1 = 3*ramp + delta1 + delta2;
    if t>=t1 & t< t1 + ramp
        g = [g -1+(t-t1)/ramp];
    end
    
    t1 = 4*ramp + delta1 + delta2;
    if t>=t1
        g = [g 0];
    end
end
end


function g = fMCWaveform(tvec, n, intervals, ramp)
% generates waveform with all moments up to n nulled
% if ramp exists then make trapezoidal else sinusoidal lobes

% set n = -1,0,1,2... (if n = -1, then no moment is nulled)
% length(intervals) must be equal to n+1

% intervals (required for n>-1) are given in fraction of max. time
% e.g. [0.5] splits the wave in 2 intervals
% e.g. [0.1 0.9] splits the wave in 3 intervals 0-0.1*T, 0.1*T-0.9*T, 0.9*T-1*T

%solve:

% m(0,1)     m(0,2) .. m(0,p)     C(2)    0
%   .           .         .
%   .     +                    *   .  =   .
%   .           .         .
% m(n,1)     m(n,2) .. m(n,p)     C(p)    0

% where n is max order of nulled gradient moment (e.g. velocity, n = 1)
% p is n+2

% rewritten as

% m(1,1)      M(1,1)   .. M(1,p)       C(1)       0
%   .           .          .
%   .      +                       *    .     =   .
%   .           .          .
% m(n+1,1)    M(n+1,1) .. M(n+1,p)     C(p-1)     0

% i.e. M*C = - m, C = -m\M


T = max(tvec);

if nargin<3
    if n>-1
        error('need intervals')
    else
        intervals = [];
    end
end

if length(intervals) ~= n+1
    error('length(intervals) should be equal to the max. grad. nulling order + 1')
end

intervals = [0 intervals 1];


if n>-1
    % generate waveforms and moments
    for i = 1:n+2
        t = tvec(tvec>=T*intervals(i) & tvec< T*intervals(i+1));
        
        if exist('ramp','var')
            s = sprintf('g%d = ftrapez(t,T*intervals(i),T*intervals(i+1),ramp*T);', i);
        else
            s = sprintf('g%d = fsin(t,T*intervals(i),T*intervals(i+1));', i);
        end
        eval(s);
        
        for j = 0:n % moments
            s = sprintf('m(%d+1,%d) = sum(g%d.*t.^%d);', j, i, i, j);
            eval(s);
        end
    end
    
    M = m(:,2:end);
    m = -m(:,1);
    C = inv(M)*m; % solve for coeficients C
    
    g = g1; % compose final waveform
    for i = 2:n+2
        s = sprintf('g = [g C(%d)*g%d];', i-1, i);
        eval(s);
    end
    g = [g 0];
    
else
    if exist('ramp','var')
        g = ftrapez(tvec,0,T,ramp*T);
    else
        g = fsin(tvec,0,T);
    end
    
end
end



function g_lobe = fsin(x,t1,t2)
g_lobe = sin(pi*(x-t1)/(t2-t1));
end

function g_lobe = ftrapez(x,t1,t2,ramp)
N1 = find(t1<=x & x<=t1+ramp);
N2 = find(t1+ramp<x & x<t2-ramp);
N3 = find(t2-ramp<=x & x<=t2);

rampUP = (x-t1)/ramp;
rampDown = 1-(x-t2+ramp)/ramp;
g_lobe = [rampUP(N1) ones(1,numel(N2)) rampDown(N3)];
end

function mat = frepvec(vec,dims)
%      if size(vec,1) == 1
%          vec = vec';
%      end
l = length(vec);

N = prod(dims);


for n = 1: length(dims)
    if n == 1
        s = ['[' num2str(dims(n))];
    else
        s = [s ',' num2str(dims(n))];
    end
end
s = [s ']);'];

s1 = [sprintf('mat = reshape(repmat(vec,%d,1),', N/l) s];
eval([s1])



end





function gij = fouter2(g)
tim = size(g,1);
dir = size(g,2);

gij = zeros(tim,dir,dir);

for i = 1:dir
    for j = 1:dir
        gij(:,i,j) = g(:,i).*g(:,j);
    end
end

end

function res = fVectorEulerRot(v, alpha, beta, gamma);
% rotates a vector arround the Euler angles alpha, beta, gamma (ZYZ rotations)
% angles or vector components can be vectors but not both
% size(v) = 3 x N
% ZYZ*v
if size(v,1) > size(v,2)
    v = v';
end

v1 = v(1,:);
v2 = v(2,:);
v3 = v(3,:);

N = length(alpha);
sA = sin(alpha);
cA = cos(alpha);
sB = sin(beta);
cB = cos(beta);
sG = sin(gamma);
cG = cos(gamma);

res(:,1) = conj(v2).*(cA.*sG + cB.*cG.*sA) - conj(v1).*(sA.*sG - cA.*cB.*cG) - cG.*sB.*conj(v3);
res(:,2) = conj(v2).*(cA.*cG - cB.*sA.*sG) - conj(v1).*(cG.*sA + cA.*cB.*sG) + sB.*sG.*conj(v3);
res(:,3) = cB.*conj(v3) + cA.*sB.*conj(v1) + sA.*sB.*conj(v2);
end


function [f PS] = PowerSpectrum(F,dt,N0,thresh,fmax)
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
end


function [f PS] = fnormPS(F,dt,mode)
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


end



function projection = fPSTprojection(PS,n)
% calculate projection of power spectrum tensor along (n.x,n.y,n.z)
% PS(frq,m,n): frq channels, m,n = 1,2,3

projection = PS(:,1,1)'.*n.x.^2 + PS(:,2,2)'.*n.y.^2 + PS(:,3,3)'.*n.z.^2 ...
    + PS(:,1,2)'.*n.x.*n.y + PS(:,1,3)'.*n.x.*n.z + PS(:,2,3)'.*n.y.*n.z ...
    + PS(:,2,1)'.*n.x.*n.y + PS(:,3,1)'.*n.x.*n.z + PS(:,3,2)'.*n.y.*n.z;

end


function vout = fVectorRot(v, Ax, Ay, Az)
%rotate vectors v(n,3) around ZYX

v1 = v(:,1);
v2 = v(:,2);
v3 = v(:,3);

sAx = sin(Ax);
cAx = cos(Ax);
sAy = sin(Ay);
cAy = cos(Ay);
sAz = sin(Az);
cAz = cos(Az);

vout(:,1) = conj(v3).*(sAx.*sAz + cAx.*cAz.*sAy) - conj(v2).*(cAx.*sAz - cAz.*sAx.*sAy) + cAy.*cAz.*conj(v1);
vout(:,2) = conj(v2).*(cAx.*cAz + sAx.*sAy.*sAz) - conj(v3).*(cAz.*sAx - cAx.*sAy.*sAz) + cAy.*sAz.*conj(v1);
vout(:,3) = cAy.*sAx.*conj(v2) - sAy.*conj(v1) + cAx.*cAy.*conj(v3);

end

function [V, L] = fB(q, dt)
for i = 1:3
    for j = 1:3
        B(i,j) = dt*sum(q(:,i).*q(:,j));
    end
end

[V L] = eig(B);
end

function gout = fscaleG(g, b,dt)
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


function A = findRotation1(g,h,N)

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
end


function A = findRotation(g, h, Ax0, Ay0, Az0, dAx, dAy, dAz, N, Nlevels)
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
end

function R = fRotXYZ(phiX,phiY,phiZ)
% Rotation matrix around XYZ

Rx = [1 0 0
    0 cos(phiX) -sin(phiX)
    0 sin(phiX) cos(phiX)];

Ry = [
    cos(phiY) 0 sin(phiY)
    0 1 0
    -sin(phiY) 0 cos(phiY)];

Rz = [
    cos(phiZ) -sin(phiZ) 0
    sin(phiZ) cos(phiZ) 0
    0 0 1];


R = Rz*Ry*Rx;


end



function [px,py] = patchXY(x,y1,y2,n)
l = length(x);
j = 1;
for i = 1:n:(l-n)
    
    px(1,j) = x(i);
    py(1,j) = y1(i);
    
    px(2,j) = x(i);
    py(2,j) = y2(i);
    
    px(3,j) = x(i+n);
    py(3,j) = y2(i+n);
    
    px(4,j) = x(i+n);
    py(4,j) = y1(i+n);
    
    j = j + 1;
end

end


function UDSR = fUDSR1000()
% electrostatic repusion - uniform orientation by Malcolm Levitt

par = [ -0.109126926451332  -0.286205517025812  -0.951933672031446  -1.935068937515610   2.830284436109598
    0.398992412811417   0.074033366179763  -0.913960674871199   0.183464298019276   2.723735656628742
    0.498509199332389   0.118526770730422  -0.858743257790824   0.233428226687057   2.603608306098054
    0.582686503939233   0.168924953445543  -0.794947041148336   0.282171724876262   2.489716608243222
    0.667745812933608   0.107572081582141  -0.736575710007953   0.159725038629462   2.398789770127735
    0.803489420467277   0.096707167404686  -0.587411674185764   0.119782795433759   2.198653166995878
    -0.084686219298883   0.544929492484748  -0.834194277420540   1.724970701675982   2.557466450063248
    0.232826055566875  -0.886993142697387   0.398792167246278  -1.314098168589806   1.160596955915261
    0.863902181931041   0.025591071656991  -0.503009062648213   0.029613989558686   2.097873167747283
    0.949612590046770   0.052279168033649  -0.309035301243038   0.054997635513283   1.884974841194313
    0.931812759204895   0.171331961207177  -0.319953654224894   0.181838472204446   1.896476896489660
    0.248969125520356  -0.260044352080549   0.932947645631113  -0.807152949466738   0.368280830842357
    0.942433498290561   0.245450167186556  -0.227097592959296   0.254782923032574   1.799892693248315
    0.991285989647086   0.077409513627156  -0.106582615512121   0.077931836734830   1.677581774717299
    0.972939081901837   0.197123526673431  -0.120548156932572   0.199900296805373   1.691638374449277
    0.180756822975933  -0.631380979781237   0.754112080077048  -1.291966296466439   0.716495258297370
    0.990423781845750   0.094030623603129   0.101088941923288   0.094656071927474   1.469534417406337
    0.973581489087901   0.211897777656026   0.085078880633349   0.214305501662101   1.485614471000505
    0.645907669021471  -0.020951263473736   0.763127989041225  -0.032425563473124   0.702656712347281
    0.967176956236385   0.159994597360157   0.197409382099442   0.163939671804406   1.372081734279568
    -0.168375581601341  -0.580411336804408   0.796726015408892  -1.853143212105257   0.648938069585121
    0.945017531919741   0.110285994720804   0.307861760751093   0.116177068877096   1.257851506641806
    0.898846310454690   0.181294078299397   0.399008480305324   0.199026194539388   1.160361062342419
    0.908348185559608   0.064862053669236   0.413154314735278   0.071285602795964   1.144881211197398
    0.855234139002418   0.135127334103697   0.500315071792586   0.156704966049306   1.046833699402882
    0.858451700100438   0.018250326749303   0.512569609095400   0.021256382810444   1.032621570069048
    0.799214258639749   0.084360059148499   0.595096588132867   0.105164333967747   0.933410497254184
    0.976781232592675  -0.212140304573821   0.029915127110897  -0.213861796184981   1.540876735971028
    0.655997062790186   0.092551433442980   0.749067477453332   0.140160056657861   0.724142964754867
    0.476506588214086  -0.002548613153390   0.879167205916808  -0.005348486073953   0.496684636424801
    -0.868390339483484   0.037699460466646   0.494446123427298   3.098206860879761   1.053598830361806
    0.064562958741558  -0.011267986364851   0.997850017207906  -0.172786787986983   0.065585877114347
    -0.055051692247717  -0.739206363507604  -0.671225195691079  -1.645133146726597   2.306656749952357
    -0.133724258696708   0.988372243948104   0.072374926797925   1.705277183569472   1.498358065714590
    0.864771390556963  -0.363335397604217   0.346637895969213  -0.397757161124098   1.216811944275040
    0.219343357516187  -0.804938886796531  -0.551327380089441  -1.304758706886153   2.154750761329024
    0.782873701792079  -0.417401558275636  -0.461394306630932  -0.489826959012310   2.050362473733570
    0.486651407237353   0.241093014217386  -0.839669319630936   0.459970588236760   2.567470384366406
    0.775627717990564   0.342038990457061  -0.530481830124147   0.415330515532470   2.129965190575951
    0.327875761120637  -0.276086930160680  -0.903478551081327  -0.699859010054577   2.698613126870505
    0.805919969652866   0.230916449957243  -0.545133557629566   0.279049291173331   2.147344743539092
    0.882066492248590   0.212758165507889  -0.420348267823226   0.236683309433469   2.004625436961395
    0.898319304584050   0.288018002607137  -0.331765063238769   0.310264109068104   1.908970322691929
    0.923268981041298   0.382514946115607  -0.035450030255795   0.392777117683074   1.606253786289102
    0.902648730881007   0.360150164821185  -0.235620727904359   0.379637628034087   1.808653556237283
    0.939861158107303   0.315412234972832  -0.131057718241397   0.323784221620796   1.702232152179483
    0.964066800574132   0.264463234761591  -0.025187327968485   0.267734325315367   1.595986318670346
    0.853928622829264   0.393969242049776  -0.339991387290387   0.432261351004627   1.917704066025155
    0.941834627248505   0.328917515672955   0.068998570993909   0.335989262072525   1.501742890079002
    0.906816413371741   0.389238097550962   0.161795228157443   0.405452957469626   1.408286746804960
    0.929652800955955   0.224918333278670   0.291817430990496   0.237376617791185   1.274669902117137
    0.858427344325716   0.446199933995375   0.252998247852472   0.479352127213548   1.315018250730219
    0.901534327023955   0.336552511422191   0.271971072456024   0.357288916279520   1.295355604869954
    0.876374134993293   0.295994756509222   0.379941416054903   0.325719651333859   1.181063364416722
    0.810391888067436   0.366288245980186   0.457272247803488   0.424506775991924   1.095870765848341
    0.584355897277603   0.050612121296773  -0.809917649205602   0.086396202621489   2.514808028213348
    0.715367634546124  -0.197657183705779  -0.670209508417617  -0.269575845171639   2.305287369217792
    0.650101378644991   0.211120825950015   0.729928896766024   0.314006327110894   0.752578406569310
    -0.788323478296067   0.217101791460273   0.575684727704263   2.872857545101817   0.957355013894313
    0.572904120060350   0.152332072246981   0.805342044713187   0.259881379617234   0.634544096145002
    0.562085685494939   0.046223017151478   0.825786361504656   0.082050199253656   0.599201246897065
    0.330338297180433   0.152559233506103   0.931451710872852   0.432645839142913   0.372413943410999
    0.140726196304790   0.085328322893086   0.986364646054307   0.545070312266337   0.165326672139195
    -0.011678952917332   0.997241905975050   0.073295177372083   1.582507045082023   1.497435364413002
    -0.426940986244832   0.885950541462084   0.181143678739644   2.019860671253062   1.388647082511131
    0.396664970295136   0.338838906500002  -0.853138380793317   0.706938539332231   2.592768209268893
    0.580079424063331   0.279108377867676  -0.765249224228570   0.448458638785269   2.442224722160040
    0.643312715183393   0.336687178107658  -0.687597625491299   0.482164193614628   2.328971537009603
    0.729479897983468   0.148974054572122   0.667582061998652   0.201449505117858   0.839839863267667
    0.692373552453575   0.392395976738671  -0.605511569915986   0.515604721834467   2.221204844954017
    0.725012466061034   0.451753327867823  -0.519880615926815   0.557229929134539   2.117507516649729
    0.728333904315397   0.528346903096478  -0.436324734358540   0.627582464780769   2.022306350713630
    0.796148888977731   0.421812859774816  -0.433844278408873   0.487215351112134   2.019551479687231
    0.794381259080798   0.501989652153245  -0.342001175951589   0.563563846479768   1.919841992610008
    0.848759954263181   0.469258030746971  -0.243728210551503   0.505047708217003   1.817004468195247
    0.892021643732243   0.429680627825934  -0.140256711726531   0.448895083425833   1.711517011357845
    0.868184490153991   0.494286537787704  -0.044005790731445   0.517565347698036   1.614816332856984
    0.799030410684171   0.599428708909943  -0.047282404089087   0.643625969159678   1.618096366258071
    0.830463686606994   0.537602806527366  -0.145990710803873   0.574511861871029   1.717310665753137
    0.831732186046521   0.552595386697765   0.053478119774572   0.586425493428817   1.517292683730595
    0.893147163376584   0.445684244325108   0.060445834535714   0.462850704239759   1.510313623128052
    -0.330064190556188   0.868552349199598   0.369695072745564   1.933957721242985   1.192115505142071
    0.475261178286047  -0.399789031915454  -0.783770082597098  -0.699362146658561   2.471509591200379
    0.840927418907709   0.409239431257101   0.354068021760766   0.452912706032102   1.208878981787526
    0.736671871071982   0.435832129327586   0.517073408151189   0.534237279210125   1.027368072399781
    0.767647562515012   0.320320660998930   0.555078277274449   0.395309932688207   0.982339257515571
    0.714616307273394   0.268192632148552   0.646062106488354   0.359030140202717   0.868382367331994
    0.688289160339205   0.385932822665018   0.614258812024186   0.511031097859569   0.909349981251539
    0.003580667975678  -0.884950422982122   0.465671480424355  -1.566750169696713   1.086403095683500
    0.561829572937730   0.272111104962740   0.781218969002011   0.451033477757566   0.674180216896348
    0.487353921654565   0.206949713232673   0.848326571103859   0.401565239550151   0.557979634977692
    0.278441436937335   0.267665243335188   0.922402126898343   0.765667961162171   0.396541912329653
    0.254014262398318   0.064148422251648   0.965070844249718   0.247366562256075   0.265082674941274
    0.159822054201367   0.178078561491535  -0.970950532688438   0.839374903767417   2.899968100361401
    0.332072113387690   0.265842536622406  -0.905017048034123   0.675080267077234   2.702216124055907
    0.776761053323799  -0.243683907337421  -0.580741267126825  -0.303994100606879   2.190435270985874
    0.954475919990771  -0.264931550820198   0.137065646818508  -0.270751679229811   1.433297835374960
    0.559636978104050   0.407152320328200  -0.721826461685861   0.628964103365413   2.377234139469321
    0.613028062117535   0.466486865219817  -0.637641435003881   0.650475620255233   2.262228951868315
    0.548267080040894   0.559273942375618  -0.621784421100325   0.795335944270153   2.241815380840951
    0.649969105200895   0.530853319541682  -0.543815148203804   0.684862978281745   2.145772920561529
    0.655269541472623   0.614244924406582  -0.439687389924946   0.753094202459215   2.026046910881071
    0.094870874227704   0.463549843904683  -0.880977331966728   1.368922459302384   2.648720120907349
    0.904440643059410  -0.144764617662664  -0.401285844075348  -0.158713600477841   1.983716574017192
    0.722233780904637   0.601955733172041  -0.340628332675812   0.694815285667243   1.918381441711998
    0.702121294889152   0.670461726123755  -0.239805673549095   0.762336566756242   1.812962005737638
    0.755495620646728   0.638868196710060  -0.145168159098514   0.701950368995138   1.716479256660851
    -0.541845821973200   0.813429351745500  -0.211508380280976   2.158424002669130   1.783914325054281
    0.402683706215550   0.744334006008082  -0.532740762705905   1.074890183515862   2.132632209294095
    0.757739224319304   0.650395186004938   0.053078903055558   0.709313930163221   1.517692468269681
    0.707688921535098   0.690331388935552   0.150396023174822   0.772983007013181   1.419827486779667
    0.787628205518276   0.595512307536691   0.158135705785842   0.647382884008102   1.411994011324231
    0.671083663052920   0.661188300807507   0.335375532886281   0.777970865852524   1.228792505860750
    0.731584097949161   0.572518576397365   0.370144819371887   0.664024545881860   1.191631419211896
    0.693070292782523   0.539193262955852   0.478460232877685   0.661167751244616   1.071895956739253
    0.776582414318511   0.478989063908159   0.409254481255389   0.552674442428393   1.149159492705393
    0.597831887052948   0.596754233178128   0.535239591217518   0.784496048795383   1.006004964635376
    0.648740199471992   0.494304134327540   0.578618679594551   0.651101968470605   0.953762284191843
    0.908167250803775   0.090557476045579  -0.408694981740378   0.099386006909504   1.991820041427242
    0.630443834306427   0.328845343169367   0.703136766255919   0.480785193953904   0.790996968338279
    0.529980740023778   0.381827134889652   0.757184557598613   0.624322473004122   0.711804270811953
    0.457495323682000   0.319287845370821   0.829911622166819   0.609317909973771   0.591847074152740
    0.388335334978777   0.236201342991871   0.890732615982910   0.546456791187176   0.471841875953378
    0.213981187477109   0.178009561825177   0.960481466403540   0.693887944404322   0.282069486645473
    -0.935044635253515   0.261486809387890  -0.239407975221718   2.868906970525220   1.812552374862275
    -0.812463741654531   0.220472374343626  -0.539717148743676   2.876611123519474   2.140897410852688
    0.273217946783699   0.364081935564240  -0.890391092582971   0.927020180439326   2.668999949075605
    0.330602544043835   0.443751015283741  -0.832938769842293   0.930491257539970   2.555193701496798
    0.375176330488646   0.529748089328766  -0.760663975020242   0.954577201603519   2.435131673899703
    0.522503292613081  -0.295601011775369  -0.799756432325409  -0.514847297456474   2.497685708463135
    0.473330703388525   0.647189108574599  -0.597582047063113   0.939329443308397   2.211278407776656
    0.581322962862643   0.619363254168346  -0.527686244121064   0.817069694800930   2.126670717978975
    0.503433010278841   0.700851778233053  -0.505333542434254   0.947883255412248   2.100564774195579
    0.575994048482665   0.700333225339632  -0.421620954884761   0.882511227330493   2.006028515388218
    0.642064130138761   0.690205674090895  -0.333721111471674   0.821517429032623   1.911044566996819
    0.606145053397993   0.757490337523042  -0.242480025569884   0.895932410815559   1.815717679714063
    0.665355667466157   0.732192999600776  -0.145585875366679   0.833186442059356   1.716901458262704
    0.815246680446633  -0.465897002443231   0.343966325583126  -0.519184724968099   1.219658609819469
    0.862029679289900  -0.503483086316175  -0.058391898555293  -0.528621670541374   1.629221458670255
    0.621287433305034   0.782039739244177  -0.049150498066078   0.899451869870431   1.619966635800036
    0.781602299635899   0.573754878281563  -0.244751271400136   0.633231103617231   1.818059480267324
    0.669556144070919   0.741091222849541   0.049783223602483   0.836065577053527   1.520992516688769
    0.613262208431460   0.776007066755677   0.147385535433494   0.902010899758980   1.422871909977319
    0.348516487996039   0.927011608364551  -0.138512582649647   1.211186921331686   1.709755688563486
    0.386110888555020   0.848176966814745   0.362648886258033   1.143606253075964   1.199687620969944
    0.493266245297553   0.791496725727125   0.360861946474401   1.013487448633424   1.201604377227524
    0.583219966747173   0.728452556574190   0.359459793590261   0.895672195840762   1.203107397687015
    0.635388299778859   0.635592290629813   0.438524969182004   0.785558662189280   1.116839571039247
    0.542401790215854   0.700479783407478   0.463819330136261   0.911905623765189   1.088494871289125
    0.506820837617601   0.630228582135578   0.588170530387210   0.893508311318374   0.942001483726986
    0.444760632259627   0.581393822676217   0.681299642555280   0.917770254331748   0.821259699023147
    0.489725703091184   0.484477308079328   0.724879627033380   0.780010842503898   0.759936686910266
    0.419498960864366   0.424615373486877   0.802323130934362   0.791459368329302   0.639619167855210
    0.348597020717518   0.349703571555561   0.869590437613107   0.786982796038579   0.516424066665703
    0.235298118711668   0.378968601155128   0.894995862934168   1.015157399037718   0.462374398490144
    0.041169822647704   0.164663279997973  -0.985490258664928   1.325794785605788   2.971034941920877
    0.096728152051601   0.273590152209106  -0.956970267675481   1.230961281914827   2.847171472224329
    0.211854001940503   0.461714300957323  -0.861358105640902   1.140613608519123   2.608733415829919
    0.495909934002981   0.346128338247297  -0.796409763136618   0.609359647347286   2.492131463050868
    0.006628564940044  -0.267006560596241  -0.963671914463321  -1.545975949007255   2.871221810037199
    0.365399038196405   0.686983907060590  -0.628121528308741   1.081964320506051   2.249933048707552
    0.298274570288285   0.772864651494020  -0.560100447411325   1.202474508371391   2.165303372932480
    0.357232157437771   0.818133294277277  -0.450603038701878   1.159105448616860   2.038237054527236
    0.052470601621499   0.900835321983158  -0.430978605771769   1.512615455765250   2.016373316601404
    0.314830946447731   0.877258078257766  -0.362353058894661   1.226232028426190   1.941587615627883
    0.538859711995486   0.772503668598732  -0.335958766502145   0.961716744342386   1.913419304420078
    0.499304745679922   0.830879106621625  -0.245631189227179   1.029689247803343   1.818967103572861
    0.778376324735086  -0.149266545701595   0.609794879796652  -0.189466396320889   0.914994569043169
    0.564498524858988   0.812072417441225  -0.147918235059753   0.963344122529851   1.719259346661990
    0.409473603723613   0.911440658582995  -0.040091068025113   1.148558162481544   1.610898142282354
    0.517114683975459   0.854539086403215  -0.048531983537865   1.026606924772570   1.619347382216430
    0.613205048351666  -0.263930622989480   0.744526826195279  -0.406445455936958   0.730970542453506
    0.569497232456622   0.820475805935133   0.049923482444055   0.964036986683811   1.520852083224880
    0.513610559065566   0.843452363191838   0.157455722799966   1.023830359834076   1.412682621384203
    0.432204093108399   0.780516311746803   0.451656848722029   1.065085278773888   1.102174804220959
    -0.374360877548754   0.402635279508123  -0.835307587093726   2.319821167858127   2.559488344600970
    0.737842101756706   0.626644227161175   0.250810776165659   0.704083545916218   1.317278614989921
    -0.162366551918535   0.538862018283273   0.826598347487922   1.863457988268510   0.597759990943955
    0.458388655157131   0.710877729731819   0.533416061057944   0.998070731535610   1.008162281364595
    0.352575120637282   0.756389588196095   0.550967853124070   1.134610774406628   0.987272769093872
    0.394657870502473   0.671462204843843   0.627203055410844   1.039426279147465   0.892839424894177
    0.329762795663857   0.619631482995137   0.712259309435268   1.081728157035365   0.778084578727817
    0.376803186185095   0.524934460870323   0.763192748047008   0.948217103948933   0.702556498651788
    -0.115824399691132  -0.876213950190708  -0.467796774173769  -1.702221723876638   2.057592656352649
    0.306669828451638   0.456198355375875   0.835366193276810   0.978952659727110   0.581997701580514
    0.165232405180870   0.293265984722713   0.941643942519020   1.057707015933288   0.343315104740737
    0.098684391992089   0.201705612280730   0.974462024273705   1.115785958008505   0.226483642410941
    -0.021937644885683   0.261057306948949  -0.965073998315897   1.654633182124139   2.876522017838876
    0.157515161494783   0.366834453463136  -0.916854109251671   1.165212927904302   2.730924057598337
    0.149213883470105   0.555261406649525  -0.818180901308116   1.308271121492315   2.529036320577189
    0.203571563604680   0.636803491090983  -0.743666546394204   1.261386215923308   2.409334405387938
    0.141165312807714   0.719222045861744  -0.680288176588677   1.376985042256835   2.318952066168102
    0.190307461881617   0.787064227447879  -0.586781877552512   1.333555679598203   2.197875176088683
    0.958906176745552  -0.275174020138041  -0.069124545859610  -0.279457051409754   1.639976039873989
    0.789996327989390   0.202702693223248   0.578634098479619   0.251168407477226   0.953743379063716
    0.250069399307263   0.841745310464786  -0.478466433347885   1.282016537813446   2.069703758024678
    0.254145644466183   0.710237720913209  -0.656484859833705   1.227161508658236   2.286945689982450
    0.172715066925053   0.891008200819014  -0.419849844266182   1.379328595562049   2.004076196729305
    0.101565706238827   0.938475627280218  -0.330072574318450   1.462991777976469   1.907176784109685
    0.285696192741324   0.928533346155579  -0.237072795845992   1.272303667401253   1.810147962703989
    0.231100673487963   0.963716489043904  -0.133577720665911   1.335439030053512   1.704774509362472
    0.760675355185381  -0.648085109823435   0.036860472578914  -0.705646181888391   1.533927502088529
    0.798590938630267  -0.325888518067038   0.506013029999919  -0.387451922343940   1.040240273198986
    0.290648524361478   0.941962056912404   0.168020589880672   1.271508127065484   1.401974958996249
    -0.853212392887771   0.164104091246519   0.495074196316954   2.951576352856996   1.052876101306111
    -0.928222512431221   0.342802952984822  -0.144530629420750   2.787818744338602   1.715834931736694
    -0.685341479661646   0.668858104973516  -0.287985922687937   2.368365931460071   1.862919318206141
    0.341244152148586   0.901378665356713   0.266587565846840   1.208890515858003   1.300945602801503
    0.274625099242015   0.887554994924438   0.369901589413988   1.270722120822828   1.191893231983626
    0.321150033609078   0.826171622803646   0.462928834257361   1.200051117442142   1.089499761106668
    0.244186485700027   0.790412471512431   0.561801642111414   1.271162103716661   0.974234322007691
    -0.590848903947566  -0.116927821993206   0.798264027216487  -2.946219053053174   0.646388841449870
    0.284640424481120   0.710165650651680   0.643929015797304   1.189593167008478   0.871173692978471
    0.554416930006989   0.791693891205731   0.256598227486383   0.959874389398264   1.311295392825772
    0.261573423985616   0.556542738496506   0.788568021220269   1.131437887120597   0.662319438050467
    0.146071461608617   0.583854369712558   0.798609543563733   1.325643785513458   0.645814968954140
    0.592167669455054  -0.682111891630964   0.429023097918745  -0.855865663526222   1.127385317055486
    0.190944592536984   0.484372705577226   0.853770018607285   1.195291222259565   0.547612477461656
    0.050417834945679   0.316827708185830   0.947142146272197   1.412986256652231   0.326589007244224
    -0.199313950682787   0.335615213306444  -0.920671699174298   2.106701004663825   2.740594149673956
    -0.824285934930524   0.547705664590031   0.143426644845731   2.555118766525442   1.426873330240240
    0.160109870732247   0.913500823848855   0.374006783523463   1.397288077986805   1.187470722688838
    0.024800079885128   0.718966205937604  -0.694602440794285   1.536315915004181   2.338663449207491
    -0.111475652169225   0.842835282127279  -0.526499635493582   1.702295758730163   2.125274391595364
    0.088046744275414   0.641746061305035  -0.761846286085306   1.434448856607123   2.436954956147843
    0.074010974020271   0.789781721138294  -0.608906568104175   1.477358533394469   2.225477752962228
    0.446936682055207  -0.233413172630509  -0.863577381058607  -0.481289618397056   2.613118398864625
    0.125889831774278   0.847647629792635  -0.515407844296895   1.423357350541943   2.112279836165043
    -0.541815293291902   0.276627469864836  -0.793670858020621   2.669535333873189   2.487615941411334
    0.042776189741646   0.971831081368749  -0.231763989603293   1.526808644620093   1.804686982889961
    0.164006405485073   0.959480452260199  -0.229126953216837   1.401499988559077   1.801977007216510
    0.107087829401447   0.985742606049362  -0.129783324865327   1.462583994459668   1.700946780048850
    0.048269031935503   0.998420533147463  -0.028749600789319   1.522488547705427   1.599549889504995
    0.172740444569441   0.984460938479548  -0.031581631020558   1.397097536659323   1.602383210089906
    0.112427318172078   0.991136784797127   0.070767025871153   1.457846420390125   1.499970100867199
    0.979150795183713   0.007286907434390  -0.203003993239406   0.007441931260985   1.775221150138437
    0.172209875897506   0.970004387182467   0.171567034975058   1.395091904125983   1.398376257264243
    0.051200675682682   0.983601462993768   0.172935400673653   1.518788976908799   1.396987128207667
    0.108834610998795   0.955519384432069   0.274112629083924   1.457384092832825   1.293129458946988
    0.209282052031154   0.857273299473708   0.470407815312513   1.331354444266327   1.081043465373513
    0.172879744780004   0.738149837843583   0.652109968284715   1.340736144269968   0.860432067300442
    0.096041173100935   0.876312678457872   0.472072221851049   1.461635075574510   1.079156393338335
    0.019145854004747   0.826070020024370   0.563242184403269   1.547623440210526   0.972491973877031
    0.132932594041121   0.813420486151480   0.566282648639872   1.408804128841128   0.968807689947317
    -0.016066695413209   0.690340360332012   0.723306330813693   1.594065711556998   0.762217812449829
    0.028782828309228   0.604314293536502   0.796225962539529   1.523203388059907   0.649765054267335
    0.074363240524602   0.508353333045691   0.857931813864025   1.425543958886008   0.539565847665131
    0.119866625182807   0.403978271233983   0.906881220743977   1.282355528625028   0.434973472976271
    -0.557074531976086  -0.824025145861235   0.103201379893848  -2.165260595437816   1.467410871829680
    -0.138774912206184   0.243696872179351  -0.959871532150093   2.088454986869109   2.857340090923571
    0.222240280229685   0.263631112934168  -0.938673475781927   0.870381845546203   2.789559116917963
    -0.024785356164463   0.459178712867984  -0.887998083201029   1.624721576587630   2.663769587981490
    -0.028451369376762   0.636642695006513  -0.770633893931114   1.615456308811001   2.450631570741951
    -0.922704966160878   0.340759194048832   0.180274005594937   2.787824530801756   1.389531313097030
    -0.090702858713201   0.711027083473289  -0.697290096006477   1.697677010638654   2.342406220663620
    -0.157712530521754   0.774576401251903  -0.612501556194003   1.771661894628276   2.230017698901446
    0.006034488250746   0.850219298303475  -0.526394082170286   1.563698879188355   2.125150241600366
    -0.065498692861945   0.898678663770078  -0.433689497818200   1.643551014126217   2.019379696919290
    -0.182362889784450   0.884604027884393  -0.429207980214923   1.774100212671349   2.014412022254740
    -0.133764509098523   0.933180711569345  -0.333587793032902   1.713169054644414   1.910903144399206
    -0.197304784845242   0.952744270567560  -0.230974407192357   1.775000786799442   1.803875378394001
    -0.077133913357160   0.969162664101826  -0.234038650485899   1.650217119345029   1.807025965530920
    -0.137244966504843   0.981671371438011  -0.132231379286836   1.709703410048190   1.703416117363671
    -0.196641551173160   0.979979956959885  -0.031167038824925   1.768825217812244   1.601968413688839
    -0.074777072183331   0.996766567670869  -0.029407465877283   1.645675707725248   1.600208032913940
    -0.253585374299152   0.964843483573472   0.069076118502460   1.827809056021861   1.501665157105480
    -0.072243706179097   0.981824813893235   0.175501799822744   1.644245018581020   1.394380870773679
    -0.011320801218990   0.961867669984157   0.273280853517054   1.582565387012941   1.293994255919532
    -0.191220925879054   0.965894699927408   0.174590910869094   1.766241945601918   1.395306044263089
    -0.134840061415390   0.951012391614645   0.278197032394106   1.711643317213630   1.288879795839767
    -0.079424240234137   0.926651138043069   0.367436332480493   1.656298406438227   1.194545299622706
    -0.018442031263478   0.884002557087993   0.467118154801081   1.591655267028870   1.084767652729097
    -0.016022777552478   0.942107395904423  -0.334928238853773   1.587802066291642   1.912325396002945
    -0.097563906651182   0.829476113470143   0.549955144808580   1.687879483211499   0.988485796177741
    -0.058630954122996   0.765697430523998   0.640523111298550   1.647219158606389   0.875617065817375
    0.101547237162746   0.675079254700691   0.730723038159721   1.421493016511572   0.751415849143240
    -0.090536829607366   0.618643367277305   0.780437996645758   1.716112119913105   0.675430282250615
    -0.043742092731938   0.527247226222976   0.848585288444013   1.653569920889019   0.557490812068711
    0.003410491456571   0.426262970283187   0.904592863510087   1.562795587994140   0.440372674752371
    -0.016218335531126   0.224149762161399   0.974419750269661   1.643025365502283   0.226671824578560
    -0.075327945334659   0.147917117686888  -0.986126881768701   2.041822713643202   2.974827461372236
    0.032871078534876   0.553886817930412  -0.831942837651048   1.511519650620466   2.553396357881437
    -0.141973947484438   0.443289617219414  -0.885063677653281   1.880747408622542   2.657427122982951
    -0.200087125483125   0.526337432989223  -0.826398238653868   1.934074320813901   2.543477191161753
    0.449835397663214  -0.511167082053303  -0.732363522599460  -0.849132252940872   2.392582932666979
    -0.258282925871505   0.601538319974936  -0.755937550201726   1.976363159047760   2.427881417851807
    -0.341483704173637   0.792097891296571  -0.505934591016841   1.977833331713624   2.101261441791847
    -0.489345431600086   0.863910001172271  -0.119166096045197   2.086166953569411   1.690246278032369
    -0.227202653591721   0.824128144332912  -0.518836926133127   1.839802519816139   2.116286192867415
    -0.298167517286644   0.857071812039832  -0.420147641480798   1.905591015580000   2.004404339006152
    -0.250912257683936   0.910857535307824  -0.327691301863363   1.839597738791774   1.904655239926529
    0.854000722996968   0.143169816490252  -0.500185134492023   0.166101403540452   2.094608890486433
    -0.259473537796388   0.957140347341682  -0.128669493952805   1.835526118308064   1.699823530757989
    -0.565771353690634  -0.019867378866769  -0.824322790295180  -3.106491514938836   2.539801253005122
    -0.319300361416728   0.947285583488640  -0.026406486205753   1.895904974599387   1.597205882848892
    0.367654790932123   0.046729579816951   0.928787543561383   0.126423905249649   0.379668503027320
    0.226410247432331   0.935321705388482   0.271867076503966   1.333298114298471   1.295463672808843
    -0.305233580563776   0.937106005506318   0.169306809491551   1.885678736521064   1.400670044634238
    -0.351168447865456   0.898079175437967   0.264829220193745   1.943538911070110   1.302769512001448
    -0.242739375964617   0.928274526646035   0.281751661106704   1.826564442500017   1.285177083799904
    -0.219690892734632   0.895367090784539   0.387367634669562   1.811407002912099   1.173021744059393
    -0.134439924106635   0.881782400940304   0.452090371715827   1.722095015069584   1.101688836545176
    0.140030332078668  -0.244730582777614   0.959426103434375  -1.051082742370253   0.285836589691261
    -0.215761222699343   0.832192589066931   0.510786246370595   1.824478967280608   1.034697233970483
    -0.178755356452169   0.770777848165173   0.611512903639427   1.798683245467889   0.912825065360713
    -0.254287656281472   0.699945974333892   0.667393003317127   1.919266789883512   0.840093762486275
    -0.135829359672683   0.700207495296411   0.700899314154062   1.762401061061334   0.794138759985244
    -0.209080588740692   0.625395561682616   0.751775032067128   1.893432237270484   0.720046552885497
    0.424944381168807   0.837884657501606  -0.342595349700600   1.101434481626192   1.920474367274568
    -0.113943483689140   0.443049429265539   0.889226678497819   1.822521204069809   0.475144387282849
    -0.130736168440852   0.243955742684296   0.960933738544317   2.062750703491996   0.280440037050724
    0.662015796669544  -0.012413959935189  -0.749387068582529  -0.018749560464891   2.417932226774313
    -0.082588022172016   0.354454353701550  -0.931418987210230   1.799712536687610   2.769088786992526
    -0.257855764739331   0.421997471012440  -0.869153921379752   2.119291243453292   2.624285232850056
    -0.792363690825302  -0.606477430842565  -0.065915911131585  -2.488306377840880   1.636760064584384
    -0.315415591656393   0.667380434179082  -0.674623124874168   2.012298864089102   2.311250351210446
    -0.425051563177504   0.626785575798290  -0.653047479599736   2.166703492050702   2.282397864381732
    -0.271048000899224   0.751768721605184  -0.601146213848712   1.916836116507260   2.215730974027105
    -0.382572549986588   0.715546782009655  -0.584492127194531   2.061784116525321   2.195050321240413
    -0.453698729571017   0.742520332284947  -0.492768727628993   2.119283269035886   2.086065082786473
    -0.413038066247264   0.813446338226857  -0.409517534003170   2.040635142502791   1.992721481702424
    -0.366863601749128   0.873248132971559  -0.320700477039471   1.968521060180837   1.897265260592045
    -0.480485623890285   0.819724445882129  -0.311745726607439   2.100973512315122   1.887826093590694
    -0.430973871791237   0.874963010145736  -0.220683603174279   2.028476164372405   1.793311624916062
    -0.315990484570361   0.921498449516484  -0.225811029845098   1.901140445309069   1.798571815997078
    -0.377917307669844   0.917564071203941  -0.123509853045113   1.961493577191793   1.694622372516693
    -0.374381107796388   0.924036851871805   0.077425335052507   1.955741503740749   1.493293425599737
    -0.434822770942181   0.900315409518279  -0.019008452175200   2.020725045078033   1.589805923849220
    -0.825182032880613   0.502820949742631   0.257382410259400   2.594329792388861   1.310483957129154
    -0.485309801892488   0.869788335552995   0.089121532319778   2.079732683188140   1.481556393979327
    0.438519530781811   0.439011019090736  -0.784200195256169   0.785958244476796   2.472202450170497
    -0.520552973784120   0.825786835014516   0.217026967451552   2.133249394991390   1.352028515962424
    -0.424142780170337   0.805541121017003   0.413770956424020   2.055443164726911   1.144203970750195
    -0.877000681511245   0.199917103823716  -0.436924428508580   2.917466937734038   2.022972941386540
    -0.401263404020980   0.753886661281046   0.520233200147793   2.059917319766599   1.023672338466292
    -0.370032166406389   0.691267748082513   0.620665043550486   2.062279654864611   0.901205720244262
    -0.477378955762826   0.674371700401052   0.563322414158168   2.186795985986518   0.972394874517967
    -0.441824721919256   0.608721467832038   0.658975788403020   2.198644490784617   0.851340065820170
    -0.327777283631065   0.622291986807749   0.710854932802963   2.055595321841919   0.780083327820578
    -0.279597055013901   0.544101103028081   0.791062245661604   2.045478007019536   0.658252828193661
    -0.231350116512517   0.452353849474448   0.861308956446694   2.043550952509929   0.532955980542481
    -0.181078540387152   0.352644553549483   0.918069921663456   2.045175557228833   0.407612414622653
    -0.085729053255347   0.128211028828240   0.988034645908101   2.160174786742323   0.154850199477522
    0.221859768034931   0.087259657904441  -0.971166306782668   0.374725835800176   2.900871519279737
    -0.316102285590289   0.317539083812478  -0.894006865351787   2.353926968169541   2.677006049440919
    -0.314262133250690   0.501649801763115  -0.805969470883204   2.130442834464250   2.508107712258621
    -0.430433898878334   0.482263699031374  -0.762986489587445   2.299467931551517   2.438717013983542
    -0.371517847668304   0.571261816079689  -0.731870498348749   2.147414755749774   2.391859151727023
    -0.531919092266907   0.453691788165181  -0.715000587855304   2.435397814308248   2.367421274156305
    -0.599277785804708   0.693898366588802   0.399213216572798   2.283155512854468   1.160137771177644
    -0.491677266825438   0.657963098933047  -0.570384103678818   2.212549097471388   2.177769739335783
    -0.557703667150158   0.674811481380207  -0.483317788049143   2.261463566890716   2.075236922224850
    -0.623725538504239   0.679813951943744  -0.385771231898162   2.313193232151319   1.966839941533394
    -0.587636399509917   0.751672310991031  -0.299453166389044   2.234323417986014   1.874915795057156
    -0.645846383848687   0.736831168912922  -0.199905670224716   2.290485893653162   1.772057973607310
    -0.540559651531847   0.841140368005302  -0.016677663135881   2.141987890500257   1.587474763160797
    -0.596645326327097   0.794866983837671  -0.110457379005822   2.214695944891245   1.681479560033823
    -0.310918763125585   0.155554394889294   0.937620580494627   2.677700603343778   0.355074529455338
    -0.809367677614890   0.586589743334852   0.028921193727647   2.514447927650235   1.541871099763680
    -0.787244455394507  -0.614436881455177   0.052091133201745  -2.478861421230927   1.518681606685731
    -0.676868273603758   0.729547239539929   0.098031451424216   2.318755802578744   1.472607176020109
    -0.610013758873640   0.769954292442307   0.187226070659293   2.240807702376234   1.382458808505783
    -0.764119674323056  -0.280208869751499  -0.581037100902355  -2.790111341894619   2.190798710893757
    -0.609135766536177   0.737474902338056   0.291692280233569   2.261176182110218   1.274800745458668
    -0.687442424406568   0.650592704599900   0.322725960913292   2.383727704294815   1.242188182072628
    0.388540924306233  -0.065589597811287   0.919094094637873  -0.167233366191138   0.405021070884700
    -0.577782608713844   0.643433214402878   0.502156306015373   2.302487789684919   1.044705868668078
    -0.504217248014560   0.731270821124719   0.459356019854103   2.174446072123354   1.093526261228064
    -0.547560003521984   0.584325360735157   0.598950678557704   2.323724386063816   0.928606225652225
    -0.508565559942047   0.517369243276970   0.688251507337142   2.347613567478788   0.811720182478564
    0.596051760062173   0.436822024287702   0.673727555042823   0.632440450588249   0.831554873170758
    -0.345629112072667   0.456130948167710   0.820051873366244   2.219233283131687   0.609294672040683
    -0.296491170381833   0.359061265700265   0.884967792272315   2.261036392914841   0.484371487205304
    -0.245436510474825   0.256773470476857   0.934787839129935   2.333624180402802   0.363135021219159
    -0.198333244548748   0.146762651826118   0.969084438083561   2.504535190577046   0.249303765291438
    -0.255231051232482   0.226812331961762  -0.939900673771878   2.415081022379144   2.793135618649222
    -0.226600227862065  -0.875541117755810  -0.426708434240252  -1.824050945000030   2.011646419922652
    -0.429230290625778   0.295738942140360  -0.853404848656844   2.538288603183847   2.593279206181812
    0.042637244750974   0.926800202448811   0.373126587233997   1.524823955349519   1.188419615726719
    -0.592171930165494   0.358944073578048  -0.721451008154587   2.596664524473803   2.376691788023011
    -0.615080657129296   0.458664675596740  -0.641328699330640   2.500850243596455   2.267025073518734
    -0.500472230784520   0.547114656815354  -0.670964304945787   2.311700193325530   2.306304847137254
    -0.568262031929270   0.577407440316973  -0.586241341884227   2.348212070862042   2.197207799311458
    -0.641258190932925   0.600647564287935  -0.477504383305957   2.388883090505749   2.068608487729791
    -0.725057638246477   0.509633014055262  -0.463212275534004   2.528932262111148   2.052412687618557
    -0.710034374870380   0.593451895991160  -0.379033024480620   2.445396121550834   1.959547453250601
    -0.768213604791383   0.576487167058700  -0.278406902985655   2.497819206647580   1.852931360596214
    -0.737718567843040   0.648812201123755  -0.186585214667620   2.420228279015810   1.758481493251362
    -0.694171799774780   0.713282079434681  -0.096717048935439   2.342617416240387   1.667664798901242
    -0.780044018025888   0.620335766764236  -0.081944288483431   2.469750246728999   1.652832600994388
    -0.733104855872462   0.680065954175929   0.008219991909962   2.393708722813125   1.562576242313685
    -0.757025309575522   0.642093648710587   0.120948034079145   2.438156202734860   1.449551454596762
    -0.693394107929418   0.688782875143113   0.211619380015996   2.359530687896952   1.357564758038175
    0.465577109338385   0.775327812574854  -0.426737317688536   1.030016125730725   2.011678357128857
    -0.764445536649477   0.599907370160639   0.236080852090209   2.476212838260023   1.332465615909213
    -0.755218940154077   0.556579047506700   0.346213974745514   2.506478446217566   1.217263848629414
    -0.672129465018994   0.603529838309801   0.428944887512402   2.409918564223998   1.127471898877152
    -0.649783979881363  -0.755702675165575  -0.081818373469105  -2.280975801977289   1.652706261743641
    -0.646687195395056   0.549518145944185   0.528985329276854   2.437247604857931   1.013391863646695
    -0.612597599917424   0.487863690234086   0.621862686070318   2.469059719520886   0.899677337572393
    -0.569989874132360   0.419338634886518   0.706588036042281   2.507309856203639   0.786131510982290
    -0.396678425548240   0.536878509043803   0.744585584893687   2.207130076227981   0.730882522323382
    -0.460761694811234   0.444163225097746   0.768386419756003   2.374534805777085   0.694480279082255
    -0.842295169363069  -0.489266024850868  -0.226180469082255  -2.615356301257734   1.798951067112491
    -0.406798111311091   0.362343136252679   0.838583775328561   2.413928452095214   0.576117991665376
    -0.082789630059821  -0.205315508675734   0.975187889102291  -1.954085160357130   0.223228143712550
    -0.191325359343160   0.131441492836581  -0.972685838713149   2.539641461533697   2.907330326653677
    -0.483432510833765   0.192230002432820  -0.854014422380348   2.763126166382653   2.594449777147667
    0.976294787985061   0.042794408398236   0.212172395855470   0.043805449382618   1.356998892658024
    -0.483844316578754   0.376586467163118  -0.789985639150510   2.480210750859688   2.481581901477583
    0.464995164227643   0.883372565965711   0.058586747170869   1.086263230314090   1.512176012157424
    -0.786126852831069   0.410347170693371  -0.462190188950246   2.660511386187802   2.051259760126003
    0.498547426828200   0.501738851872050  -0.706900691557943   0.788588664240037   2.355903077900945
    -0.845638188394645   0.393212447141760  -0.360943244487843   2.706344033356751   1.940075449541605
    0.352538721889175   0.933565430771075   0.064591300017608   1.209725266301287   1.506160029375239
    -0.838673966839436   0.475450935593342  -0.265654635173602   2.625861031671246   1.839679219697781
    -0.815432197159650   0.552035976977102  -0.174145375931496   2.546480052610526   1.745834142839179
    -0.443389868772476   0.843718212809850   0.302564045522396   2.054650024023686   1.263414682527826
    0.405592981058889   0.899273694116657   0.163710588495508   1.147092126552077   1.406345503267985
    -0.906954314001619   0.418549356672299  -0.047437415012389   2.709225729608837   1.618251551316248
    -0.850378360719915   0.522176766724893  -0.064715283450262   2.590904621686356   1.635556867590864
    -0.879648561652580   0.446863479879558   0.162884739409433   2.671563729017796   1.407182588706961
    -0.874491922210625   0.398748611582204   0.276158329133181   2.713778996079772   1.291001636082993
    -0.735841721236882   0.505646493335030   0.450398251623084   2.539531145898773   1.103584981411699
    -0.412769964434358  -0.210823222149582  -0.886100742276704  -2.669380225580352   2.659659859009896
    -0.811896213573240   0.454511938745377   0.366392461608252   2.631244933356226   1.195667434175200
    -0.787968319729122   0.398976567410708   0.468960153702245   2.672889262905161   1.082683253719072
    -0.706776748562075   0.447989661611231   0.547514283632379   2.576655378846181   0.991405501044421
    -0.669164022208897   0.383642104732160   0.636426151928006   2.621025243126889   0.880940280470915
    -0.754249629782682   0.336540089525847   0.563780333210217   2.721910039693709   0.971840547453417
    -0.623712682711038   0.313621590219017   0.715977644604841   2.675683426708220   0.772772834223626
    -0.519100994168406   0.345017309050790   0.781982873411384   2.554991734916182   0.672955576569570
    -0.358209458644344   0.258900923582973   0.897026362770687   2.515751607003774   0.457801473735165
    -0.411819765934831   0.145380317767869   0.899593821450045   2.802230587680328   0.451957754722958
    -0.305963291708745   0.115902865017898  -0.944961898706711   2.779483298159067   2.808274456504906
    -0.369656322401222   0.208948868726229  -0.905369854571508   2.627115330699722   2.703046257856736
    -0.529152696903284   0.086821234298390  -0.844073158342857   2.978965722296687   2.575630688397293
    -0.713777568546343   0.422782979095882  -0.558369174674678   2.606841165408415   2.163215008577104
    -0.688841084300655   0.155392527860449  -0.708061524773545   2.919721038821720   2.357545614701334
    -0.731779102815247   0.246549679077392  -0.635383821347182   2.816619219834314   2.259301796657423
    -0.773098396742108   0.138939949517114  -0.618881700636692   2.963772563148079   2.238114523108246
    -0.767881303463550   0.324016021213556  -0.552604670436339   2.742298815143902   2.156282529882388
    -0.836942380392286   0.307054572545652  -0.453039668662802   2.789962820769527   2.040968371930565
    -0.885652657525080  -0.379511011999001  -0.267564500616686  -2.736752772906546   1.841660809944488
    -0.892469128189153   0.286676557812949  -0.348303612427166   2.830785837740699   1.926557112293915
    -0.893234612736762   0.370439175130664  -0.254768020241005   2.748464323181303   1.828404126625849
    -0.578942887823344   0.810397411746164   0.089895315083073   2.191117380480643   1.480779492765286
    -0.113360481138270  -0.088115847320055   0.989638822382777  -2.480843147442430   0.144077198419794
    -0.316981816768878   0.824195684450197   0.469280301708448   1.937952696075801   1.082320740333157
    -0.879089557746195   0.449462861379335  -0.158696835829515   2.668960059058525   1.730166948654824
    0.220681671483374   0.919097261320121  -0.326435022791915   1.335150191125814   1.903325847439910
    -0.958579284046188   0.272690351404030   0.082254048220402   2.864440562417090   1.488449243617381
    -0.982391145064270   0.158741952801260   0.098532383104005   2.981390076420134   1.472103807332181
    -0.953094077718520   0.230742197582647   0.195881896234985   2.904064855395658   1.373639638988953
    -0.935481827041286   0.174368579638711   0.307358666235843   2.957312962471377   1.258380238051516
    -0.911451692292402   0.288946931909753   0.292857445111539   2.834596714888842   1.273582379768620
    -0.887505124665852   0.230443058550179   0.399037154232397   2.887549977388575   1.160329791053230
    -0.856238782418809   0.345434388936168   0.384096642034594   2.758130895988911   1.176567112475217
    -0.827005405963480   0.284833049984441   0.484698042232214   2.809901991903130   1.064778412382872
    -0.809903716686633   0.091694113071323   0.579351498941049   3.028856630449292   0.952863492167504
    -0.711955229848724   0.267394408198958   0.649322709563573   2.782315725343303   0.864102799519702
    -0.661155952678963   0.193051247171952   0.724985532409114   2.857500139543478   0.759782937886290
    -0.570739070867834   0.239442456769337   0.785445238626593   2.744363814494290   0.667381160788649
    -0.508379788634069   0.162022035063930   0.845753421903879   2.833066905403195   0.562820706862571
    -0.162144381037241   0.032758333726719   0.986223144764662   2.942244130971998   0.166184268764326
    -0.007816030665285   0.049816608399064  -0.998727798347605   1.726423690705204   3.091145225638435
    -0.419453220257478   0.102343709955765  -0.901989335329483   2.902275443397824   2.695151414221974
    -0.591192004542902   0.173518390952077  -0.787644197443203   2.856104268928874   2.477772378580624
    -0.631933748633568   0.068047824124322  -0.772029294113813   3.034324090395699   2.452824104198803
    -0.724048478130819   0.049499536702600  -0.687970636860813   3.073333914968102   2.329485386761300
    0.472787794490889  -0.611300034157148  -0.634652636975464  -0.912478097650259   2.258355276073890
    -0.803682817178285   0.029758534459793  -0.594313350850526   3.104581851763973   2.207207912685019
    -0.847084852117857   0.113295261383920  -0.519241212790758   3.008634466083879   2.116759180783848
    -0.870789622265743   0.003260159541046  -0.491644999073571   3.137848759283527   2.084774150251967
    -0.906107334111967   0.089501971253489  -0.413471759870302   3.043135686241080   1.997060057559528
    -0.927535552557458   0.176336447495410  -0.329519735412964   2.953721770334764   1.906591182026543
    -0.204370066741525   0.694206957322795  -0.690151850120459   1.857101707084945   2.332495193501863
    -0.950637283708680   0.062480292751817  -0.303948955978847   3.075962407102015   1.879631321452052
    -0.964358462297036   0.148423509847704  -0.219050719972763   2.988881868616581   1.791637782265249
    -0.964319467033664   0.231906706161704  -0.127699824342631   2.905586844825657   1.698845795098942
    -0.950067071563043   0.310435810735997  -0.031657020458479   2.825777427709098   1.602458637242619
    -0.996693134618529   0.081051960618266   0.005777117293682   3.060450328469356   1.565019177365436
    -0.992416377012296   0.039244263950007   0.116488722135708   3.102069094639420   1.454042531494742
    -0.921677821229058   0.382276727006583   0.066139986722037   2.748425735804459   1.504608023361476
    -0.970549648950547   0.117552261591609   0.210273262010811   3.021060496224970   1.358941864295886
    -0.947324627276778   0.054068799213910   0.315678025061104   3.084579248781070   1.249625194323449
    -0.906363278526073   0.109498878804227   0.408063233923476   3.021364072639131   1.150464701607175
    -0.975961328653975  -0.000117879625971   0.217943733739609  -3.141471870499381   1.351089267392237
    0.139393958975688  -0.542924767031562   0.828131645056997  -1.319478294185687   0.595030068055317
    -0.998338239884531  -0.039455239840782   0.042000509917883  -3.102092296149320   1.528783458614518
    -0.251367436823795  -0.915632913402052  -0.313736799880897  -1.838724243303587   1.889922317378962
    -0.740776711391030   0.143711961780283   0.656198701539369   2.949971259831353   0.855026250645116
    -0.752240890040623   0.018694731933126   0.658622919695967   3.116745713160798   0.851809107361325
    -0.679863306134805   0.068731153657064   0.730110891227094   3.040839417948220   0.752312109279877
    -0.598047289344066   0.116270260731213   0.792982134841454   2.949571563041258   0.655107999075902
    -0.404313578233443   0.032403642593775   0.914046243032990   3.061618761536890   0.417646085410436
    -0.278525412185390   0.052651434413673   0.958984578197763   2.954760920858192   0.287398353907075
    -0.122847421393513   0.033382807979323  -0.991863951955297   2.876257864857775   3.013943880856693
    -0.350364325746813   0.006304273575875  -0.936592278090452   3.123601114850316   2.783572044124563
    -0.460372719858537  -0.004936641523034  -0.887711996303152  -3.130869923704477   2.663147826446319
    -0.593503845424195  -0.126761790636796  -0.794785904442981  -2.931172152481821   2.489451049814016
    -0.644489434719432   0.258758506316061  -0.719498022196106   2.759799189424276   2.373875578891272
    -0.663436728827398  -0.039161966515820  -0.747206830282898  -3.082632131582526   2.414645579734788
    -0.749732013802439  -0.058080873077264  -0.659187772688725  -3.064278227460949   2.290534455053144
    -0.881859941919853  -0.110250093483926  -0.458440791950170  -3.017218007486280   2.047036297180274
    -0.922717704062091  -0.023952549905749  -0.384731482937640  -3.115639784363569   1.965713222973209
    -0.015643244019726   0.991248014694114  -0.131082654388428   1.586576379137628   1.702257305320363
    -0.523028589613906   0.778771137157678   0.346347528324694   2.162211758534870   1.217121487286108
    -0.960093464115781  -0.053610534100770  -0.274493079686156  -3.085811710106692   1.848858820113472
    -0.980978873390798   0.031836412526078  -0.191486012017561   3.109150321467589   1.763472278960548
    -0.982800627033794  -0.086084368243777  -0.163378116790625  -3.054224755442784   1.734910141117387
    -0.997031118162989  -0.002009030268325  -0.076973457841264  -3.139577643714313   1.647845998189168
    0.072353438766014  -0.046216437029391  -0.996307693861209  -0.568432566819461   3.055632423644588
    -0.987403219046904   0.116739464200269  -0.106802530466494   3.023910181964106   1.677802952140230
    -0.986716677439866  -0.146664943330225   0.069854082629244  -2.994033662300419   1.500885309144154
    -0.435635728150566  -0.316608545291302  -0.842603430688038  -2.513122902625850   2.572895709034071
    -0.980307279904769  -0.098035120164425   0.171425646214503  -3.041919565055732   1.398519772234102
    -0.951410782532537  -0.074062271187243   0.298885099774835  -3.063904641838949   1.267272191341125
    -0.315739342115758  -0.602087954171512  -0.733347641491995  -2.053778610925971   2.394029349615483
    -0.908395482505696  -0.137979171677764   0.394688985843743  -2.990851601165776   1.165066995298695
    -0.867281829069875  -0.085559450491181   0.490409838193382  -3.043258416850767   1.058236364122151
    -0.815010225549669  -0.032864253431168   0.578513848663876  -3.101290756389598   0.953890810900028
    -0.747241686933897  -0.102531396399575   0.656595137097844  -3.005231015276400   0.854500763753881
    -0.677448324234091  -0.051181328698104   0.733787598413273  -3.066185746585358   0.746915943330599
    -0.860247464578499  -0.345726348681905  -0.374763380699531  -2.759456532388461   1.954937868544006
    -0.597583478947777   0.000615592148453   0.801806464637806   3.140562518131923   0.640484262249931
    -0.505800898946855  -0.057593754927487   0.860725513748666  -3.028214520570262   0.534103190125823
    -0.915671443702225  -0.017950561094085   0.401526567669799  -3.121991450290059   1.157613252233398
    -0.330940133218460  -0.043883333386388   0.942630829793004  -3.009759719747855   0.340371141225101
    0.293747224618185   0.955257835800141  -0.034569280749601   1.272467877445715   1.605372496500844
    -0.385277504348566  -0.102824449922775  -0.917054184408435  -2.880787250861683   2.731425508486109
    -0.492026071735958  -0.112155174965582  -0.863325872113492  -2.917476327722266   2.612619800152385
    -0.613498182553536  -0.232980604553211  -0.754546233113338  -2.778657464567010   2.425758713580515
    -0.685955136574846  -0.147379365259615  -0.712562189077036  -2.929957217141589   2.363939673235526
    -0.698430898100357  -0.257845732952857  -0.667615052688847  -2.787935428852074   2.301797101591999
    -0.765147691140893  -0.167110704325414  -0.621790176218332  -2.926565940616273   2.241822729180756
    -0.829401435560007  -0.194878624230830  -0.523560484097968  -2.910815785964803   2.121820965592006
    -0.823943021147971  -0.305882330557125  -0.477026097561230  -2.786120651393077   2.068064213818677
    -0.880771077265700  -0.224972389594026  -0.416689012781201  -2.891513129714536   2.000596339817438
    -0.926058027309778  -0.139885453338179  -0.350497631945631  -2.991671336636282   1.928898715761754
    -0.913739669644551  -0.261573591046167  -0.310900422298008  -2.862781991301423   1.886936584308561
    -0.954368225765734  -0.173514898449849  -0.243051166762686  -2.961745845180067   1.816306433498852
    -0.934210123729407  -0.294305083327064  -0.201583636859090  -2.836403691963737   1.773770807605725
    -0.942257155952560  -0.322920798296545  -0.088733365109835  -2.811427114151752   1.659646548403567
    -0.969704554184277  -0.204584767195008  -0.133484645656483  -2.933665492388495   1.704680593298221
    -0.991838076120760  -0.117568107448137  -0.049345423983728  -3.023607602607944   1.620161798531778
    -0.972471481625234  -0.231989265922665  -0.021913418772191  -2.907413125662250   1.592711499742533
    -0.962300056372336  -0.256851716938626   0.089475119510835  -2.880759014276331   1.481201388512611
    -0.959287015449635  -0.210205878403530   0.188631679932964  -2.925875006892956   1.381027700532678
    -0.891711058402302  -0.257575748041374   0.372164106737838  -2.860391160873453   1.189456800742071
    -0.938836018642069  -0.184337542828055   0.290872137552609  -2.947712143834090   1.275658063926830
    -0.854708498852562  -0.207244840409096   0.475944280470947  -2.903709452695509   1.074758920611249
    0.216500272252258   0.651656508809691   0.726960402388371   1.250038506142958   0.756911327883689
    -0.732434142711027  -0.220935173706441   0.643993692213261  -2.848627272896011   0.871089155055471
    -0.644523627745025  -0.285259732456124   0.709377317312481  -2.724917919979497   0.782181962901239
    -0.666452505331318  -0.169714321312756   0.725973902615625  -2.892239243754831   0.758346861108835
    -0.573780753767809  -0.234000912599300   0.784868918678998  -2.754360669044870   0.668311764971295
    -0.495882168415619  -0.176748601777631   0.850212212813410  -2.799198710321941   0.554408054977141
    -0.415636503088480  -0.104712244056884   0.903483061957971  -2.894796220146062   0.442969002743076
    -0.223114070369717  -0.065703677556839   0.972575518074851  -2.855203426868907   0.234737117497208
    -0.044648875844192  -0.064103901729782  -0.996943914003626  -2.179174772508158   3.063392348264806
    -0.304867145760329  -0.202811474279058  -0.930550121883128  -2.554574631358399   2.766708704435759
    -0.516198457223451  -0.218164214715368  -0.828217077931720  -2.741725769338700   2.546715016160436
    -0.535727406540154  -0.321765306419334  -0.780681262402674  -2.600721906716629   2.466551545923850
    -0.628537304053682  -0.338522866207613  -0.700245047442330  -2.647552985911619   2.346537016189900
    -0.685491681054090  -0.387005717688225  -0.616707166880895  -2.627634688469995   2.235349133184119
    -0.754608421696578  -0.382509168690467  -0.533153698076200  -2.672442233640847   2.133120226928351
    -0.792566804930080  -0.424887890143479  -0.437387860522383  -2.649490853951388   2.023488217602729
    -0.042619204733676   0.786836793448592  -0.615687797396868   1.624908689767237   2.234054846176452
    -0.744253664129848  -0.526671429130031  -0.410735546510679  -2.525750601958632   1.994056979968133
    -0.822743510212770  -0.459538511100438  -0.334540689929287  -2.632213368501457   1.911914121685473
    -0.144460242886902   0.623845048866448  -0.768081111100824   1.798349913605308   2.446635455956345
    -0.771383245057421  -0.561495896862802  -0.299483300120258  -2.512381107192633   1.874947378259276
    -0.899264771329903  -0.408801718480283  -0.155576431417310  -2.714923685102103   1.727007290002566
    -0.850877032064875  -0.513108048644067  -0.112820240742289  -2.598945025365146   1.683857285869726
    -0.901212175446936  -0.431718503755619  -0.037891269986768  -2.694851559180682   1.608696669699451
    -0.851696499545617  -0.524007162971105   0.005437445789704  -2.590034556995753   1.565358854211082
    -0.936955060973166  -0.348344259131620   0.027775724055608  -2.785645068228834   1.543017030046003
    -0.895122743761097  -0.435674518794318   0.094567368975477  -2.688625034078537   1.476087435077245
    -0.916711651321088  -0.302048521001187   0.261546246949100  -2.823303898326779   1.306172459322874
    -0.523563943828606   0.753156542026387  -0.398291375660449   2.178266257187036   1.980449667585476
    -0.221722004714494  -0.300718408775927  -0.927581689798080  -2.206124767752788   2.758683625896072
    -0.795497780386855  -0.429447916368265   0.427501776050747  -2.646577051286121   1.129068839161944
    -0.832302742979124  -0.322443768581596   0.450890408117454  -2.771985193570114   1.103033671307790
    -0.706107146625813  -0.334630132377274   0.624047571895861  -2.699035313425334   0.896884283253567
    -0.786908161575569   0.495986345631317  -0.367114546419846   2.579192876954000   1.946701389196102
    -0.504218556532034   0.051733926717181   0.862025085525458   3.039348247468486   0.531544815566633
    -0.469571223883240  -0.290496673303976   0.833735298820986  -2.587578157198698   0.584957961433165
    -0.823275492878189  -0.081040440134772  -0.561827295428887  -3.043472160594528   2.167389341441588
    -0.435817898171915  -0.399639614952493   0.806443387841724  -2.399471011432270   0.632683911017239
    -0.306130251653639  -0.157166565193798   0.938926482642631  -2.667284426684624   0.351299039315071
    -0.276940852179436  -0.268067921132078   0.922736882352740  -2.372473397402749   0.395674289710389
    -0.274476525835983  -0.092820571543194  -0.957103431329866  -2.815492671618455   2.847630710705010
    -0.193108854169428  -0.191153277529625  -0.962376950540199  -2.361283619710463   2.866415034547996
    -0.871904698082033   0.487271752312803   0.048460673339211   2.631973390422013   1.522316665577249
    -0.458083452068436  -0.418449451503828  -0.784257360486473  -2.401380497650450   2.472294582049329
    -0.491708961938940  -0.510046454361582  -0.705744225015061  -2.337891175508288   2.354269396578212
    0.456727104165480   0.877244006002246  -0.147794811322138   1.090774478204624   1.719134551275897
    -0.611743584741182  -0.470641068730479  -0.635819762945634  -2.485826086468254   2.259866475381468
    -0.705090708249702  -0.482118237179054  -0.520008748502653  -2.541844139845243   2.117657519677094
    -0.661355771328236  -0.590354129961212  -0.462699194908068  -2.412857489189342   2.051833839011572
    -0.690675664521280  -0.627876480390098  -0.358801131284490  -2.403785907052709   1.937779511511268
    -0.711567547491419  -0.657049790628807  -0.248912028619930  -2.396007626554610   1.822353092643924
    -0.787167199743279  -0.588420838529512  -0.184739590863873  -2.499682580163975   1.756603212594387
    -0.723221557088109  -0.677655564380789  -0.133167246073519  -2.388709785607312   1.704360334521503
    0.519652318784637  -0.853589886100248  -0.036684791523391  -1.023942637306516   1.607489351544629
    -0.722110484894321  -0.691653512650016  -0.013109769051452  -2.377734338305523   1.583906471396086
    -0.640779507020692   0.767703389541004  -0.005755785747659   2.266323335407988   1.576552144323668
    -0.770644022049108  -0.613722090061572   0.171618726980566  -2.469065626224671   1.398323787017228
    -0.841191442728033  -0.527456370160871   0.119108078054549  -2.581531396247705   1.451404809727077
    -0.695198430735408  -0.679738276958675   0.233784124226344  -2.367438282127732   1.334828476885881
    -0.879195228774124  -0.421104516611280   0.222905217055134  -2.694913491927597   1.346002668141246
    -0.751562367115701  -0.592066893238218   0.290879360331200  -2.474346577083548   1.275650514727162
    -0.747684833663585  -0.538491881144474   0.388579314234485  -2.517429737087241   1.171707088065931
    -0.668312314274774  -0.595002857828979   0.446464164031161  -2.414158912806048   1.107986432305034
    -0.716949442307529  -0.491748595324112   0.494132387292804  -2.540399256934561   1.053959738037521
    -0.612328455817901  -0.396049055399697   0.684250690838287  -2.567476029357966   0.817220634902296
    -0.573404399873234  -0.501337767959098   0.647972095559528  -2.423149069926531   0.865877381037418
    -0.544423148036633  -0.346286313518995   0.763995566054743  -2.575078895596689   0.701313161754431
    -0.508019517599980  -0.453790166888739   0.732113825967366  -2.412517477901271   0.749376356753413
    0.058540127106452   0.756159595305351   0.651763546039519   1.493532753246805   0.860888923996643
    -0.391134834900928  -0.223875621605989   0.892688773861056  -2.621733714393132   0.467519843392930
    -0.926699211996422  -0.341714567363501   0.156396051538349  -2.788318140844888   1.413755586326645
    -0.250008468811975  -0.407313989920496  -0.878402572364937  -2.121299609696733   2.643305742524476
    -0.356270483894656  -0.415591828523725  -0.836872017915389  -2.279490447548770   2.562340069254451
    0.264925712371587   0.544653442118434  -0.795717911645169   1.118078470920512   2.490988308309456
    -0.388159736186914  -0.513811356624232  -0.765068564905964  -2.217777493950801   2.441944126193301
    -0.350521455601612  -0.682509286818570  -0.641339054299155  -2.045246891763712   2.267038569520823
    -0.386261902722077  -0.745820699677497  -0.542727580318237  -2.048662615843231   2.144477513757738
    -0.459965247129298  -0.665559296183176  -0.587760831203849  -2.175521299578630   2.199084671718233
    -0.560781193543194  -0.639903540707845  -0.525402618529670  -2.290391942144970   2.123984564378869
    -0.523405809574348  -0.763794601092701  -0.377709102148021  -2.171568077818359   1.958117196919707
    -0.622737634345269  -0.716511356860789  -0.314339488864459  -2.286288951572901   1.890557120851006
    0.307947418948058   0.019771542362730  -0.951197914881981   0.064116276985294   2.827891286335697
    -0.132093088057585  -0.671747864949882  -0.728907553824686  -1.764959882342396   2.387521203443512
    -0.641240877117543  -0.740635892470443  -0.200670407130785  -2.284390593827990   1.772838526840622
    -0.563481492386513  -0.812483871473856  -0.149527811234975  -2.177166426333770   1.720887024609735
    -0.713508508064006  -0.692426650013431   0.107009080322434  -2.371188267165692   1.463581961064895
    -0.642931235711774  -0.764692469603260   0.043414894658380  -2.269908007487047   1.527367782105969
    -0.465069989194342  -0.762909347171040   0.449087110870299  -2.118245195873613   1.105052964352228
    -0.610917480139399  -0.735531568602788   0.292870524376982  -2.263907217525428   1.273568700733749
    -0.549486884543963  -0.805823990905514   0.220707633296360  -2.169267263937234   1.348256391055279
    -0.669455702025222  -0.653864900376019   0.352547521735973  -2.367975520668152   1.210504303701770
    -0.819398008902589  -0.473117952648805   0.323645339666550  -2.617958830133012   1.241216667507725
    -0.582147231779539  -0.642582566515573   0.498188965897158  -2.306888632889271   1.049287493879367
    -0.670871904945320  -0.443688582250077   0.594198055479783  -2.557272133181003   0.934528093866306
    0.338159024494732  -0.702853203265981   0.625816146173556  -1.122363829090449   0.894618851878795
    -0.529774950287569  -0.600500786087902   0.598946832327892  -2.293702002075317   0.928611028712346
    -0.466960858617112  -0.556156619409885   0.687486269830997  -2.269233986973250   0.812774468476205
    -0.397452203819408  -0.504240119231778   0.766663973221788  -2.238313023538294   0.697167328732443
    -0.322826931521648  -0.443804485126871   0.835954754317236  -2.199681081136717   0.580926121869375
    -0.243112325589077  -0.376303211342187   0.894031481705330  -2.144392799314208   0.464531660451068
    -0.194544144822348  -0.182227388243598   0.963818320374363  -2.388873105386290   0.269822146004487
    -0.160530131226000  -0.079826252599880  -0.983797665358291  -2.680134254231216   2.961335740943497
    -0.139542093568676  -0.394274161484411  -0.908336881178033  -1.910960386294862   2.710086580611513
    -0.285878849621043  -0.954922650917033  -0.079975084275886  -1.861678552754919   1.650856911014436
    -0.282488420554898  -0.509066817866411  -0.813050593259586  -2.077404627183797   2.520169264327669
    -0.242150309859538  -0.682277996297934  -0.689826038362286  -1.911842536028205   2.332045065278578
    -0.421635125120269  -0.599589842952274  -0.680232196748495  -2.183670674023394   2.318875692257662
    0.194500500989228  -0.839935837857419   0.506633342169358  -1.343241219576897   1.039520940380823
    -0.416703563793610  -0.797568588733457  -0.436167956396876  -2.052255856958019   2.022132119540346
    -0.623375245654979  -0.546201565297331  -0.559524041638373  -2.422083253135815   2.164607751104924
    -0.444751808979075  -0.834206581919852  -0.326029457398311  -2.060605839973966   1.902896809211134
    -0.547495256292497  -0.793443500001120  -0.265887864791131  -2.174795733532841   1.839921150070653
    -0.465824100972103  -0.859132641563632  -0.211893867663571  -2.067633480393736   1.784308752349914
    -0.377711785778452  -0.913212905548596  -0.152892105824806  -1.962978064897330   1.724290454705085
    -0.479715497486407  -0.872376560931023  -0.093979664872410  -2.073558909834989   1.664914885253811
    -0.679532779997508   0.352634029387410  -0.643338512936119   2.662911250075684   2.269647397886736
    -0.486695950221698  -0.873112605554599   0.028309540078531  -2.079323245309157   1.542483003999263
    -0.596062674650893  -0.682719658794317  -0.422614665367632  -2.288532700112850   2.007124671638519
    -0.048250932078208   0.010643333846985   0.998778537513800   2.924486360946589   0.049430986819886
    -0.459302499447721  -0.873475369655670   0.161499203112484  -2.054896409537853   1.408586716874703
    0.175445149299345  -0.688530060850746  -0.703665655615078  -1.321294656254893   2.351339768806485
    0.381340428621184  -0.830284885648282  -0.406456007660153  -1.140244856473165   1.989368174566496
    -0.517549110607933  -0.782860333655206   0.345358677461849  -2.154935304514799   1.218175375025786
    -0.385790253305371  -0.740182052891729   0.550723532302338  -2.051267566289009   0.987565506624785
    -0.488723442859953  -0.692596000933550   0.530528204613071  -2.185297505784547   1.011572755588039
    -0.316827594515405  -0.695157486042823   0.645272302948302  -1.998433053717760   0.869416657638545
    -0.357433596078431  -0.603720951979455   0.712574372634507  -2.105350475225110   0.777635615033015
    -0.284189257916323  -0.546500707422574   0.787764839575616  -2.050328423517331   0.663624444314363
    -0.207074711428829  -0.480887980442698   0.851978177040005  -1.977408272358370   0.551044372539341
    -0.162600917653824  -0.295532562383320   0.941393353571865  -2.073790223197594   0.344058778865075
    -0.001745073625791  -0.110297138299611   0.993897125461663  -1.586616574395608   0.110535896191672
    -0.078195809790042  -0.175898160674160  -0.981297738916547  -1.989110138094865   2.947887468624645
    -0.060249549561737  -0.481095036758025  -0.874595653650533  -1.695381906516920   2.635397558433846
    -0.173262880745945  -0.497753409053004  -0.849836171230508  -1.905768389543224   2.586470699824953
    -0.207666911272584  -0.594722276338544  -0.776646552808422  -1.906743125184054   2.460121065489076
    -0.756744347344323  -0.383605113393108   0.529325145577909  -2.672428304603453   1.012991384825171
    -0.166055631398342  -0.753385211983671  -0.636267435631604  -1.787740197890499   2.260446622979677
    -0.275847182472218  -0.756880965035531  -0.592485895772230  -1.920287988846835   2.204937521394643
    -0.197223135916919  -0.821836255846114  -0.534498085343566  -1.806321098017023   2.134710099746689
    -0.308246356204217  -0.815860294748103  -0.489240394223904  -1.932034906665465   2.082014908008351
    -0.336166408668126  -0.861795910471094  -0.379868338215325  -1.942718876150041   1.960450288286106
    -0.141873188773803  -0.923704438480473  -0.355868105677689  -1.723196920509261   1.934639146027355
    -0.360119882732546  -0.893709531216252  -0.267576052504563  -1.953842858631513   1.841672798972124
    -0.271070222371547  -0.941941223457137  -0.198160707748818  -1.851003162091212   1.770277386214927
    -0.179462513363213  -0.976202966880615  -0.121741421672507  -1.752603631904444   1.692840493051369
    -0.389301102478860  -0.920387286922301  -0.036495118585236  -1.970950798560129   1.607299551509134
    0.022048870991530  -0.996560126341383  -0.079885930391059  -1.548674957780988   1.650767470960617
    -0.297003101931939  -0.953948779720702   0.042081862033289  -1.872624707586078   1.528702034510289
    -0.227194791590957  -0.964808113254983   0.132392716080991  -1.802064862969378   1.438013768395302
    -0.332738756419274  -0.927605684998406   0.169801687668860  -1.915206882332194   1.400167895765276
    -0.357437413550080  -0.893360069871266   0.272298147173712  -1.951392815307668   1.295015701939345
    0.175990193328259  -0.031922518532445   0.983874181317321  -0.179437103694720   0.179829621393281
    -0.297863378116080  -0.879774704631591   0.370504894795838  -1.897250471912583   1.191243783802150
    -0.192620003490132  -0.899584970351075   0.391975018813594  -1.781732139947650   1.168018897484231
    -0.347794619696766  -0.811917935366843   0.468858367462526  -1.975510990741784   1.082798494773212
    -0.160468440570081  -0.804918876155307   0.571275312253106  -1.767576203418666   0.962737488730981
    -0.275591363948809  -0.775215348302650   0.568410559255224  -1.912365138865153   0.966223638881581
    -0.127131815280679  -0.672290379624575   0.729289480938594  -1.757691935939052   0.753513410627248
    -0.243958476030739  -0.642270193678505   0.726617685089577  -1.933801566409576   0.757410289741300
    -0.491717970129071  -0.723426313532159  -0.484631619626099  -2.167772640620983   2.076738303616913
    -0.293159210377627   0.765032252104653   0.573396311995855   1.936735277382484   0.960150962908621
    -0.127272894395555  -0.404574459054632   0.905605387259169  -1.875579203878084   0.437991381838097
    -0.048326204909353  -0.319305019446174   0.946418978294253  -1.721004438655277   0.328835722363293
    -0.026035972512980  -0.375600308730019  -0.926415962846721  -1.640003914087177   2.755575528858408
    0.018613109294328  -0.559144282146108  -0.828861402108767  -1.537520052284403   2.547865737308296
    -0.096351349257883  -0.579955456653094  -0.808930210707010  -1.735428416541169   2.513126492675631
    0.025459530684669  -0.795426342714263  -0.605515272816077  -1.538799848364362   2.221209497802478
    0.119047461656912  -0.892619983266886   0.434807161102102  -1.438210161442450   1.120972201959953
    -0.086466288568767  -0.814245170905284  -0.574045627627777  -1.676591797409291   2.182234456717930
    -0.005305778411305  -0.863938350929148  -0.503569632234982  -1.576937635028030   2.098521886601159
    0.730891096215978   0.278279485354751  -0.623184349533026   0.363793489271099   2.243604125208687
    0.666628342390823   0.220381594350612  -0.712066433693336   0.319280784263473   2.363233328676783
    -0.031762190630124  -0.918820082014439  -0.393396517693478  -1.605351022113165   1.975119413694894
    0.040136028299403  -0.979323042672629  -0.198281308556236  -1.529835808377152   1.770400428507505
    0.132619722021371  -0.984541775349734  -0.114409361164933  -1.436900300700744   1.685456763284801
    -0.162342768808178  -0.957045256440325  -0.240227397564824  -1.738826023774332   1.813396428454707
    -0.070186079735733  -0.984239314038917  -0.162317241571963  -1.641985795015185   1.733834912859869
    -0.196249637139157  -0.980553672650932   0.000758270082644  -1.768328091898516   1.570038056639589
    -0.088386694282201  -0.995174202781166  -0.042615705939665  -1.659379195097914   1.613424942340158
    -0.133998064898977  -0.972170001318565   0.192171816715124  -1.707767261779109   1.377421592398993
    -0.980097089563054   0.197986723950808  -0.014524192550977   2.942267786752446   1.585321030046010
    -0.250102526960179  -0.933531370076737   0.256842183238625  -1.832559284728874   1.311042977568794
    -0.150549121815036  -0.942170556060526   0.299415439169928  -1.729246501865180   1.266716400037773
    -0.117709601578009  -0.868339254961163   0.481800153580131  -1.705532199813617   1.068088461250781
    -0.555470674405470  -0.420564532511930  -0.717340786424802  -2.493542771870015   2.370774360160382
    -0.237999540333437  -0.843533782372375   0.481463370151486  -1.845793748051586   1.068472753289808
    -0.201940938563961  -0.728722353009367   0.654357386721048  -1.841129101271456   0.857463828235732
    -0.530615943197789  -0.574243755037500  -0.623450744345336  -2.316727797987847   2.243944797124486
    -0.084651402041102  -0.754949017982527   0.650296947847439  -1.682458549354573   0.862821070058297
    -0.051731580539668  -0.605343860177744   0.794281218788393  -1.656047370492354   0.652972744951177
    -0.091181250691670  -0.508136083386391   0.856436629461305  -1.748349287926245   0.542469071962556
    -0.011284642007099  -0.427497908494355   0.903945902744032  -1.597187150785778   0.441887940910499
    -0.065009568038773   0.337721294434749   0.938998446936271   1.760964971908677   0.351089852700920
    0.895628093734424  -0.443132180821922   0.038525161046373  -0.459456897391350   1.532261629614275
    -0.237286466773849   0.018489978128310  -0.971263740389189   3.063827196348179   2.901280552887855
    0.182407269321219  -0.961578300510730   0.205170075999578  -1.383328093679434   1.364158849914543
    0.096957383557795  -0.628149396742908  -0.772028238567183  -1.417650914260557   2.452822443457393
    0.994722029524314  -0.042215277094501  -0.093519807307677  -0.042413818803515   1.664452993418542
    -0.330810994553501  -0.310834240941634  -0.891036565209735  -2.387318079591147   2.670419930883942
    -0.652820099672530   0.516490867901617  -0.554132746584951   2.472259859994165   2.158117095664381
    0.105831199329623  -0.840205178393218  -0.531841156220233  -1.445497437679459   2.131569549815060
    0.077844374122374  -0.901522395945064  -0.425673141068317  -1.484662284921282   2.010501985749315
    0.167724991919897  -0.923874773016364  -0.343982166499126  -1.391207126458598   1.921950914131219
    0.936553678313331  -0.124500803130766   0.327668670546276  -0.132160178214032   1.236961367501118
    0.150264479143138  -0.961249816811023  -0.231126320415085  -1.415729316511356   1.804031516432532
    0.947394819748969  -0.006504092655912   0.320001175453376  -0.006865131981551   1.245065598807996
    0.242279551173453  -0.958810753187934  -0.148265163320240  -1.323289797938525   1.719610142932466
    0.112576579086177  -0.993634533336780   0.004114365303844  -1.457979632682870   1.566681949882967
    -0.806532160099170  -0.154775302694630   0.570570311531850  -2.951995418023945   0.963596194652853
    0.092570273257410  -0.988185696191792   0.122146528198688  -1.477391908314583   1.448344007960572
    0.000045802083694  -0.999211795518706   0.039696166042241  -1.570750488581349   1.531089727911870
    0.223671593233393  -0.974165182811984  -0.031196393659020  -1.345104689151755   1.601997782804194
    0.072427759653119  -0.968243434737062   0.239288258634238  -1.496132127373997   1.329163579188736
    -0.038131682461361  -0.961669322284192   0.271547582147456  -1.610427119089756   1.295795656414065
    0.042778470490493  -0.933679843670892   0.355544585087710  -1.525011280823964   1.207299667805818
    0.078627048003855  -0.837810801534475   0.540269236727726  -1.477222187491860   1.000039298708810
    -0.070884000978605  -0.920440728428397   0.384401253712142  -1.647655558602224   1.176237171032547
    -0.040833284840212  -0.825897976357746   0.562339021851888  -1.620197177126804   0.973584598124184
    0.034251860019961  -0.771675058041675   0.635094020505271  -1.526439064702553   0.882666093190873
    -0.019714791503239  -0.987469449911329   0.156573984070795  -1.590758638230656   1.413575434381881
    -0.009794312520270  -0.694061941651102   0.719848659506815  -1.584906972940871   0.767212061223590
    0.064550431917511  -0.623008453585284   0.779547117562843  -1.467553873768493   0.676853892774915
    0.024516089054694  -0.528736071614844   0.848432158720286  -1.524462165575850   0.557780183030345
    -0.787775167314278  -0.272474098962803   0.552420266787316  -2.808594704135396   0.985531362609691
    0.028979471688636  -0.226555898352009   0.973566954627447  -1.443574065901315   0.230435779492975
    0.039394972200691  -0.157640493338853  -0.986710449435491  -1.325908124012279   2.978380578750695
    0.054749250800012  -0.456849036363406  -0.887857802528459  -1.451524134901387   2.663464618430133
    0.133664805283580  -0.528947992397447  -0.838061895785277  -1.323278732256321   2.564517390323607
    0.210250945537370  -0.593149051979486  -0.777154258842142  -1.230148780718568   2.460927430094888
    0.283482885522147  -0.650622151194793  -0.704505691950532  -1.159888900611873   2.352522727181708
    0.265302567109698  -0.735322085995224  -0.623631283478175  -1.224534686801863   2.244175727407979
    0.140301098297744  -0.766968908650827  -0.626158362540505  -1.389867459266019   2.247412634446381
    0.187288262338278  -0.872373545365451  -0.451539039438255  -1.359318257754580   2.039285808976252
    0.277311737109345  -0.886093046183500  -0.371399130271792  -1.267492379611419   1.951311809071846
    0.259807855115920  -0.929320669313064  -0.262417552781028  -1.298188726950614   1.836323033803823
    0.367889797975812  -0.881935825889791  -0.294696955460025  -1.175602842514109   1.869934702000030
    0.332741568565785  -0.940513350718148  -0.068685410891277  -1.230751617477946   1.639535858691510
    0.350392528965758  -0.918474712518734  -0.183382872996935  -1.206344537198394   1.755222911865503
    0.407758973277530   0.187233754586336  -0.893686824819052   0.430459691849220   2.676292266402968
    0.312098061648994  -0.948758223098083   0.049524075142386  -1.252992099988413   1.521251985204846
    0.203375392150468  -0.975224602824809   0.087036911209189  -1.365200858008250   1.483649149028672
    0.839020724764904   0.251507539134068   0.482481275464345   0.291239534596811   1.067311008056147
    0.290810255401778  -0.941854577078929   0.168342956456824  -1.271319199874659   1.401647934278395
    0.944453023570496   0.275263446804859   0.179550887276161   0.283596967732132   1.390266426553837
    0.154889790271132  -0.934000353986520   0.321951070230897  -1.406457112626159   1.243006765089052
    0.344479432114480  -0.867514850501568   0.358820156916146  -1.192803116022356   1.203792759155604
    -0.425308429941063  -0.653107291251372   0.626548964994247  -2.148020005075603   0.893678950621661
    0.884581282320973  -0.260351688819004  -0.386953424970598  -0.286239469007772   1.968121664228724
    0.267357732710913  -0.776395771417827   0.570727122961751  -1.239157020511169   0.963405239315280
    -0.628920033209666  -0.546618771758436   0.552872055894355  -2.426091958096932   0.984989269438096
    0.151893193864856  -0.779396628863730   0.607839907026009  -1.378323243428872   0.917458893415782
    0.950847598430320  -0.062718534054040  -0.303241207696787  -0.065865247060821   1.878888513463945
    0.107881104933088  -0.706559662707401   0.699381948747049  -1.419281443011808   0.796263908382578
    0.254114583067576  -0.546244216430366   0.798149757056135  -1.135371554036726   0.646578538181117
    0.222334441917799  -0.457974728649774   0.860712811485495  -1.118837575756447   0.534128141530895
    0.105077522870167  -0.445451161651682   0.889118651683023  -1.339141008893373   0.475380474412472
    0.065975600680430  -0.339727657375452   0.938206980857126  -1.378982470541773   0.353383969477562
    0.107103647720106   0.066355627729487  -0.992031118117609   0.554667661377621   3.015263810662936
    0.090502073715232  -0.350037537132742  -0.932353525896313  -1.317787285938640   2.771665058503513
    0.171642376792242  -0.423200735721359  -0.889629154071521  -1.185487616791602   2.667328811784091
    0.321657208571626  -0.547307251425760  -0.772652193882014  -1.039463605944123   2.453804727022209
    0.387199437852026  -0.599066490814330  -0.700853718624992  -0.996999755990380   2.347389970094780
    -0.052375551505500  -0.958748269729522  -0.279389614866305  -1.625371176246364   1.853954677057089
    0.371597895564385  -0.687899611973676  -0.623465418334158  -1.075514390496774   2.243963565253435
    0.359927503629032  -0.763139833847093  -0.536721330046798  -1.130092767372431   2.137342825079334
    0.456089880432462  -0.764643227513884  -0.455305123608613  -1.032973433865097   2.043511217244309
    0.544559691481608  -0.752846108177551  -0.369699174756140  -0.944580570705842   1.949481563235782
    0.464837268903539  -0.821470456335011  -0.330321968399088  -1.055858130987284   1.907440997400618
    0.546488937762235  -0.800376816448694  -0.246468644242950  -0.971714617276134   1.819831120969980
    0.554568025273584  -0.820318409858997  -0.139756974031939  -0.976332562089606   1.711012302546363
    0.849560660566845   0.503004641049242   0.158849032418045   0.534572798594544   1.411271552947533
    0.440220422572811  -0.891450933327135  -0.107336913594509  -1.112101211059626   1.678340424929291
    0.417535386038387  -0.908618269010545   0.008777392896328  -1.140047359591213   1.562018821189094
    0.396223208579389  -0.909133463322469   0.128388141392990  -1.159792543102084   1.442052828511515
    0.301024792867596   0.157254694969006  -0.940561021406493   0.481404890314161   2.795074757910927
    0.265916000327583  -0.920843861702854   0.285193378488976  -1.289669931756800   1.281588147198284
    0.374716767466803  -0.894227989686532   0.244833916442174  -1.173985167369431   1.323447934942486
    0.555071121735036  -0.788873178632872   0.263771032999224  -0.957641421272718   1.303866715280909
    0.452538631504968  -0.834552850656408   0.314213822829093  -1.073920502451062   1.251167905933384
    0.309031916447688  -0.828113085287979   0.467683645845664  -1.213625750973335   1.084127989861559
    -0.461076791701681   0.262935199019578   0.847510043167997   2.623325372736294   0.559519886217974
    0.381698179386679  -0.758368998526128   0.528377669784949  -1.104500565300049   1.014107750781441
    0.223504652672365  -0.710739411300230   0.667004617269182  -1.266119743444892   0.840615170265153
    0.405483846857388  -0.619564541017890   0.672103138998045  -0.991288057676788   0.833750863929217
    0.293457659941662  -0.632627261500682   0.716704506632763  -1.136467051165794   0.771731111736407
    0.360002076632820  -0.543685626150342   0.758158588118322  -0.985926350691275   0.710311814173309
    0.349325036980205  -0.454336582204458   0.819481719509075  -0.915328405493848   0.610290230863546
    0.057866153383003  -0.947702095896818  -0.313866604985993  -1.509812605407607   1.890059027933620
    0.111284901653614  -0.131026358258335   0.985113071736186  -0.866690098000649   0.172765808971247
    0.126205560288092  -0.243179107041532  -0.961735971278423  -1.092078752779775   2.864065817548819
    0.212016500217405  -0.318004321096555  -0.924079139143116  -0.982764073685436   2.749415498004318
    0.293859071506476  -0.387843473783380  -0.873627086311426  -0.922401667287651   2.633403491923197
    0.248613255698018  -0.489115797754889  -0.836036593384422  -1.100537779602044   2.560815672769942
    0.545021902463028  -0.444858988621440  -0.710669829159991  -0.684554565908814   2.361246182214945
    0.372371404397101  -0.452712475716146  -0.810179579795917  -0.882465369471252   2.515254732266973
    0.550100816670943  -0.530285670446952  -0.645124948530584  -0.767059425296770   2.271983127431620
    0.557106860163998  -0.617562431269802  -0.555201395749635  -0.836818830413028   2.159401418237625
    0.642024939601504  -0.604416430990466  -0.471682896529472  -0.755234609059818   2.061994680140065
    0.558997478595967  -0.074433117017591   0.825821730166019  -0.132375971945543   0.599138531740776
    0.058296205796948  -0.716130726452358  -0.695527379058868  -1.489571009770332   2.339949967578350
    0.621586272160673  -0.732211324373456  -0.278382978503836  -0.866931628603243   1.852906451370166
    0.797774295454212   0.532064102078601   0.283661708363728   0.588188785594135   1.283185804715706
    0.654709287470920  -0.736793877709677  -0.168791382092980  -0.844319751080023   1.740399655004605
    0.628954306431283  -0.774770317373466  -0.064400588029629  -0.888904859984616   1.635241514329610
    0.591021897384986  -0.805787710644665   0.037406980437889  -0.937957991113809   1.533380617039016
    0.496748273411577  -0.863818906787821   0.084012196367621  -1.048916458180055   1.486684988180437
    0.581031026154992  -0.800456271677461   0.147216520056329  -0.942914976662174   1.423042789356838
    0.480049867050714  -0.854219514038332   0.199652565675243  -1.058805897432769   1.369792991852444
    0.647981695371081  -0.734417827727717   0.201866730244786  -0.847843030546120   1.367532810656050
    0.629452585705949  -0.710533490491583   0.314533942902738  -0.845833244706483   1.250830688689801
    0.713041622621644  -0.698620059906777  -0.059174794506992  -0.775182499678019   1.630005710799276
    0.528279713525195  -0.762479786915268   0.373557383574982  -0.964889061272640   1.187955242295637
    0.652639206451201  -0.585522827657270   0.480858694935406  -0.731244703929174   1.069162524006118
    0.557176252997457  -0.637403685848734   0.532232246638779  -0.852456998281161   1.009561215749371
    0.451739092015924  -0.677636220681170   0.580293843811538  -0.982812510048066   0.951706875482249
    0.513453647151006  -0.584860554920433   0.627935891254429  -0.850321798963606   0.891898154613298
    0.458944158118484  -0.522939131581471   0.718265218696919  -0.850481855232987   0.769490555508990
    0.472036704946927  -0.427685767740412   0.770886653960971  -0.736144012202942   0.690564363311399
    0.411909826937947  -0.358859015102286   0.837586116021372  -0.716678110400607   0.577946750105304
    0.292609940740739  -0.371618586239455   0.881066994581621  -0.903789880616613   0.492683001174784
    -0.396260656870134  -0.914599713066152   0.080528607812140  -1.979643505197770   1.490180427930583
    0.245801232049128  -0.206763618608001  -0.947011383428550  -0.699352543777270   2.814596294480980
    0.404694994520663  -0.340467854004705  -0.848707017644694  -0.699418573263409   2.584331970415491
    0.716654078334567   0.695789275834702  -0.047793468589900   0.770627103758533   1.618608009211025
    0.388645339123732   0.886992279882233  -0.249398267449311   1.157831293935944   1.822855165043836
    0.595823849909971  -0.346027698025628  -0.724747385009111  -0.526148581727866   2.381464018118063
    0.751220898559481  -0.346699290395710  -0.561664280159061  -0.432387970161347   2.167192298869582
    0.633296080707614  -0.528086161555086  -0.565739410094782  -0.695054502931684   2.172126011549267
    0.637145820766299  -0.671827620021989  -0.377747074185399  -0.811887384620369   1.958158207134540
    0.884326017401069  -0.383725376706717  -0.265936703406754  -0.409400340217460   1.839971812680638
    0.719091522578238  -0.587870833047628  -0.370560745098560  -0.685332222724815   1.950409000211603
    0.773413115064350  -0.574603754477249  -0.267698858397096  -0.638970887241852   1.841800254481284
    0.989258504500377   0.145784928625509  -0.010694197460899   0.146314758497716   1.581490728108135
    0.702213936534104  -0.659512074156001  -0.268215233310353  -0.754049806134566   1.842336229549582
    0.742693358454005  -0.649841771371684  -0.161592844794635  -0.718818623720052   1.733100824734397
    0.794451663128524  -0.603990116007870  -0.063580615891347  -0.650035471066270   1.634419858184701
    0.674052920727303  -0.736655460503696   0.054693624574344  -0.829745737237670   1.516075397098051
    0.790444992455612  -0.593899058452066   0.149935393658725  -0.644363104360546   1.420293399356658
    0.719396710327687  -0.680471331575627   0.139381275906756  -0.757598863835442   1.430959762727189
    0.720827853678411  -0.646255624624692   0.250521202700138  -0.730903761217185   1.317577738379412
    -0.570334131961711  -0.820934621177728  -0.028023662702464  -2.177981727103800   1.598823658744495
    0.757954565768660  -0.568195494208465   0.320404052083676  -0.643272724499418   1.244640331457908
    0.452232765142564   0.852402428265158   0.262479763827942   1.083017012822203   1.305205148769314
    0.739018535796594  -0.522814591902147   0.424872341113480  -0.615704812711497   1.131975465929536
    -0.822819045821922  -0.521360456685013   0.226168282563885  -2.576813611417373   1.342654097187847
    0.612090034084439  -0.537408570631230   0.580118796787871  -0.720520613101077   0.951921796994561
    0.554874762820276  -0.482223581317095   0.677919180438832  -0.715460062116939   0.825867921333714
    0.577991295520093  -0.379490303029280   0.722435583432703  -0.580978249919999   0.763477970168831
    -0.865012880322943  -0.369739163108326   0.339183826470781  -2.737659019500807   1.224747170180603
    0.450174001663372  -0.239602529165243   0.860194161944851  -0.489109049741280   0.535146041904053
    0.027337995568261   0.106964865865670   0.993886890681555   1.320573577881240   0.110628638630206
    0.341819221595257  -0.167340763377888   0.924746824087608  -0.455260180251511   0.390426503400185
    0.462665877049773  -0.696397504131932  -0.548608059048152  -0.984385142867961   2.151494811047443
    -0.113836940864978  -0.990386521881300   0.078585546828686  -1.685236043721488   1.492129667708950
    0.744444999430625   0.165983412806497  -0.646723394888297   0.219374397445840   2.274076975976882
    0.629098572433469  -0.241311355937905  -0.738920709993691  -0.366273953341240   2.402263461177029
    0.681225009493464  -0.628663514157330   0.375119544156242  -0.745292985778943   1.186270594289369
    0.688556590376797  -0.314131015522896  -0.653614203437497  -0.428011694184158   2.283146437154029
    0.638805780583645  -0.424487582700137  -0.641667723067254  -0.586502342755127   2.267467014470637
    0.705318364381223  -0.442843461189531  -0.553548257831486  -0.560655612870732   2.157415110571483
    0.717981419813809  -0.518136858268083  -0.464797673086072  -0.625113971344594   2.054202413972442
    0.790266722906618  -0.488868560192013  -0.369440167716266  -0.553992578882919   1.949202822169862
    0.835021270064543  -0.480877699160596  -0.267387578226454  -0.522500452984629   1.841477197546586
    0.553269021309184   0.536944412572526   0.636854840498777   0.770425509121244   0.880384380036225
    0.981356207271382  -0.144848386043715   0.126328696306339  -0.146542146253370   1.444129182232259
    0.879882907342024  -0.447365860567541  -0.160218151799386  -0.470375159486395   1.731707982632356
    0.818374409782960  -0.550955205747533  -0.163437103107236  -0.592533538313613   1.734969931092501
    -0.456679119727705  -0.847721965415189   0.269836341061277  -2.064932868487480   1.297573262908111
    -0.406681769086411  -0.833303259969586   0.374453756312856  -2.024826794194703   1.186988725993718
    0.831997988872182  -0.553011843489386   0.044240789698143  -0.586625474222438   1.526541092675404
    0.861477265985216  -0.488458363142339   0.138799667387727  -0.515801945975887   1.431547080169719
    0.815001653123669  -0.529053355363233   0.236378621251143  -0.575773160611431   1.332159173580917
    0.872115271427666  -0.424311431500401   0.243669371158398  -0.452814759211502   1.324648853883914
    0.917584612679007  -0.315235790291656   0.242208329943277  -0.330916807862843   1.326155017508334
    -0.359039700094521  -0.336300806013858   0.870627510265134  -2.388884651855653   0.514319866839814
    0.791519140370106  -0.417480079288817   0.446327048054195  -0.485359077890233   1.108139662640371
    0.703870453493298  -0.480622214048262   0.523037926026867  -0.599107691018492   1.020384911835398
    0.728174709590436  -0.382196322785476   0.568935464848295  -0.483344291730867   0.965585505168169
    0.645312518620371  -0.434830897112659   0.628087449506861  -0.592945412349335   0.891703400707350
    0.671848782533875  -0.326049833452703   0.665064447638885  -0.451820553718687   0.843216212096160
    0.516217594245071  -0.315180151745837   0.796354737122409  -0.548145699494105   0.649552173666629
    0.456265235582651  -0.863655931887268  -0.214290611354400  -1.084769668355351   1.786761839010415
    0.352826467933012  -0.274721030959280   0.894450579224244  -0.661573522989585   0.463595308605981
    0.225728989219266  -0.149712449956521   0.962617580222834  -0.585626404487989   0.274290639196337
    0.158442486388652  -0.133617515600274  -0.978285407246675  -0.700602428690087   2.932816992773038
    0.366060194453126  -0.166776824026930  -0.915524671980187  -0.427500664726474   2.727606622330636
    0.647007119538881   0.722366169692383   0.244067417223100   0.840374490013857   1.324238416013122
    0.648590701153767  -0.128712507546756  -0.750175441332154  -0.195904261579936   2.419123688040788
    0.571912053724642  -0.067039525863283  -0.817570978433364  -0.116687486294974   2.527976312236451
    0.731534086579975  -0.079279714064659  -0.677179892723807  -0.107953276671947   2.314719551015486
    0.794524485199531  -0.134758427571453  -0.592090372001512  -0.168010022336914   2.204446630365240
    0.420808947352431  -0.802151279921787   0.423642719691933  -1.087662571363020   1.133333351209781
    0.846857041889136  -0.200754990616607  -0.492473942808457  -0.232762202402267   2.085726344569327
    0.824160218982307  -0.310347812477497  -0.473755389138181  -0.360139807857907   2.064346544162629
    0.179847517528647  -0.356265547373878   0.916913153026596  -1.103305650790169   0.410520673924263
    0.848653587398230  -0.376561559006239  -0.371468007875383  -0.417616254371438   1.951385993892177
    0.917714330778816  -0.281169330031151  -0.280613996324175  -0.297299692303686   1.855230075217194
    0.927395158627017  -0.334791523210033  -0.166891748573367  -0.346442352249260   1.738472685138114
    0.918328019397424  -0.390826433834301  -0.062676529945160  -0.402365862360771   1.633513965315514
    0.295862599959849  -0.828246916640390  -0.475891129378019  -1.227708001589946   2.066773298972467
    0.315805871953071   0.623511573097670  -0.715192260481914   1.101965746284011   2.367695473795597
    0.943566192141133  -0.329112144723757   0.036986987490338  -0.335601840638510   1.533800900845933
    0.914551426418956  -0.379223737126972   0.140659324733575  -0.393076201044626   1.429668998009423
    0.491228548501539  -0.726228016082314   0.480923467710006  -0.976078017365584   1.069088648291147
    0.280959664854905  -0.096306212261916  -0.954875269448519  -0.330224688527275   2.840035865993076
    0.859096675394559  -0.286194532218602   0.424317796061149  -0.321571374546779   1.132587963454342
    -0.635181199670441  -0.754339014855147   0.165883978287792  -2.270649129957631   1.404141985308421
    0.888332604274813  -0.169357279612939   0.426829352346383  -0.188385605596408   1.129812528052092
    0.827172874620718  -0.211896041657764   0.520466236197381  -0.250776292950157   1.023399447592813
    0.753284393201340  -0.266583785218477   0.601245131720730  -0.340140771372600   0.925737893248752
    0.699094248140022  -0.205673505480457   0.684810660957405  -0.286127206818998   0.816452483438938
    0.234100797658886   0.969809099120891   0.068316380157320   1.333938849092799   1.502426694508256
    0.636283385860250  -0.139933144641428   0.758658136388833  -0.216476549080506   0.709545376763695
    0.548670893437261  -0.195352046718252   0.812894721681564  -0.342050775506403   0.621691071011825
    0.465938093457528  -0.130775159694791   0.875099737556792  -0.273630512441888   0.505154454901390
    0.284912859013201  -0.050618827028598   0.957215961588174  -0.175829565525780   0.293573318536697
    0.191116817335276  -0.024516436544098  -0.981261079667799  -0.127583070410781   2.947697119459968
    0.393275798262152  -0.051378020501649  -0.917983902642415  -0.129905491806517   2.733763303023880
    0.489533650408080   0.001364213940548  -0.871983339312435   0.002786755124826   2.630035615201891
    0.475817528991103  -0.119835727782629  -0.871342112749986  -0.246721220780002   2.628727272161526
    0.738376829935112   0.042427600004414  -0.673052416809300   0.057397518212914   2.309124551985869
    0.803339419127078  -0.021270897564229  -0.595141434109053  -0.026471909896895   2.208237960067573
    -0.019304705558076  -0.652977583912654  -0.757131166477057  -1.600351834731021   2.429706648989152
    0.860388643288889  -0.088401027697105  -0.501912981304130  -0.102386185205166   2.096605435076942
    0.913100957091991  -0.028272186842317  -0.406752167306878  -0.030952937897305   1.989692340724356
    0.939653849403507  -0.173694209317034  -0.294755771700745  -0.182785898075433   1.869996252196322
    0.554166368326900  -0.184732760150077  -0.811648596094830  -0.321767830042397   2.517765166655911
    0.840814622420000   0.312319885333041  -0.442139186173618   0.355653976843433   2.028778567467865
    0.035601062981032   0.368206365888704  -0.929062234962575   1.474408145578756   2.762666022309057
    0.975048363839868  -0.109200505990911  -0.193276842028541  -0.111530210039118   1.765297198329524
    0.957963211050493  -0.222600405364544  -0.180984932536850  -0.228316653945505   1.752784156902564
    0.983829836643833  -0.159691937164877  -0.081100787504796  -0.160913240543671   1.651986283019373
    0.999645678328188   0.026300122280570   0.004101386083230   0.026303376456167   1.566694929213092
    0.995493699469522  -0.093241433871689   0.017272212540279  -0.093391043810085   1.553523255338050
    0.993059198193241  -0.025045423993646   0.114918038709230  -0.025215128998568   1.455623835428758
    0.953063322229011  -0.195462036827395   0.231224773718316  -0.202283233851808   1.337459942704185
    0.429489346053000   0.591483784657969  -0.682411777531464   0.942753137465898   2.321853331306474
    0.971774922110781  -0.075827122125341   0.223391468742147  -0.077871723565443   1.345503838200919
    0.911725141682115  -0.243072135439603   0.331169447560549  -0.260546512045256   1.233253636355551
    0.904796074161305  -0.052432065796455   0.422604948691573  -0.057884306424771   1.134478703070446
    0.848702870570814  -0.097830692799038   0.519742814314273  -0.114764323434753   1.024246444093209
    0.795157221811587  -0.033123343568293   0.605498007190557  -0.041632275325127   0.920404850625503
    0.718667090105937  -0.085416706283823   0.690088110233968  -0.118299381678474   0.809185535712338
    0.730231378256247   0.030284535007144   0.682528373878612   0.041448766754389   0.819579800887150
    0.553010694571918  -0.690844367455758  -0.465739446089483  -0.895759598352134   2.055266360493745
    -0.575284995866519  -0.707638984753100   0.410236810620964  -2.253392289778706   1.148082613030831
    0.453792475757916   0.110667932276559   0.884208684480724   0.239204229348582   0.485999173629562
    0.968368987115544   0.126580163259934  -0.215032479086077   0.129977865202951   1.787521413587213];

UDSR.tri = [490         512         546
    388         207         387
    395         365         393
    420         421         450
    233         265          64
    233         507         265
    604         603         605
    988         917         990
    987         990         917
    922         746         923
    955         923         746
    469         470         439
    469         468         437
    469         439         468
    953         129         626
    536         363         537
    567         538         363
    537         363         538
    536         537         589
    972         941         979
    827         344          94
    157          94         156
    699         729         726
    692         690         691
    721          33         719
    721         689          33
    721         691         689
    690         689         691
    545         546         511
    545         511         542
    510         511         512
    546         512         511
    510         512         514
    207         386         387
    395         396         365
    420         450         206
    420         206         481
    421         417         450
    388         387         417
    416         417         387
    416         450         417
    356         357         608
    357         389         608
    395         393         392
    395         392         221
    359         221         392
    388         390         389
    388         417         390
    421         359         390
    421         390         417
    854         823         112
    110         112          80
    823          80         112
    50          81         823
    77          80          81
    823          81          80
    953         197         161
    192         193         197
    162         161         197
    162         197         193
    100         131         132
    133         132         131
    265          34          64
    180         179         175
    895         145         146
    175         146         145
    530         234         964
    243         274         275
    304         305         274
    275         274         305
    180         181         179
    243         275         621
    245         621         275
    604         573         603
    603         573         602
    603         602         633
    528         521          31
    616         802         646
    336         335         337
    745         342         311
    374         311         342
    434         587         496
    434         496         358
    497         358         496
    619         620         590
    650         590         620
    10         804         977
    978         804         982
    10         236         804
    236         982         804
    978         104         804
    977         804         104
    10         977         119
    39         980          41
    42          41         980
    985          17         987
    988          97         917
    985         987         986
    987         917         986
    896         928         929
    930         929         928
    931         898         929
    931         929         930
    962         932         930
    931         930         932
    699         731         729
    772         732         906
    701         906         732
    656         700         698
    162         109         161
    133         131         109
    161         109         131
    133         134         132
    877         309         165
    877         144         598
    877         165         200
    877         200         144
    233         232         507
    233         234         232
    498         747         407
    499         747         593
    498         654         747
    654         593         747
    436         345         376
    286         348         316
    439         502         468
    319         318         316
    953         626         191
    953         191         197
    192         197         191
    158         626         128
    129         128         626
    537         538         539
    567         568         538
    567         600         568
    536         589         504
    536         504         534
    474         472         589
    470         504         472
    589         472         504
    571         404         444
    604         404         573
    571         573         404
    541         444         543
    541         543         542
    545         542         543
    908         979         940
    941         940         979
    908         940          57
    942          57         940
    844         875         847
    718         686         685
    811         779         716
    685         716         779
    811         812         779
    780         779         812
    718         685         717
    685         779         717
    780         717         779
    780         749         717
    1         160         716
    811         716         160
    466         810         513
    466         156         810
    344         810          94
    156          94         810
    827         634         344
    827           2         634
    1         582         594
    562         594         582
    562         593         594
    654         594         593
    596         627         597
    628         659         519
    628         630         659
    597         627         659
    519         659         627
    562         531         593
    499         593         531
    729         815         726
    471         816         814
    757         787         759
    721         719         722
    628         519         658
    690         628         658
    690         658         689
    510         514         479
    542         511         509
    510         509         511
    581         321         320
    581         320         352
    581         352         353
    581         353         354
    207         354         353
    384         353         352
    384         386         353
    207         353         386
    357         355         389
    207         355         354
    388         355         207
    388         389         355
    581         323         321
    581         354         323
    490         482         512
    733         482         453
    733         514         482
    512         482         514
    733         452         481
    733         453         452
    125         480         206
    125         479         480
    514         480         479
    733         480         514
    733         481         480
    481         206         480
    389         391         608
    389         390         391
    359         392         391
    359         391         390
    254         484         452
    420         481         484
    481         452         484
    110         111         112
    110         142         111
    77          78          80
    110          80          78
    77          81          44
    173         143         142
    173         174         143
    939         111         143
    142         143         111
    110         876         142
    110          78         876
    107         876          78
    895         214         174
    895         146         214
    939         143         214
    174         214         143
    598         171         169
    530         964         415
    204         415         964
    953         161         989
    161         131         989
    953         989         129
    100         989         131
    48         980          73
    39          73         980
    100         132         101
    269         270         267
    64         267         238
    269         267          34
    64          34         267
    265         264          34
    151         150         181
    179         181         150
    146         147         214
    939         147         113
    939         214         147
    895         208         145
    180         175         210
    175         145         210
    180         210         211
    218         155         185
    151         181         182
    234         235         964
    233         235         234
    233          64         235
    64         238         235
    507         230         262
    199         230         231
    507         232         230
    231         230         232
    260         262         273
    199         273         230
    262         230         273
    418          82         299
    304         302         301
    304         274         302
    304         713         305
    336         713         335
    336         306         713
    305         713         306
    245         246         276
    245         276         621
    553         276         216
    246         216         276
    553         182         213
    181         213         182
    180         213         181
    180         211         213
    243         621         244
    693         491         695
    664         636         633
    603         636         605
    603         633         636
    602         670         633
    897         610         609
    360         609         610
    604         605         315
    776         554         761
    528          31         487
    205         487          31
    521         522          31
    521         761         522
    616         646         688
    897         612         610
    579         610         612
    897         644         612
    901         612         644
    776         761         552
    521         552         761
    490         546         515
    698         700         873
    643         641         998
    675         998         641
    434         464         587
    401         400         429
    335         368         337
    350         368         508
    367         368         350
    367         337         368
    336         337         338
    340         341         178
    745         280         219
    745         311         280
    619         590         588
    619         588         557
    650         620         649
    648         799         649
    646         802         618
    648         618         802
    648         649         618
    620         618         649
    493         524         494
    524         526         494
    236         750         982
    985         750          14
    236          14         750
    985         986         750
    985          14         885
    985         885          17
    18          17         885
    10         119          11
    42          11         119
    46          13          45
    48          75          45
    48          43         980
    42         980          43
    48          45          43
    13          43          45
    42          43          11
    13          11          43
    42         119         292
    42         292          41
    917          97          28
    917          28         986
    896         929         865
    898         865         929
    958         203         928
    930         928         203
    962         930         203
    172         963         965
    962         963         932
    172         932         963
    965         963         995
    172         965         966
    967         902         966
    868         869         836
    172         966         933
    902         933         966
    818         820         674
    922         890         202
    891         202         890
    922         923         890
    923         924         890
    891         890         892
    924         892         890
    955         925         923
    923         925         924
    967         904         902
    935         902         904
    903          32          63
    62         296         999
    999         296          30
    113         147         148
    854         114         116
    113         148         114
    988         990         756
    988         756         991
    965         995          19
    729         731         697
    89         769         767
    89         736         769
    772         766         732
    735         732         766
    735         766         769
    767         769         766
    772         763         762
    699         763         731
    772         906         763
    906         731         763
    748         762         793
    906         730         731
    701         730         906
    701         700         730
    731         730         697
    656         730         700
    656         697         730
    132         134         102
    132         102         101
    72         101         102
    166         309         167
    598         169         167
    877         598         167
    877         167         309
    136         166         167
    136         167         169
    165         309         163
    162         163         109
    530         144         201
    144         200         201
    231         201         200
    530         201         234
    234         201         232
    231         232         201
    499         435         747
    407         747         435
    436         435         467
    499         467         435
    436         376         435
    407         435         376
    436         378         345
    345         378         176
    319         316         317
    316         348         317
    439         440         502
    319         287         318
    319         320         287
    321         287         320
    192         191         225
    129         331         128
    158         128         127
    158         127         189
    827         251          95
    827          94         251
    95         251         127
    189         127         251
    157         251          94
    157         189         251
    941         312         940
    942         940         312
    224         257         569
    286         572         284
    7         284         572
    316         318         477
    255         572         477
    286         316         477
    286         477         572
    165         451         200
    165         198         451
    231         200         451
    199         231         451
    199         451         198
    537         539         505
    474         505         506
    474         589         505
    537         505         589
    538         568         525
    541         525         444
    571         444         525
    571         525         568
    541         539         525
    538         525         539
    602         573         570
    600         570         568
    571         570         573
    571         568         570
    596         597         737
    536         534         535
    536         535         363
    474         473         472
    545         547         546
    546         547         515
    548         515         547
    548         547         578
    444         404         574
    444         574         543
    977         104         976
    104         945         976
    194         951         950
    908         910         878
    908          57         910
    879         910          96
    57          96         910
    913         880         912
    818         905         997
    471         905         816
    844          40         875
    908         878         130
    908         130         979
    685         655         716
    624         655         686
    685         686         655
    1         716         655
    1         655         582
    624         582         655
    718         687         686
    519         687         658
    844         845         812
    844         847         845
    813         845         847
    780         812         845
    780         781         749
    813         781         845
    780         845         781
    33         852         719
    718         717         975
    717         749         975
    749         852         975
    33         975         852
    654         684         594
    1         684         160
    1         594         684
    597         659         740
    630         740         659
    596         625         627
    596         516         625
    562         582         751
    516         751         625
    624         751         582
    624         625         751
    563         425         516
    562         425         531
    562         751         425
    516         425         751
    534         503         501
    469         503         470
    534         504         503
    470         503         504
    563         501         532
    563         532         425
    531         425         532
    563         533         501
    534         501         533
    534         533         535
    565         535         533
    757         726         840
    726         815         840
    813         847         285
    813         285         814
    848         285         846
    848         471         285
    471         814         285
    36         952         785
    759         787         790
    542         509         540
    506         540         509
    541         542         540
    541         540         539
    506         505         540
    539         540         505
    321         323         322
    291         322         325
    354         108         323
    357         108         355
    354         355         108
    401         429         431
    427         458         459
    423         458         426
    427         426         458
    427         424         426
    396         426         424
    254         423         422
    254         422         484
    399         424         428
    399         428         400
    429         400         428
    427         428         424
    365         396         397
    396         424         397
    399         397         424
    365         397         350
    367         397         399
    367         350         397
    356         447         330
    356         608         447
    418         299          65
    298          65         299
    298         327          65
    330          65         327
    418         332         508
    418          65         332
    330         332          65
    362         332         447
    330         447         332
    350         508         364
    365         350         364
    362         364         332
    508         332         364
    393         365         364
    393         364         362
    490         483         482
    453         482         483
    608         391         361
    362         447         361
    608         361         447
    392         361         391
    393         362         361
    393         361         392
    359         595         221
    421         595         359
    221         595         422
    484         422         595
    420         595         421
    420         484         595
    112         111         177
    854         112         177
    939         113         177
    939         177         111
    854         177         114
    113         114         177
    107          78          79
    77          79          78
    50         796          49
    18          49         796
    50          49          81
    81          49          44
    142         876         140
    173         142         140
    173         140         171
    169         171         140
    598         144         170
    598         170         171
    530         170         144
    530         415         170
    75         141          79
    107         141         106
    107          79         141
    72          73          71
    72          71         101
    39          71          73
    270         268         267
    239         268         379
    239         238         268
    238         267         268
    34         264         266
    269          34         266
    269         266         298
    295         266         264
    295         327         266
    298         266         327
    265         507         263
    265         263         264
    507         262         263
    204         419         415
    204         208         419
    895         419         208
    895         174         419
    204         297         208
    218         248         155
    151         182         183
    224         258         257
    260         273         258
    301         270         300
    301         300          82
    299          82         300
    269         300         270
    269         298         300
    298         299         300
    304         449         713
    335         713         449
    304         301         449
    301          82         449
    379         268         271
    270         271         268
    301         271         270
    301         302         271
    243         272         274
    274         272         302
    379         271         272
    302         272         271
    340         339         402
    336         339         306
    336         338         339
    402         339         338
    340         178         308
    340         308         339
    306         339         308
    185         183         215
    182         215         183
    218         185         215
    218         215         216
    553         215         182
    553         216         215
    745         219         279
    219         248         279
    246         247         216
    246         278         247
    278         279         247
    218         216         247
    218         247         248
    248         247         279
    621         241         244
    553         213         241
    211         241         213
    211         244         241
    553         241         276
    621         276         241
    211         210         240
    211         240         244
    692         691         377
    692         377         693
    693         377         491
    695         491         696
    656         696         697
    605         636         398
    893         639         398
    602         601         670
    602         570         601
    662         601         600
    600         601         570
    662         740         631
    630         631         740
    662         631         601
    670         601         631
    664         633         632
    633         670         632
    360         638         609
    959         638         639
    761         554         523
    493         523         524
    556         524         523
    556         523         554
    493         522         523
    761         523         522
    619         557         555
    556         555         557
    556         554         555
    459         458         457
    205         457         487
    205         459         457
    487         457         455
    616         688         583
    616         583         614
    644         614         583
    901         644         583
    893         249         639
    959         639         249
    643         611         641
    959         641         611
    959         611         638
    609         638         611
    648         677         799
    676         799         677
    741         770         771
    89         770         736
    89         768         770
    771         770         768
    738         706         736
    735         769         706
    736         706         769
    676         707         709
    738         709         707
    650         649         679
    799         679         649
    741         771         773
    335         333         368
    335         449         333
    82         333         449
    418         333          82
    418         508         333
    508         368         333
    337         369         338
    367         369         337
    367         399         369
    399         400         369
    340         402         372
    340         372         341
    402         403         372
    405         372         403
    463         464         833
    434         833         464
    405         833         433
    374         433         358
    434         358         433
    434         433         833
    342         310         341
    278         310         279
    178         310         278
    178         341         310
    745         310         342
    745         279         310
    311         343         280
    903         280         343
    622         591         559
    903         671          32
    903         343         671
    465         671         343
    622         588         927
    590         927         588
    463         495         464
    494         526         495
    587         464         495
    587         495         526
    493         494         492
    587         526         527
    587         527         496
    559         496         527
    986         984         750
    194         984          28
    986          28         984
    982         750         984
    10        1000         236
    10          11        1000
    236        1000          14
    13        1000          11
    46          76          44
    46          45          76
    77          44          76
    75          76          45
    77          76          79
    75          79          76
    977           9         119
    119           9         292
    977         976           9
    974           9         976
    194          28         954
    97         954          28
    955         954          97
    194         954         951
    955         746         954
    746         951         954
    868         866         898
    898         866         865
    217         865         866
    217         866         956
    924         894         892
    868         836         867
    868         867         866
    956         866         867
    29          19         996
    995         996          19
    870         900         933
    931         932         900
    172         900         932
    172         933         900
    818         997         819
    818         819         820
    922         202         888
    853         881         886
    853         819         881
    997         881         819
    955         926         925
    955          97         926
    988         991         926
    988         926          97
    968         296          93
    62          93         296
    968         366         296
    967         366         904
    968         904         366
    967          30         366
    30         296         366
    147         149         148
    175         179         149
    117         148         149
    175         149         146
    146         149         147
    117         149         150
    179         150         149
    151         916         150
    117         916         118
    117         150         916
    118          85         115
    117         118         115
    117         115         148
    148         115         114
    116         114         115
    116         115          85
    18          20          17
    18         796          20
    990         758         756
    962         203         961
    958         961         203
    29          60          61
    29          61          19
    999          30          61
    999          61          60
    767         766         765
    748         765         762
    772         762         765
    772         765         766
    748         795         828
    748         793         795
    792         795         793
    832           8         798
    828         798           8
    835         803         801
    771         768         803
    792         793         764
    134         135         102
    136         106         135
    136         135         166
    166         135         134
    309         394         163
    166         394         309
    133         109         394
    109         163         394
    133         394         134
    166         134         394
    165         196         198
    165         163         196
    162         196         163
    198         196         228
    162         193         196
    193         228         196
    157         981         189
    157         156         188
    157         188         981
    313         981         188
    176         347         346
    348         346         347
    286         346         348
    286         284         346
    284         314         346
    345         176         314
    176         346         314
    7         283         284
    284         283         314
    468         502         229
    436         467         408
    436         408         378
    437         408         467
    378         408         229
    437         468         408
    468         229         408
    176         378         410
    176         410         347
    349         347         410
    378         229         410
    319         351         320
    319         317         351
    352         320         351
    126         440         441
    126         441         473
    472         473         441
    470         472         441
    470         441         439
    439         441         440
    502         440         668
    224         259         258
    291         260         259
    260         258         259
    192         225         223
    255         569         223
    626         190         191
    191         190         225
    158         190         626
    39          70          71
    101          71          70
    974         942         973
    942         312         973
    225         190         282
    255         253         572
    255         223         253
    225         253         223
    225         282         253
    7         572         253
    7         253         282
    318         256         477
    224         569         256
    255         256         569
    255         477         256
    565         737         599
    662         600         599
    662         599         740
    597         599         737
    597         740         599
    563         516         564
    563         564         533
    596         564         516
    596         737         564
    565         564         737
    565         533         564
    474         506         475
    474         475         473
    404         575         574
    604         575         404
    604         315         575
    543         574         544
    545         543         544
    545         544         547
    547         544         578
    978         982         983
    982         984         983
    194         950         983
    194         983         984
    879          96         946
    945         946          96
    948         949         882
    978         983         949
    950         949         983
    950         882         949
    879         912         911
    879         911         910
    878         910         911
    878         911         846
    912         880         911
    848         846         911
    848         911         880
    942         943          57
    57         943          96
    974         943         942
    974         976         943
    945         943         976
    945          96         943
    816         905         817
    36         817         952
    818         817         905
    818         674         817
    674         952         817
    2         970         634
    811         843         812
    844         812         843
    844         843          40
    878         846          83
    878          83         130
    875         130          83
    847         875          83
    847          83         285
    846         285          83
    972         979         227
    979         130         227
    875         227         130
    875          40         227
    686         687         657
    627         625         657
    519         627         657
    519         657         687
    624         686         657
    624         657         625
    749         781         673
    749         673         852
    658         687         635
    33         635         975
    718         635         687
    718         975         635
    33         689         635
    689         658         635
    654         561         684
    466         513         561
    466         561         498
    498         561         654
    501         294         532
    469         437         294
    469         294         503
    501         503         294
    499         531         500
    531         532         500
    499         500         467
    437         467         500
    437         500         294
    532         294         500
    757         755         787
    757         840         755
    754         755         840
    754         785         755
    719         753         722
    36         785         753
    754         722         753
    754         753         785
    787         788         790
    674         820         788
    934         788         820
    934         790         788
    325         322         324
    323         324         322
    323         108         324
    321         290         287
    321         322         290
    291         290         322
    291         259         290
    423         426         329
    423         329         422
    395         329         396
    396         329         426
    395         221         329
    221         422         329
    295         328         327
    356         330         328
    330         327         328
    901         580         612
    901         550         580
    551         580         550
    521         520         552
    528         520         521
    552         520         550
    551         550         520
    528         518         520
    551         520         518
    18         885          47
    18          47          49
    46          44          47
    44          49          47
    876         137         140
    136         137         106
    136         169         137
    169         140         137
    107         137         876
    107         106         137
    173         171         411
    171         170         411
    173         411         174
    174         411         419
    415         411         170
    415         419         411
    75          74         141
    72          74          73
    48          74          75
    48          73          74
    264         263         293
    295         264         293
    964         235         237
    204         964         237
    204         237         297
    239         297         237
    239         237         238
    238         237         235
    219         280         187
    903         187         280
    903          63         187
    151         183         152
    151         152         916
    257         258         164
    199         198         164
    199         164         273
    273         164         258
    198         228         164
    257         164         228
    178         277         308
    245         277         246
    178         278         277
    246         277         278
    305         306         307
    306         308         307
    275         305         307
    245         275         307
    245         307         277
    308         277         307
    210         209         240
    145         209         210
    145         208         209
    208         297         209
    240         209         222
    239         379         222
    239         222         297
    297         222         209
    244         240         242
    243         244         242
    243         242         272
    240         222         242
    379         272         242
    379         242         222
    664         637         636
    664         665         637
    667         637         665
    893         637         667
    893         398         637
    636         637         398
    698         667         666
    667         665         666
    695         666         665
    695         696         666
    656         698         666
    656         666         696
    691         184         377
    721         184         691
    721         722         184
    697         696         728
    729         728         815
    729         697         728
    491         728         696
    664         663         665
    664         632         663
    695         665         663
    693         695         663
    670         631         712
    670         712         632
    630         712         631
    639         607         398
    639         638         607
    360         607         638
    605         398         607
    605         607         315
    360         315         607
    619         555         617
    646         618         617
    619         617         620
    620         617         618
    646         585         688
    646         617         585
    555         585         617
    554         585         555
    776         585         554
    776         688         585
    453         483         485
    488         485         483
    487         455         486
    528         487         486
    488         486         485
    455         485         486
    528         486         518
    488         518         486
    458         456         457
    457         456         455
    423         456         458
    254         456         423
    616         615         802
    616         614         615
    643         615         614
    643         998         615
    959         642         641
    959         249         642
    249         672         642
    675         641         642
    893         669         249
    893         667         669
    873         672         669
    249         669         672
    698         873         669
    698         669         667
    643         614         613
    643         613         611
    897         613         644
    644         613         614
    897         609         613
    609         611         613
    648         802         645
    648         645         677
    802         615         645
    998         677         645
    998         645         615
    675         640         998
    998         640         677
    676         640         707
    676         677         640
    738         707         705
    738         705         706
    738         736         708
    738         708         709
    739         709         708
    741         739         708
    741         708         770
    736         770         708
    799         678         679
    739         678         709
    676         678         799
    676         709         678
    835         837          16
    489         774          16
    870         871         839
    870         933         871
    935         871         902
    902         871         933
    489          16         806
    837         806          16
    32         841         704
    32         704          63
    968          93         704
    63         704          93
    303         936         841
    968         936         904
    968         704         936
    841         936         704
    303         947          12
    303          12         936
    935         904          12
    904         936          12
    935         872         871
    839         871         872
    935          12         872
    947         872          12
    338         369         370
    401         403         370
    401         370         400
    400         370         369
    402         338         370
    402         370         403
    401         432         403
    401         431         432
    463         432         431
    463         833         432
    405         403         432
    405         432         833
    405         433         373
    405         373         372
    341         372         373
    342         341         373
    374         342         373
    374         373         433
    311         375         343
    374         375         311
    374         358         375
    465         343         375
    497         375         358
    497         465         375
    559         591         529
    497         496         529
    559         529         496
    465         448         671
    303         809         947
    744         809         715
    744         808         809
    947         809         808
    493         492         460
    493         460         522
    31         522         460
    205          31         460
    205          59         459
    205         460          59
    492          59         460
    492         461          59
    429         461         431
    494         495         462
    494         462         492
    463         462         495
    492         462         461
    463         431         462
    431         461         462
    559         527         558
    557         588         558
    622         559         558
    622         558         588
    526         212         527
    556         557         212
    557         558         212
    527         212         558
    556         212         524
    524         212         526
    924         925         138
    924         138         894
    896         138         928
    896         894         138
    892         894         909
    217         909         865
    896         865         909
    896         909         894
    217         956         864
    860         859         830
    860         858         859
    857         859         858
    832         834         801
    956         867         834
    931         900         899
    868         899         869
    870         869         899
    870         899         900
    868         898         899
    931         899         898
    857         858         606
    857         606         856
    934         820         821
    853         821         819
    820         819         821
    887         888         863
    202         863         888
    886         881         883
    913         850         880
    913         883         850
    881         850         883
    997         850         881
    925         926          35
    925          35         138
    958          35         991
    991          35         926
    958         928          35
    928         138          35
    22         409          20
    990         409         758
    22         758         409
    987         409         990
    987          17         409
    17          20         409
    22          20          51
    796          51          20
    86          55          85
    116          85          55
    50         823          52
    854          52         823
    29          58          60
    118          88          85
    86          85          88
    121         120         371
    121         371         152
    916         152         371
    120          88         371
    118         916         371
    118         371          88
    19          61         851
    965          19         851
    967         851          30
    30         851          61
    967         966         851
    965         851         966
    828         829         798
    828         795         829
    792         826         795
    795         826         829
    860         826         858
    860         829         826
    828           8         797
    767         765         797
    748         828         797
    748         797         765
    89         767         720
    89         720         768
    767         797         720
    8         720         797
    759         764         727
    757         759         727
    699         726         727
    757         727         726
    102         135         105
    72         102         105
    106         141         105
    106         105         135
    72         105          74
    141          74         105
    407         376         250
    313         188         250
    416         414         446
    125         206         446
    416         446         450
    450         446         206
    416         586         414
    416         387         586
    387         386         586
    313         252         981
    313         283         252
    7         252         283
    7         282         252
    502         380         229
    502         668         380
    349         410         380
    229         380         410
    349         380         381
    668         381         380
    317         382         351
    349         382         347
    349         381         382
    348         347         382
    348         382         317
    443         414         412
    414         586         412
    318         287         289
    224         289         259
    287         290         289
    259         289         290
    224         256         289
    318         289         256
    192         223         226
    569         226         223
    192         226         193
    193         226         228
    257         226         569
    257         228         226
    95         127          66
    128          66         127
    128         331          66
    331         159          66
    941          56         312
    39         723          70
    39          41         723
    41         907         723
    331          98         159
    974         973           6
    907           6         973
    974           6           9
    292           9           6
    41           6         907
    41         292           6
    565         599         566
    567         363         566
    567         566         600
    600         566         599
    363         535         566
    565         566         535
    506         509         478
    506         478         475
    510         478         509
    510         479         478
    125         476         479
    479         476         478
    475         478         476
    473         475         334
    126         334         443
    126         473         334
    475         476         334
    574         575         576
    578         544         576
    574         576         544
    579         578         576
    879          37         912
    879         946          37
    913         912          37
    948          37         946
    104         800         945
    945         800         946
    978         800         104
    978         949         800
    948         946         800
    948         800         949
    951         918         950
    950         918         882
    513         969         937
    344         969         810
    513         810         969
    344         634         969
    513         937         778
    160         684         778
    513         778         561
    684         561         778
    811         842         843
    811         160         842
    160         778         842
    937         842         778
    36         783         817
    816         817         783
    781         782         673
    813         814         782
    813         782         781
    673         782         783
    814         816         782
    816         783         782
    852         673         784
    719         852         784
    36         784         783
    673         783         784
    36         753         784
    719         784         753
    787         755         786
    787         786         788
    952         786         785
    785         786         755
    674         786         952
    674         788         786
    108         288         324
    356         288         357
    357         288         108
    356         328         288
    295         326         328
    325         324         326
    324         288         326
    328         326         288
    295         293         326
    325         326         293
    612         580         623
    579         612         623
    579         623         578
    548         578         623
    551         518         517
    488         517         518
    490         515         517
    490         517         483
    488         483         517
    552         550         584
    776         552         584
    901         583         584
    901         584         550
    776         584         688
    688         584         583
    14        1000          15
    14          15         885
    885          15          47
    13          15        1000
    46          15          13
    46          47          15
    263         261         293
    262         261         263
    291         325         261
    325         293         261
    291         261         260
    260         261         262
    999          60          91
    185         153         183
    183         153         152
    121         153         122
    121         152         153
    219         187         186
    155         186          92
    219         186         248
    155         248         186
    62          92         124
    62         124          93
    92         186         124
    187         124         186
    63          93         124
    63         124         187
    491         377         694
    377         184         694
    815         728         694
    491         694         728
    630         629         712
    692         629         690
    690         629         628
    628         629         630
    693         663         661
    692         693         661
    632         661         663
    632         712         661
    692         661         629
    712         629         661
    453         454         452
    453         485         454
    254         452         454
    254         454         456
    455         454         485
    455         456         454
    707         921         705
    675         921         640
    707         640         921
    672         920         642
    675         642         920
    675         920         921
    873         702         672
    701         702         700
    700         702         873
    735         706         734
    706         705         734
    701         734         702
    701         732         734
    735         734         732
    927         652         592
    622         592         591
    622         927         592
    682         592         652
    653         592         682
    653         591         592
    741         710         739
    741         773         710
    650         679         680
    682         652         714
    744         715         714
    682         714         715
    489         775         774
    744         775         808
    489         808         775
    773         805         774
    774         805          16
    771         805         773
    835          16         805
    835         805         803
    771         803         805
    837         838         806
    836         869         838
    836         838         837
    839         806         838
    870         838         869
    870         839         838
    489         806         807
    489         807         808
    839         807         806
    839         872         807
    947         808         807
    947         807         872
    653         560         591
    591         560         529
    653         448         560
    465         560         448
    497         560         465
    497         529         560
    303         841         777
    303         777         809
    715         809         777
    459          59         430
    59         461         430
    427         459         430
    429         430         461
    427         430         428
    429         428         430
    892         909         862
    217         862         909
    217         864         862
    830         862         864
    956         944         864
    956         834         944
    832         798         944
    832         944         834
    857         889         859
    891         889         202
    202         889         863
    857         856         889
    856         863         889
    801         834         647
    835         801         647
    836         647         867
    867         647         834
    835         647         837
    836         837         647
    86         794          55
    54          55         794
    756         758         992
    792         764         789
    759         790         789
    759         789         764
    853         855         821
    853         886         855
    887         855         886
    887         863         855
    856         855         863
    856         606         822
    856         822         855
    821         855         822
    934         821         822
    887         886         884
    886         883         884
    997         849         850
    880         850         849
    471         849         905
    997         905         849
    848         849         471
    848         880         849
    54          23          51
    54         794          23
    25          23         794
    22          51          23
    54          51          53
    796          53          51
    50          53         796
    50          52          53
    54          84          55
    54          53          84
    52          84          53
    116          55          84
    854         116          84
    854          84          52
    120          58          87
    86          88          87
    120          87          88
    768         720         660
    801         803         660
    768         660         803
    832         801         660
    832         660           8
    8         660         720
    764         760         727
    762         763         760
    699         760         763
    699         727         760
    762         760         793
    793         760         764
    407         250         281
    466         281         156
    156         281         188
    188         281         250
    466         498         281
    498         407         281
    313         250         220
    313         220         283
    314         283         220
    345         314         220
    345         220         376
    376         220         250
    282         103         252
    158         103         190
    190         103         282
    158         189         103
    189         981         103
    981         252         103
    351         382         383
    384         383         752
    752         383         381
    381         383         382
    384         352         383
    352         351         383
    586         385         412
    386         385         586
    384         385         386
    384         752         385
    440         442         668
    126         442         440
    126         443         442
    443         412         442
    159          38          66
    2         971         970
    941         971          56
    972         971         941
    972         970         971
    100         413         989
    331         413          98
    129         989         413
    129         413         331
    100         101          99
    101          70          99
    100          99         413
    98         413          99
    443         445         414
    443         334         445
    476         445         334
    125         445         476
    125         446         445
    414         445         446
    575         577         576
    360         610         577
    579         577         610
    579         576         577
    360         577         315
    315         577         575
    948         914          37
    883         914         884
    913         914         883
    913          37         914
    746         139         951
    951         139         918
    922         139         746
    922         888         139
    634         970         957
    634         957         969
    937         969         957
    551         517         549
    551         549         580
    548         549         515
    515         549         517
    548         623         549
    580         549         623
    62         999         123
    999          91         123
    62         123          92
    122         123          91
    121         122          90
    122          91          90
    121          90         120
    120          90          58
    60          58          90
    60          90          91
    184         725         694
    754         725         722
    722         725         184
    815         694         725
    754         840         725
    815         725         840
    705         921         703
    921         920         703
    705         703         734
    702         734         703
    672         703         920
    672         702         703
    739         710         711
    21         711         710
    739         711         678
    679         678         711
    679         711         680
    21         680         711
    21         710         742
    774         775         742
    773         774         742
    773         742         710
    650         680         651
    927         651         652
    650         651         590
    590         651         927
    21         681         680
    652         681         714
    652         651         681
    680         681         651
    841         683         777
    32         683         841
    32         671         683
    671         448         683
    653         406         448
    715         777         406
    448         406         683
    777         683         406
    653         682         406
    682         715         406
    891         892         861
    892         862         861
    891         861         889
    859         889         861
    830         859         861
    830         861         862
    864         944         831
    860         830         831
    830         864         831
    860         831         829
    798         829         831
    798         831         944
    995         994         996
    756         992         960
    958         991         960
    991         756         960
    958         960         961
    858         825         606
    792         789         825
    792         825         826
    858         826         825
    882         918         915
    884         914         915
    948         882         915
    948         915         914
    888         919         139
    918         139         919
    887         919         888
    887         884         919
    918         919         915
    884         915         919
    25          24          23
    22          23          24
    22          24         758
    758          24         992
    58          69          87
    29         996          69
    29          69          58
    668         442         438
    752         381         438
    668         438         381
    412         438         442
    752         438         385
    412         385         438
    95          66         791
    66          38         791
    827          95         791
    827         791           2
    907         724         723
    937         957         874
    40         843         874
    843         842         874
    937         874         842
    40         938         227
    40         874         938
    957         938         874
    970         938         957
    972         938         970
    972         227         938
    92         123         154
    122         154         123
    155          92         154
    185         155         154
    185         154         153
    122         153         154
    714         681         743
    744         714         743
    744         743         775
    21         743         681
    21         742         743
    775         743         742
    25          26          24
    992          24          26
    962         961         168
    995         168         994
    962         168         963
    995         963         168
    606         825         824
    934         822         824
    606         824         822
    789         824         825
    934         824         790
    790         824         789
    25          27          26
    994          26          27
    996         994          27
    996          27          69
    87          69         195
    86         195         794
    86          87         195
    69          27         195
    25         794         195
    25         195          27
    56         971           3
    38           3         791
    2           3         971
    2         791           3
    38           4           3
    56           3           4
    70          68          99
    98          99          68
    70         723          68
    723         724          68
    994         993          26
    961         993         168
    994         168         993
    992          26         993
    961         960         993
    992         993         960
    159          67          38
    38          67           4
    159          98          67
    98          68          67
    724           4          67
    724          67          68
    312          56           5
    56           4           5
    312           5         973
    724           5           4
    907         973           5
    907           5         724];

UDSR.N = 1000;
UDSR.x = par(:,1);
UDSR.y = par(:,2);
UDSR.z = par(:,3);
UDSR.phi = par(:,4);
UDSR.theta = par(:,5);



end


