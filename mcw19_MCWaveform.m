function g = mcw19_MCWaveform(tvec, n, intervals, ramp)
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
            s = sprintf('g%d = mcw19_trapez(t,T*intervals(i),T*intervals(i+1),ramp*T);', i);
        else
            s = sprintf('g%d = mcw19_sin(t,T*intervals(i),T*intervals(i+1));', i);
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







