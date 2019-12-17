% Compute conditional spectral Granger causality for the multi-variate
% using a non parametric method. 
% 
% For more information, we refer to:
% Detto et al., 2012, Causality and Persistence in Ecological Systems: A
% Nonparametric Spectral Granger Causality Approach
% 
% and
% Dhamala et al., 2008, Estimating Granger Causality from Fourier and
% Wavelet Transforms of Time Series Data

function [CSGC, Fr] = CalculateCSGC(X, varargin)

% X is a cell containg a number of realisations for each time series
% varargin contains parameters and there value if they need to deviate from
% the default

%% Setting parameters
% Default parameters
opts = struct('mother', 'Morlet',...    % mother wavelet (need to be complex)
    'param', 15,...                     % wavelet base parameter
    'dt', 1,...                         % temporal step
    'max_iter', 30,...                   % max iterations for Wilson algorithm
    'tol', 1e-6,...                     % tolerance for Wilson algorithm
    'fextrap', 0,...                    % limit for spectra extrapolation
    'extrap_method', 'param',...        % method of extrapolation (param or interp)
    'order', 2,...                      % order of the AR model for arfit.m
    'disc', 100,...                     % normalise each trial
    'anom', 0);                         % remove periodicity anom
                       


% Switch parameters
opts = parseArgs(varargin, opts);

% Define parameters for internal use
m = opts.order;
dt = opts.dt;
[n, R] = size(X{1});
dth = pi/(opts.disc-1);                 % Interval of discretisation of angular frequency
M = length(X);
max_iter = opts.max_iter;
tol = opts.tol;
theta = 0:dth:pi;                      % Discretisation of fourier space (in angular frequency)
T = length(theta);
S = zeros(M, M, T);
Fr = theta/2/pi/dt;

% Parameters and constant of wavelet
dj = 0.1;
Cdelta = 0.776;
J = -1;
param = opts.param;
fourier_factor = (4*pi)/(param + sqrt(2 + param^2));
df = Cdelta*2*pi*log(2)/dt;             % spectral interval
sc = 1./Fr/fourier_factor;

%% Initialisation
f = zeros(M, M, T);
A = nan(M, m*M, R);
C = nan(M, M, R);

%% Wavelet transformation
for r = 1:R
    for j = 1:M
        v(:, j) = X{j}(:, r);        
    end
    
    WT = wave(v, 'scale', sc(2:end), 'dj', dj, 'dt', dt, 'param', param);
%     WT = wave(v, 'scale', sc(2:end));
    f(:, :, 2:end) = f(:, :, 2:end) + WT;
    if strcmp(opts.extrap_method, 'param')
        % fit parametric model for extrapolation
        [~, A(:, :, r), C(:, :, r)] = arfit(v, m, m);
    end
end

%% Extrapolate spectra at f = 0
if strcmp(opts.extrap_method, 'param')
    A(abs(A) == inf) = nan;
    S = AR_spectrum(nanmean(A, 3), nanmean(C, 3), Fr*dt);
    S = S*dt;
    use=find(Fr>=opts.fextrap & Fr>Fr(1));
    S(:, :, use) = f(:, :, use)/R/df;
else
    S = spectral_matrix_extrap(f/R/df,Fr,opts.fextrap);
end
%% Determine conditional spectral Granger causality
[SIGMA4, H]= wilson(S, max_iter, tol);
b = 0;

for i = 1:M
    
    for j = 1:M
        b = b + 1;
        if i~=j
            sub = find((1:M~=i) & (1:M~=j));
            sub1 = [i, sub];
            sub2 = [i, j, sub];

            [SIGMA3, G] = wilson(S(sub1, sub1, :), max_iter, tol);
            CSGC{i, j} = Granger_multicond(G, H(sub2, sub2, :), SIGMA3, SIGMA4(sub2, sub2));
        end
    end
end



    function  S = spectral_matrix_extrap(S,fr,fextrap)
        if fextrap==0
            fextrap=fr(2);
        end
        gg=find(fr<=fextrap);                   %extrapolate after fextrap
        hh=find(fr<fextrap/2 & fr>fextrap);     %retain for interpolation
        
        if isempty(hh)==1
            hh=2;
        end
        a=0;
        method='cubic';
        
        for i=1:M
            for j=1:M
                
                if i==j
                    yy(1,:)=S(i,j,hh);
                    y1=fliplr(yy);
                    y2=yy;
                    x1=-fliplr(fr(hh));
                    x2=fr(hh);
                    xi=fr(gg);
                    
                    S(i,j,gg)=interp1([x1 x2],[y1 y2],xi,method);
                    
                else
                    yy1(1,:)=S(j,i,hh);
                    yy2(1,:)=S(i,j,hh);
                    y1=fliplr(yy1);
                    y2=yy2;
                    x1=-fliplr(fr(hh));
                    x2=fr(hh);
                    xi=fr(gg);
                    
                    S(i,j,gg)=interp1([x1 x2],real([y1 y2]),xi,method)+...
                        sqrt(-1)*interp1([x1 0 x2],imag([y1 0 y2]),xi,method);
                    
                end
            end 
        end 
    end

end
            






