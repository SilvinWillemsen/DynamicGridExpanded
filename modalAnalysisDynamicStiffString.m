%{
    Dynamic stiff string...
%}

clear all;
% close all;
clc;

%% Draw settings
fs = 44100;             % Sample rate (in Hz)
k = 1/fs;               % Time step (in s)
lengthSound = 1000;   % Length of the simulation (in samples)

%{
    Where to add points along the system
    -1: Adding to the center alternating between left and right system.
    0: Undefined
    1: Right system has a single moving point.
    >2: (Expected behaviour) Selects where to add points (to left system).
%}

numFromBound = -1;

if numFromBound == 0
    error("numFromBound = 0 is undefined")
end

dispCorr = true;

%% Initialise parameters speeds for the entire simulation
Lvec = linspace(1, 2, lengthSound);   
rho = 7850;
r = 5e-4;
A = pi * r^2; 
I = pi / 4 * r^4;

Evec = linspace(sqrt(2e11), sqrt(2e11), lengthSound).^2;
E = Evec(1);
kappaSq = E * I / (rho * A);
kappaSqVec = Evec * I / (rho * A);

sig0 = 1;
sig1 = 0.005;

cVec = linspace(2940, 2250, lengthSound);
c = cVec(1);
h = sqrt((c^2 * k^2 + 4 * sig1 * k + sqrt((c^2 * k^2 + 4 * sig1 * k)^2 + 16 * kappaSq * k^2))/2);

Nfrac = Lvec(1) / h;
N = floor(Nfrac);
NPrev = N;

% calculate maximum number of points
hVec = sqrt((cVec.^2 * k^2 + 4 * sig1 * k + sqrt((cVec.^2 * k^2 + 4 * sig1 * k).^2 + 16 * kappaSqVec * k^2))/2);
Nmax = ceil(max(Lvec ./ hVec));


modesSave = zeros (lengthSound, Nmax);
sigmaSave = zeros (lengthSound, Nmax);
centDeviation = zeros (lengthSound, Nmax);

%% Initialise M and Mw (number of intervals between grid points for left and right system respectively)
if numFromBound == -1
    M = ceil (0.5 * N);
    Mw = floor (0.5 * N);
else
    M = N-numFromBound;
    Mw = numFromBound;   
end

%% Initialise states (excluding outer boundaries due to Dirichlet conditions)
uNext = zeros (M, 1);
u = zeros (M, 1);

qNext = zeros (N, 1);
q = zeros (N, 1);

u(1) = 1;
uPrev = u;

q(1:M) = u;
% q(startLoc:endLoc) = initState;
qPrev = q;

wNext = zeros (Mw, 1);
w = zeros (Mw, 1);
wPrev = zeros (Mw, 1);

%% Initialise Dxx matrices (matrix form of the \delta_{xx} operator)
eu = ones(length(u), 1);
Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));

ew = ones(length(w), 1);
Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));

%% Initialise output
out = zeros(floor(lengthSound), 1);

%% main loop
for n = 1:lengthSound  

    c = cVec(n);
    E = Evec(n);
    L = Lvec(n);
    kappaSq = E * I / (rho * A);

    % recalculate gridspacing and (fractional) number of intervals
    h = sqrt((c^2 * k^2 + 4 * sig1 * k + sqrt((c^2 * k^2 + 4 * sig1 * k)^2 + 16 * kappaSq * k^2))/2);
    Nfrac = L / h;
    N = floor (Nfrac);
    
    lambdaSq = c^2 * k^2 / h^2; 
    muSq = kappaSq * k^2 / h^4; % should always be 1

    % Calculate alpha (fractional part of Nfrac)
    alf = Nfrac - N;
    % Can only add/remove one point at a time
%     if abs(N - NPrev) > 1
%         error('Can only add/remove one grid point at a time')
%     end
    
    %% Check whether to add or remove points
    if N ~= NPrev
        
        %% Refresh M and Mw
        if numFromBound == -1
            M = ceil (0.5 * N);
            Mw = floor (0.5 * N);
        else
            M = N-numFromBound;
            Mw = numFromBound;
        end
        %% Refresh matrices
        eu = ones(M, 1);
        Dxxu = spdiags([eu -2*eu eu], -1:1, M,M);
        ew = ones(Mw, 1);
        Dxxw = spdiags([ew -2*ew ew], -1:1, Mw, Mw);
    end
    
    Dxx = zeros(M + Mw);
    Dxx(1:M, 1:M) = Dxxu;
    Dxx(M+1:end, M+1:end) = Dxxw;

    %% Calculate (quadratic) interpolator and virtual grid points
    ip = [-(alf - 1) / (alf + 1), 1, (alf - 1) / (alf + 1)];

    if numFromBound == 1
        Dxx(M, M:(M+1)) = Dxx(M, M:(M+1)) + fliplr(ip(2:end));
        Dxx(M+1, (M-1):(M+1)) = Dxx(M+1, (M-1):(M+1)) + ip;
    else
        Dxx(M, M:(M+2)) = Dxx(M, M:(M+2)) + fliplr(ip);
        Dxx(M+1, (M-1):(M+1)) = Dxx(M+1, (M-1):(M+1)) + ip;
    end
    
    Dxx = Dxx / h^2;
    Dxxxx = Dxx * Dxx;
    Id  = speye(N);         % identity matrix

    sigX = 1;
    b = (1-alf)/(alf + eps);
    Ju = [zeros(M-1, 1); 1/h; -1/h; zeros(Mw-1, 1)];
    etaVec = [zeros(1, M-1), -1, 1, zeros(1, Mw-1)];

    Amat = (1 + sig0 * k) * Id;% - b * k^2 * (1 + sigX/k) / 2 * Ju * etaVec;
    B = 2 * Id + c^2 * k^2 * Dxx - kappaSq * k^2 * Dxxxx + 2 * sig1 * k * Dxx;
    C = -(1 - sig0 * k) * Id - 2 * sig1 * k * Dxx;% + b * k^2 * (1 - sigX/k) / 2 * Ju * etaVec;

    [f, sigma, ~] = analyseQ([Amat\B, Amat\C; ...
              eye(size(B)), zeros(size(B))], k);

    modesSave(n, 1:(M+Mw)) = f(1:(M+Mw));
    sigmaSave(n, 1:(M+Mw)) = sigma(1:(M+Mw));
    %     subplot(211)
%     hold off;
%     plot(f(1:(M+Mw)))
%     hold on;
%     beta = (1:length(f))*pi; % bar 
    p = 1:(M+Mw);
    if kappaSq ~= 0
        % gives machine imprecision for lowewr values of kappaSq
        beta = sqrt(((sqrt(4 * pi^4 * kappaSq^2 * p.^4 / L^4 + c^4 + 4 * c^2 * pi^2 * kappaSq * p.^2 / L^2)) - c^2)/(2*kappaSq));
    else
        beta = p * pi / L;
    end
    expectedF = 1/(pi*k) * asin(sqrt(lambdaSq * sin(beta * h / 2).^2 + 4 * muSq * sin(beta*h / 2).^4));
%     plot(expectedF);
    
    subplot(212)
    centDeviation(n, 1:(M+Mw)) = 1200 * log2(f(1:(M+Mw))./expectedF');
%     plot (centDeviation(n, 1:(M+Mw)));
%     drawnow;
   

    NPrev = N;
    
end

figure;
modesSave(modesSave == 0) = nan;
plot(modesSave(1:n, :))
