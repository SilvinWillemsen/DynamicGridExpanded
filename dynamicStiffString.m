%{
    Dynamic stiff string...
%}

clear all;
% close all;
clc;

%% Draw settings
drawThings = false;
drawStart = 2000;
drawSpeed = 1;
plotModalAnalysis = false;

fs = 44100;             % Sample rate (in Hz)
k = 1/fs;               % Time step (in s)
if plotModalAnalysis
    lengthSound = 1000;   % Length of the simulation (in samples)
else
    lengthSound = fs * 5;
end

%{
    Where to add points along the system
    -1: Adding to the center alternating between left and right system.
    0: Undefined
    1: Right system has a single moving point.
    >2: (Expected behaviour) Selects where to add points (to left system).
%}

numFromBound = 1;

if numFromBound == 0
    error("numFromBound = 0 is undefined")
end

dispCorr = true;
%% Initialise wave speeds for the entire simulation
L = 1;                  % System length (does not really matter)
rho = 7850;
r = 5e-4;
A = pi * r^2; 
Evec = linspace(sqrt(2e11), sqrt(2e14), lengthSound*1/4).^2;
Evec = [ones(1, floor(lengthSound * 1/8)) *  Evec(1), Evec];
Evec = [Evec, ones(1, floor(lengthSound * 1/8)) * Evec(end)];
Evec = [Evec, fliplr(Evec)];
Evec = [Evec, ones(1, lengthSound - length(Evec)) * Evec(end)];

% Evec = linspace(sqrt(0), sqrt(0), lengthSound).^2;

% Evec = [Evec, ones(1, lengthSound - length(Evec)) * Evec(end)];
E = Evec(1);
I = pi / 4 * r^4;
kappaSq = E * I / (rho * A);
kappaSqVec = Evec * I / (rho * A);

sig0 = 1;
sig1 = 0.005;

halfWidth = 5;
initState = hann(2*halfWidth + 1);
cVec = linspace(300, 300, lengthSound);
c = cVec(1);
h = sqrt((c^2 * k^2 + 4 * sig1 * k + sqrt((c^2 * k^2 + 4 * sig1 * k)^2 + 16 * kappaSq * k^2))/2);
Nfrac = L / h;
N = floor(Nfrac);
NPrev = N;

hVec = sqrt((cVec.^2 * k^2 + 4 * sig1 * k + sqrt((cVec.^2 * k^2 + 4 * sig1 * k).^2 + 16 * kappaSqVec * k^2))/2);
Nmax = ceil(max(L ./ hVec));

modesSave = zeros (lengthSound, Nmax);
sigmaSave = zeros (lengthSound, Nmax);

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

    if mod(n, fs/2) == 1
        q(1) = q(1) + 1/h;
    end
    % Retrieve new c
%     T = T0 + u' * u;
%     c = sqrt(T / (rho * A));
    c = cVec(n);
    E = Evec(n);
    kappaSq = E * I / (rho * A);

    % recalculate gridspacing and (fractional) number of intervals
    h = sqrt((c^2 * k^2 + 4 * sig1 * k + sqrt((c^2 * k^2 + 4 * sig1 * k)^2 + 16 * kappaSq * k^2))/2);
    Nfrac = L / h;
    N = floor (Nfrac);
    
    lambdaSq = c^2 * k^2 / h^2; % should always be 1

    % Calculate alpha (fractional part of Nfrac)
    alf = Nfrac - N;
    % Can only add/remove one point at a time
    if abs(N - NPrev) > 1
        error('Can only add/remove one grid point at a time')
    end
    
    %% Check whether to add or remove points
    if N ~= NPrev
        
        % Note that the "next" vectors (at n+1) are appended/prepended with a 0 as these will be overwritten by calculation of scheme anyway
        if N > NPrev
        %% Add point if N^n > N^{n-1}
        
            % Create cubic interpolator 
            cubicIp = [alf * (alf + 1) / -((alf + 2) * (alf + 3)); ...
                        2 * alf / (alf + 2); ...
                        2 / (alf + 2); ...
                        2 * alf / -((alf + 3) * (alf + 2))]';
             
            qNextLeft = qNext(1:M);
            qNextRight = qNext(M+1:end);
            
            qLeft = q(1:M);
            qRight = q(M+1:end);
            
            qPrevLeft = qPrev(1:M);
            qPrevRight = qPrev(M+1:end);
            
            % Add grid points to u and w in an alternating fashion 
            if numFromBound == -1
                if mod(N, 2) == 1   
                    % if N is odd, add to left system
                    uNext = [uNext; 0]; % append uNext without interpolation (will be overwritten by calculation of scheme anyway)
                    u = [u; cubicIp * [u(end-1:end); w(1:2)]];
                    uPrev = [uPrev; cubicIp * [uPrev(end-1:end); wPrev(1:2)]];
                    

                else
                    % if N is even, add to right system
                    wNext = [0; wNext]; % prepend wNext without interpolation (will be overwritten by calculation of scheme anyway)
                    w = [fliplr(cubicIp) * [u(end-1:end); w(1:2)]; w];
                    wPrev = [fliplr(cubicIp) * [uPrev(end-1:end); wPrev(1:2)]; wPrev];
                    
%                     qNext = [qNextLeft; fliplr(cubicIp) * [qNextLeft(M-1:M); qNextRight(1:2)]; qNextRight];
%                     q = [qLeft; fliplr(cubicIp) * [qLeft(M-1:M); qRight(1:2)]; qRight];
%                     qPrev = [qPrevLeft; fliplr(cubicIp) * [qPrevLeft(M-1:M); qPrevRight(1:2)]; qPrevRight];

                end
                qNext = [qNextLeft; cubicIp * [qNextLeft(M-1:M); qNextRight(1:2)]; qNextRight];
                q = [qLeft; cubicIp * [qLeft(M-1:M); qRight(1:2)]; qRight];
                qPrev = [qPrevLeft; cubicIp * [qPrevLeft(M-1:M); qPrevRight(1:2)]; qPrevRight];

            % Otherwise, add grid points to u using interpolation
            elseif numFromBound == 1
                uNext = [uNext; 0];
                u = [u; cubicIp(1:3) * [u(end-1:end); w(1)]];
                uPrev = [uPrev; cubicIp(1:3) * [uPrev(end-1:end); wPrev(1)]];
                qNext = [qNextLeft; cubicIp(1:3) * [qNextLeft(M-1:M); qNextRight(1)]; qNextRight];
                q = [qLeft; cubicIp(1:3) * [qLeft(M-1:M); qRight(1)]; qRight];
                qPrev = [qPrevLeft; cubicIp(1:3) * [qPrevLeft(M-1:M); qPrevRight(1)]; qPrevRight];

            else
                uNext = [uNext; 0];
                u = [u; cubicIp * [u(end-1:end); w(1:2)]];
                uPrev = [uPrev; cubicIp * [uPrev(end-1:end); wPrev(1:2)]];
                qNext = [qNextLeft; cubicIp * [qNextLeft(M-1:M); qNextRight(1:2)]; qNextRight];
                q = [qLeft; cubicIp * [qLeft(M-1:M); qRight(1:2)]; qRight];
                qPrev = [qPrevLeft; cubicIp * [qPrevLeft(M-1:M); qPrevRight(1:2)]; qPrevRight];

            end
        
        else
        %% Remove point if N^n < N^{n-1}
            % Remove grid points from u and w in an alternating fashion 
            if numFromBound == -1 
                if mod(N,2) == 0
                    % if N is even, remove from left system
                    uNext = uNext(1:end-1);
                    u = u(1:end-1);
                    uPrev = uPrev(1:end-1);
                    qNext(M) = [];
                    q(M) = [];
                    qPrev(M) = [];
                else 
                    % if N is odd, remove from right system
                    wNext = wNext(2:end);
                    w = w(2:end);
                    wPrev = wPrev(2:end);
                    qNext(M+1) = [];
                    q(M+1) = [];
                    qPrev(M+1) = [];

                end
                
            % Otherwise, remove grid point from u
            else
                uNext = uNext(1:end-1);
                u = u(1:end-1);
                uPrev = uPrev(1:end-1);
                qNext(M) = [];
                q(M) = [];
                qPrev(M) = [];

            end

            
        end
        
        %% Refresh matrices
        eu = ones(length(u), 1);
        Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));
        ew = ones(length(w), 1);
        Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));

        %% Refresh M and Mw
        if numFromBound == -1
            M = ceil (0.5 * N);
            Mw = floor (0.5 * N);
        else
            M = N-numFromBound;
            Mw = numFromBound;
        end
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

    Amat = (1 + sig0 * k);
    B = 2 * Id + c^2 * k^2 * Dxx - kappaSq * k^2 * Dxxxx + 2 * sig1 * k * Dxx;
    C = -(1 - sig0 * k) * Id - 2 * sig1 * k * Dxx;

    if plotModalAnalysis
        [f, sigma, ~] = analyseQ([Amat\B, Amat\C; ...
                  eye(size(B)), zeros(size(B))], k);

        modesSave(n, 1:length(f)) = f;
        sigmaSave(n, 1:length(f)) = sigma;
    end
    % virtualGridPoints = [u_{M+1}^n; w_{)-1}^n]
    if numFromBound == 1 % use boundary value for w(2) (= 0)          
        virtualGridPoints = [ip(3) * u(end) + ip(2) * w(1) + 0; ...
                              ip(1) * u(end-1) + ip(2) * u(end) + ip(3) * w(1)]; 
    else
        virtualGridPoints = [ip(3) * u(end) + ip(2) * w(1) + ip(1) * w(2); ...
                              ip(1) * u(end-1) + ip(2) * u(end) + ip(3) * w(1)];
    end
    
    %% Calculate left system
    uNext = (2 * u + lambdaSq * Dxxu * u - (1 - sig0 * k) * uPrev) / (1 + sig0 * k);
    uNext(end) = uNext(end) + lambdaSq * virtualGridPoints(1);
    
    %% Calculate right system
    wNext = (2 * w + lambdaSq * Dxxw * w  - (1 - sig0 * k) * wPrev) / (1 + sig0 * k);
    wNext(1) = wNext(1) + lambdaSq * virtualGridPoints(2);

    %% Calculate using matrix
    qNext = B / Amat * q + C / Amat * qPrev;
    
    %% Displacement correction
    if dispCorr
        epsilon = 0; % Calculation is still defined for epsilon = 0 
        sig0 = 1; % Damping coefficient

        %% Calculate correction force
        etaPrev = qPrev(M+1) - qPrev(M);
        rForce = (1 - sig0 / k) / (1 + sig0 / k);
        oOP = (h * (1 + sig0 / k) * (1-alf)) / (2 * h * (alf + epsilon) + 2 * k^2 * (1 + sig0 / k) * (1-alf));
        
        % Here, uNext and wNext are the 'intermediate' states of u and w (schemes without connection forces)
        F = (qNext(M+1) - qNext(M) + rForce * etaPrev) * oOP;
        
        %% Add displacement correction to inner boundaries
        qNext(M) = qNext(M) + k^2/h * F;
        qNext(M+1) = qNext(M+1) - k^2/h * F;

    end

    %% save output
    out(n) = q(5);

    %% draw stuff
    if n > drawStart && drawThings && mod(n, drawSpeed) == 0

        % Grid point locations
        hLocsLeft = (0:(length(u)))*h;
        hLocsRight = (fliplr(L - (0:(length(w)))*h));
        
%         hold off;
%         
%         % Plot left system (with left outer boundary)
%         plot(hLocsLeft, [0;u], 'LineWidth' ,2, 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
%         hold on;
%         
%         % Plot right system (with right outer boundary)
%         plot(hLocsRight, [w; 0], 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
        
        plot([hLocsLeft, hLocsRight], [0; q; 0], 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'g');
        % Settings
        xlim([0, L])
%         ylim([-1, 1])
        grid on;
        xlabel("$x$ (m)", 'interpreter', 'latex')
        ylabel("displacement (m)")
        title("$\alpha = " + num2str(round(alf*100)/100) + "$", 'interpreter', 'latex')
%         legend(["$u$", "$w$"], 'interpreter', 'latex')
        set(gca, "Fontsize", 16, 'Linewidth', 2)
        
        drawnow;
    end
    
    %% Update states
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;
    
    qPrev = q;
    q = qNext;

    NPrev = N;
    
end

figure;
if plotModalAnalysis
    modesSave(modesSave == 0) = nan;
    plot(modesSave(1:n, :))
else
    spectrogram (out,512,64,512, 44100, 'yaxis');
end