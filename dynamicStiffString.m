%{
    Dynamic stiff string...
%}

clear all;
close all;
clc;

%% Draw settings
drawThings = false;
drawStart = 145530;
drawSpeed = 5;
plotModalAnalysis = false;

fs = 44100;             % Sample rate (in Hz)
k = 1/fs;               % Time step (in s)
if plotModalAnalysis
    lengthSound = 1000;   % Length of the simulation (in samples)
else
    lengthSound = fs * 3;
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
% L, rho, r, T, E, sig0, sig1
% params = [0.3, 7850, 5e-4, 00, 2e14, 1, 0.01]';
f0 = 200;
7850 * 0.0002794^2 * pi
Tuse = (f0^2 * 78.50 * pi * (0.002794)^2) * (2 * 1)^2;
% params = [1, 7850, 5e-4, 555, 2e11, 1, 0.005;
%           1, 7850, 5e-4, 0, 7e12, 1, 0.005;
%           0.3, 7850, 5e-4, 0, 7e12, 1, 0.005;
%           0.3, 7850, 5e-4, 0, 2e12, 1, 0.01;]';
rUse = 0.0002794;
E1 = 2e11;
E2 = 2e13;
T1 = (f0^2 * 7850 * pi * rUse^2) * (2 * 1)^2;
T2 = T1 + (E1 - E2) * rUse^4 * pi^3 / (4 * 1^2);
params = [1, 7850, rUse, T1, E1, 1, 0.005;
          1, 7850, rUse, T2, E2, 1, 0.005;
          ]';
%           1, 7850, 5e-4, 0, 2e11, 1, 0.01;]';
%           0.5, 7850, 5e-4, 300, sqrt(2e11), 1, 0.005;
%           0.3, 7850, 5e-4, 300, sqrt(2e11), 1, 0.005;]';
     
changeRatio = 0.5;

numParams = size(params, 1);
numChanges = size(params, 2);
chunkSize = floor(lengthSound / (numChanges*2 - 1));
paramVecs = zeros(numParams, lengthSound);

paramVecs = params(:, 1) * ones(1, chunkSize);
for i = 1:numChanges-1
    paramMat = [];
    for j = 1:numParams
        if j == 1
            paramMat = [paramMat; sqrt(1./linspace(1/(params(j, i))^2, 1/(params(j, i+1))^2, chunkSize))];
        elseif j == 2
            paramMat = [paramMat; 1./linspace(1/sqrt(params(j, i)), 1/sqrt(params(j, i+1)), chunkSize).^2];
        elseif j == 4 || j == 5
            paramMat = [paramMat; linspace(sqrt(params(j, i)), sqrt(params(j, i+1)), chunkSize).^2];
%         elseif j == 3 
%             paramMat = [paramMat; linspace(sqrt(sqrt(params(j, i))), sqrt(sqrt(params(j, i+1))), chunkSize).^4];
        else
            paramMat = [paramMat; linspace(params(j, i), params(j, i+1), chunkSize)];
        end
    end
    paramVecs = [paramVecs, paramMat];

    paramVecs = [paramVecs, params(:, i+1) .* ones(numParams, chunkSize)];

end
Lvec = paramVecs(1, :);
rhoVec = paramVecs(2, :);
rVec = paramVecs(3, :);
Tvec = paramVecs(4, :);
Evec = paramVecs(5, :);
sig0Vec = paramVecs(6, :);
sig1Vec = paramVecs(7, :);
% Lvec = [startParams(1) * ones(1, floor(0.5 * (1-changeRatio) * lengthSound)), linspace(startParams(1), endParams(1), floor(changeRatio * lengthSound))];
% rhoVec = [startParams(2) * ones(1, floor(0.5 * (1-changeRatio) * lengthSound)), linspace(startParams(2), endParams(2), floor(changeRatio * lengthSound))];
% rVec = [startParams(3) * ones(1, floor(0.5 * (1-changeRatio) * lengthSound)), linspace(startParams(3), endParams(3), floor(changeRatio * lengthSound))];
% Tvec = [startParams(4) * ones(1, floor(0.5 * (1-changeRatio) * lengthSound)), linspace(startParams(4), endParams(4), floor(changeRatio * lengthSound))];
% Evec = [startParams(5) * ones(1, floor(0.5 * (1-changeRatio) * lengthSound)), linspace(startParams(5), endParams(5), floor(changeRatio * lengthSound))];
% sig0Vec = [startParams(6) * ones(1, floor(0.5 * (1-changeRatio) * lengthSound)), linspace(startParams(6), endParams(6), floor(changeRatio * lengthSound))];
% sig1Vec = [startParams(7) * ones(1, floor(0.5 * (1-changeRatio) * lengthSound)), linspace(startParams(7), endParams(7), floor(changeRatio * lengthSound))];
% Lvec = [Lvec, ones(1, lengthSound - length(Lvec)) * endParams(1)];
% rhoVec = ([rhoVec, ones(1, lengthSound - length(rhoVec)) * endParams(2)]).^2;
% rVec = ([rVec, ones(1, lengthSound - length(rVec)) * endParams(3)]).^2;
% Tvec = sqrt([Tvec, ones(1, lengthSound - length(Tvec)) * endParams(4)]).^2;
% Evec = ([Evec, ones(1, lengthSound - length(Evec)) * endParams(5)]).^2;
% sig0Vec = [sig0Vec, ones(1, lengthSound - length(sig0Vec)) * endParams(6)];
% sig1Vec = [sig1Vec, ones(1, lengthSound - length(sig1Vec)) * endParams(7)];

% Evec = linspace(sqrt(2e11), sqrt(2e12), lengthSound*1/4).^2;
% Evec = [ones(1, floor(lengthSound * 1/8)) *  Evec(1), Evec];
% Evec = [Evec, ones(1, floor(lengthSound * 1/8)) * Evec(end)];
% Evec = [Evec, fliplr(Evec)];
% Evec = [Evec, ones(1, lengthSound - length(Evec)) * Evec(end)];

% Evec = linspace(sqrt(0), sqrt(0), lengthSound).^2;

% Evec = [Evec, ones(1, lengthSound - length(Evec)) * Evec(end)];

A = pi * rVec(1)^2;
I = pi / 4 * rVec(1)^4;
cSq = Tvec(1) / (rhoVec(1) * A);
kappaSq = Evec(1) * I / (rhoVec(1) * A);

h = sqrt((cSq * k^2 + 4 * sig1Vec(1) * k + sqrt((cSq * k^2 + 4 * sig1Vec(1) * k)^2 + 16 * kappaSq * k^2))/2);
Nfrac = Lvec(1) / h;
N = floor(Nfrac);
NPrev = N;

cSqVec = Tvec ./ (rhoVec * A);
kappaSqVec = (Evec .* pi / 4 .* rVec.^4) ./ (rhoVec * A);

hVec = sqrt((cSqVec * k^2 + 4 * sig1Vec * k + sqrt((cSqVec * k^2 + 4 * sig1Vec * k).^2 + 16 * kappaSqVec * k^2))/2);
Nmax = ceil(max(Lvec ./ hVec));

modesSave = zeros (lengthSound, Nmax);
sigmaSave = zeros (lengthSound, Nmax);

pIdx = lengthSound-1;
p = 1;
if cSqVec(pIdx) == 0
    f0 = sqrt(kappaSqVec(pIdx)) * pi / 2 * (p)^2 / L^2;
else
    f0 = sqrt(cSqVec(pIdx)) * p / (2 * Lvec(pIdx)) * sqrt(1 + kappaSqVec(pIdx) * pi^2 / (cSqVec(pIdx) * Lvec(pIdx)^2) * p^2);
end

%% Initialise M and Mw (number of intervals between grid points for left and right system respectively)
if numFromBound == -1
    M = ceil (0.5 * N);
    Mw = floor (0.5 * N);
else
    M = N-numFromBound;
    Mw = numFromBound;   
end

%% Initialise states (excluding outer boundaries due to Dirichlet conditions)
qNext = zeros (N, 1);
q = zeros (N, 1);
qPrev = q;
qMat = q;
qPrevMat = q;

%% Matrices
DxxFull = zeros(M + Mw);

eu = ones(M, 1);
ew = ones(Mw, 1);
DxxFull(1:M, 1:M) = spdiags([eu -2*eu eu], -1:1, M, M);
DxxFull(M+1:end, M+1:end) = spdiags([ew -2*ew ew], -1:1, Mw, Mw);
DxxFull = DxxFull;

DxxxxFull = DxxFull * DxxFull;

        
%% Initialise output
out = zeros(floor(lengthSound), 1);

%% main loop
for n = 1:lengthSound  

    if mod(n, chunkSize / 2) == 1
%         q = q + 1/h * hann(N);
%         qMat = qMat + 1/h * hann(N);
        div = 10;
        q(1:floor(N/div)) = q(1:floor(N/div)) + 1/h * hann(floor(N/div));
        qMat(1:floor(N/div)) = qMat(1:floor(N/div)) + 1/h * hann(floor(N/div));
    end
    
    % Retrieve new params
    A = pi * rVec(n)^2;
    I = pi / 4 * rVec(n)^4;
    cSq = Tvec(n) / (rhoVec(n) * A);
    kappaSq = Evec(n) * I / (rhoVec(n) * A);

    h = sqrt((cSq * k^2 + 4 * sig1Vec(n) * k + sqrt((cSq * k^2 + 4 * sig1Vec(n) * k)^2 + 16 * kappaSq * k^2))/2);
    Nfrac = Lvec(n) / h;
    N = floor(Nfrac);

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
                if mod(N, 2) == 0
                    qNext = [qNextLeft; fliplr(cubicIp) * [qNextLeft(M-1:M); qNextRight(1:2)]; qNextRight];
                    q = [qLeft; fliplr(cubicIp) * [qLeft(M-1:M); qRight(1:2)]; qRight];
                    qPrev = [qPrevLeft; fliplr(cubicIp) * [qPrevLeft(M-1:M); qPrevRight(1:2)]; qPrevRight];
                else
                    qNext = [qNextLeft; cubicIp * [qNextLeft(M-1:M); qNextRight(1:2)]; qNextRight];
                    q = [qLeft; cubicIp * [qLeft(M-1:M); qRight(1:2)]; qRight];
                    qPrev = [qPrevLeft; cubicIp * [qPrevLeft(M-1:M); qPrevRight(1:2)]; qPrevRight];

                end
            elseif numFromBound == 1
                qNext = [qNextLeft; cubicIp(1:3) * [qNextLeft(M-1:M); qNextRight(1)]; qNextRight];
                q = [qLeft; cubicIp(1:3) * [qLeft(M-1:M); qRight(1)]; qRight];
                qPrev = [qPrevLeft; cubicIp(1:3) * [qPrevLeft(M-1:M); qPrevRight(1)]; qPrevRight];

            else
                qNext = [qNextLeft; cubicIp * [qNextLeft(M-1:M); qNextRight(1:2)]; qNextRight];
                q = [qLeft; cubicIp * [qLeft(M-1:M); qRight(1:2)]; qRight];
                qPrev = [qPrevLeft; fliplr(cubicIp) * [qPrevLeft(M-1:M); qPrevRight(1:2)]; qPrevRight];
            end
            disp("Added Point")

        else
        %% Remove point if N^n < N^{n-1}
            diffSave(n) = (q(M) - q(M+1));
            % Remove grid points from u and w in an alternating fashion 
            if numFromBound == -1 
                if mod(N,2) == 0
                    % if N is even, remove from left system
                    qNext(M) = [];
                    q(M) = [];
                    qPrev(M) = [];
                else 
                    % if N is odd, remove from right system
                    qNext(M+1) = [];
                    q(M+1) = [];
                    qPrev(M+1) = [];

                end
                
            % Otherwise, remove grid point from u
            else
                qNext(M) = [];
                q(M) = [];
                qPrev(M) = [];

            end
            disp("Removed Point")
            
        end
        
        %% Refresh M and Mw
        if numFromBound == -1
            M = ceil (0.5 * N);
            Mw = floor (0.5 * N);
        else
            M = N-numFromBound;
            Mw = numFromBound;
        end
        
        DxxFull = zeros(M + Mw);
    
        eu = ones(M, 1);
        ew = ones(Mw, 1);
        DxxFull(1:M, 1:M) = spdiags([eu -2*eu eu], -1:1, M, M);
        DxxFull(M+1:end, M+1:end) = spdiags([ew -2*ew ew], -1:1, Mw, Mw);
        
        DxxxxFull = DxxFull * DxxFull;
        
    end
    
    NPrev = N;
%     alf = 0;
    %% Calculate (quadratic) interpolator and virtual grid points
    ip = [-(alf - 1) / (alf + 1), 1, (alf - 1) / (alf + 1)];
    DxxxxAlfMatrix = [-4, 6, ip(3) - 4, 1, ip(1), 0;
                    1, -4, (ip(3) - 2)^2 + 2, ip(3) - 4, 4 * (alf^2 + alf - 1) / (alf + 1)^2, ip(1)];
    DxxxxAlfMatrix = [DxxxxAlfMatrix; rot90(DxxxxAlfMatrix, 2)];
    Dxx = DxxFull / h^2;
    Dxxxx = DxxxxFull / h^4;
    if numFromBound == 1
        Dxx(M, M:(M+1)) = Dxx(M, M:(M+1)) + fliplr(ip(2:end)) / h^2;
        Dxx(M+1, (M-1):(M+1)) = Dxx(M+1, (M-1):(M+1)) + ip / h^2;
%         Dxxxx(M-1:M+2, M-2:M+3) = DxxxxAlfMatrix / h^4;    
        Dxxxx = Dxx * Dxx;
    else
        Dxx(M, M:(M+2)) = Dxx(M, M:(M+2)) + fliplr(ip) / h^2;
        Dxx(M+1, (M-1):(M+1)) = Dxx(M+1, (M-1):(M+1)) + ip / h^2;
        Dxxxx(M-1:M+2, M-2:M+3) = DxxxxAlfMatrix / h^4;    
%         Dxxxx = Dxx * Dxx;
    end

    Id  = speye(N);         % identity matrix

    Amat = (1 + sig0Vec(n) * k);
    B = 2 * Id + cSq * k^2 * Dxx - kappaSq * k^2 * Dxxxx + 2 * sig1Vec(n) * k * Dxx;
    C = -(1 - sig0Vec(n) * k) * Id - 2 * sig1Vec(n) * k * Dxx;

    if plotModalAnalysis
        [f, sigma, ~] = analyseQ([Amat\B, Amat\C; ...
                  eye(size(B)), zeros(size(B))], k);

        modesSave(n, 1:length(f)) = f;
        sigmaSave(n, 1:length(f)) = sigma;
    end

    %% Calculate using matrix
    qNext = B / Amat * q + C / Amat * qPrev;
    
    %% Displacement correction
    if dispCorr
        epsilon = 0; % Calculation is still defined for epsilon = 0 
        sigX = 1; % Damping coefficient

        %% Calculate correction force
        etaPrev = qPrev(M+1) - qPrev(M);
        rForce = (1 - sigX / k) / (1 + sigX / k);
        oOP = (h * (1 + sigX / k) * (1-alf)) / (2 * h * (alf + epsilon) + 2 * k^2 * (1 + sigX / k) * (1-alf));
        
        % Here, uNext and wNext are the 'intermediate' states of u and w (schemes without connection forces)
        F = (qNext(M+1) - qNext(M) + rForce * etaPrev) * oOP;
        
        %% Add displacement correction to inner boundaries
        qNext(M) = qNext(M) + k^2/h * F;
        qNext(M+1) = qNext(M+1) - k^2/h * F;

    end

%     beta = (1 - alf) / (alf);
%     jVec = zeros(M + Mw, 1);
%     jVec(M) = 1/h;
%     jVec(M+1) = -1/h;
%     etaVec = zeros(1, M + Mw);
%     etaVec(M) = -1;
%     etaVec(M+1) = 1;
%     
%     Adisp = (1 + sig0Vec(n) * k) * Id - beta * k^2 * (1 + sigX / k) / 2 * jVec * etaVec;
%     Cdisp = -(1 - sig0Vec(n) * k) * Id - 2 * sig1Vec(n) * k * Dxx + beta * k^2 * (1 - sigX / k) / 2 * jVec * etaVec;
%     qNextMat = Adisp \ (B * qMat + Cdisp * qPrevMat);
%     hold off;
%     plot(qNext);
%     hold on;
%     plot(qNextMat);
    %% save output
    out(n) = q(5);

    %% draw stuff
    if n > drawStart && drawThings && mod(n, drawSpeed) == 0

        % Grid point locations
        hLocsLeft = (0:M)*h;
        hLocsRight = (fliplr(Lvec(n) - (0:Mw)*h));
        
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
        xlim([0, Lvec(n)])
        ylim([-1000, 1000])
        grid on;
        xlabel("$x$ (m)", 'interpreter', 'latex')
        ylabel("displacement (m)")
        title("$\alpha = " + num2str(round(alf*100)/100) + "$", 'interpreter', 'latex')
%         legend(["$u$", "$w$"], 'interpreter', 'latex')
        set(gca, "Fontsize", 16, 'Linewidth', 2)
        
        drawnow;
    end
    
    %% Update states    
    qPrev = q;
    q = qNext;
%     qPrevMat = qMat;
%     qMat = qNextMat;
    NPrev = N;
    
end

figure;
if plotModalAnalysis
    modesSave(modesSave == 0) = nan;
    plot(modesSave(1:n, :))
else
    spectrogram (out,512,64,512, 44100, 'yaxis');
end