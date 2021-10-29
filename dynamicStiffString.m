%{
    Dynamic stiff string...
%}

clear all;
% close all;
clc;

%% Draw settings
drawThings = true;
drawStart = 0;
drawSpeed = 10;
plotModalAnalysis = false;

fs = 44100;             % Sample rate (in Hz)
k = 1/fs;               % Time step (in s)
if plotModalAnalysis
    lengthSound = 1000;   % Length of the simulation (in samples)
else
    lengthSound = fs * 1;
end

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
%% Initialise wave speeds for the entire simulation
startParams = [1; 7850; 5e-4; 300; 2e11; 1; 0.005];
% endParams = startParams;
endParams = [1; 7850; 5e-4; 300; 2e12; 1; 0.005];
Lvec = linspace(startParams(1), endParams(1), lengthSound);
rhoVec = linspace(startParams(2), endParams(2), lengthSound);
rVec = linspace(startParams(3), endParams(3), lengthSound);
Tvec = linspace(startParams(4), endParams(4), lengthSound);
Evec = linspace(startParams(5), endParams(5), lengthSound);
sig0Vec = linspace(startParams(6), endParams(6), lengthSound);
sig1Vec = linspace(startParams(7), endParams(7), lengthSound);

% Evec = linspace(sqrt(2e11), sqrt(2e14), lengthSound*1/4).^2;
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

%% Matrices
DxxFull = zeros(M + Mw);

eu = ones(M, 1);
ew = ones(Mw, 1);
DxxFull(1:M, 1:M) = spdiags([eu -2*eu eu], -1:1, M, M);
DxxFull(M+1:end, M+1:end) = spdiags([ew -2*ew ew], -1:1, Mw, Mw);
DxxFull = DxxFull / h^2;

DxxxxFull = DxxFull * DxxFull;

        
%% Initialise output
out = zeros(floor(lengthSound), 1);

%% main loop
for n = 1:lengthSound  

    if mod(n, fs/2) == 1
        q(1) = q(1) + 1/h;
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
                qNext = [qNextLeft; cubicIp * [qNextLeft(M-1:M); qNextRight(1:2)]; qNextRight];
                q = [qLeft; cubicIp * [qLeft(M-1:M); qRight(1:2)]; qRight];
                qPrev = [qPrevLeft; cubicIp * [qPrevLeft(M-1:M); qPrevRight(1:2)]; qPrevRight];

            % Otherwise, add grid points to u using interpolation
            elseif numFromBound == 1
                qNext = [qNextLeft; cubicIp(1:3) * [qNextLeft(M-1:M); qNextRight(1)]; qNextRight];
                q = [qLeft; cubicIp(1:3) * [qLeft(M-1:M); qRight(1)]; qRight];
                qPrev = [qPrevLeft; cubicIp(1:3) * [qPrevLeft(M-1:M); qPrevRight(1)]; qPrevRight];

            else
                qNext = [qNextLeft; cubicIp * [qNextLeft(M-1:M); qNextRight(1:2)]; qNextRight];
                q = [qLeft; cubicIp * [qLeft(M-1:M); qRight(1:2)]; qRight];
                qPrev = [qPrevLeft; cubicIp * [qPrevLeft(M-1:M); qPrevRight(1:2)]; qPrevRight];

            end
            disp("Added Point")

        else
        %% Remove point if N^n < N^{n-1}
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
        DxxFull = DxxFull / h^2;
        
        DxxxxFull = DxxFull * DxxFull;
        
    end
    
    NPrev = N;
    alf = 0;
    %% Calculate (quadratic) interpolator and virtual grid points
    ip = [-(alf - 1) / (alf + 1), 1, (alf - 1) / (alf + 1)];
    DxxxxAlfMatrix = [-4, 6, ip(3) - 4, 1, ip(1), 0;
                    1, -4, (ip(3) - 2)^2 + 2, ip(3) - 4, 4 * (alf^2 + alf - 1) / (alf + 1)^2, ip(1)];
    DxxxxAlfMatrix = [DxxxxAlfMatrix; flipud(fliplr(DxxxxAlfMatrix))];
    Dxx = DxxFull;
    Dxxxx = DxxxxFull;
    if numFromBound == 1
        Dxx(M, M:(M+1)) = Dxx(M, M:(M+1)) + fliplr(ip(2:end)) / h^2;
        Dxx(M+1, (M-1):(M+1)) = Dxx(M+1, (M-1):(M+1)) + ip / h^2;
        Dxxxx = Dxx * Dxx;
        
    else
        Dxx(M, M:(M+2)) = Dxx(M, M:(M+2)) + fliplr(ip) / h^2;
        Dxx(M+1, (M-1):(M+1)) = Dxx(M+1, (M-1):(M+1)) + ip / h^2;
%         Dxxxx(M-1:M+2, M-2:M+3) = DxxxxAlfMatrix / h^4;    
        Dxxxx = Dxx * Dxx;
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