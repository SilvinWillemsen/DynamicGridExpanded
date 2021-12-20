%{
    Dynamic stiff string...
%}

clear all;
close all;
clc;

%% Draw settings
fs = 44100;             % Sample rate (in Hz)
k = 1/fs;               % Time step (in s)
lengthSound = 1000;   % Length of the simulation (in samples)
maxEig = false;
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

dispCorr = false;
changeN = false;

Nstart = 15;
Nend = 20;

Nvec = linspace(Nstart, Nend, lengthSound);

changeWaveSpeed = true;

%% Initialise parameters speeds for the entire simulation
Lvec = linspace(1, 1, lengthSound);   
rho = 7850;
r = 5e-4;
A = pi * r^2; 
I = pi / 4 * r^4;

Evec = linspace(sqrt(0), sqrt(0), lengthSound).^2;
E = Evec(1);
kappaSq = E * I / (rho * A);
if changeWaveSpeed
    kappaSqVec = Evec * I / (rho * A);
else 
%     kappaSqVec = linspace(sqrt(sqrt((1/(2 * k * Nstart^2)))), sqrt(sqrt((1/(2 * k * Nend^2)))), lengthSound).^8;
    if changeN
        kappaSqVec = (Lvec.^2./(2 * k * Nvec.^2)).^2;
    else
        %% is this right?
        kappaSqVec = linspace(sqrt(sqrt((1/(2 * k * Nstart^2)))), sqrt(sqrt((1/(2 * k * Nend^2)))), lengthSound).^8;
    end
end

sig0 = 0;
sig1 = 0.000;


if changeWaveSpeed
    cStart = sqrt((Lvec(1)^2 - 4 * kappaSqVec(1) * k^2 * Nstart^4) / (k^2 * Nstart^2));
    cEnd = sqrt((Lvec(end)^2 - 4 * kappaSqVec(end) * k^2 * Nend^4) / (k^2 * Nend^2));
    if changeN
        cVec = sqrt((Lvec(1)^2 - 4 * kappaSqVec(1) * k^2 * Nvec.^4) ./ (k^2 * Nvec.^2));
    else
        cVec = linspace(cStart, cEnd, lengthSound);
    end
else
    cVec = linspace(0, 0, lengthSound);
end
c = cVec(1);
kappaSq = kappaSqVec(1);

h = sqrt((c^2 * k^2 + 4 * sig1 * k + sqrt((c^2 * k^2 + 4 * sig1 * k)^2 + 16 * kappaSq * k^2))/2);

Nfrac = Lvec(1) / h;
N = floor(round(Nfrac * 10000)/10000);
NPrev = N;

% calculate maximum number of points
hVec = sqrt((cVec.^2 * k^2 + 4 * sig1 * k + sqrt((cVec.^2 * k^2 + 4 * sig1 * k).^2 + 16 * kappaSqVec * k^2))/2);
Nmax = ceil(max(Lvec ./ hVec));


modesSave = zeros (lengthSound, Nmax);
sigmaSave = zeros (lengthSound, Nmax);
zSave = zeros (lengthSound, Nmax);
expectedF = zeros (lengthSound, Nmax);
centDeviation = zeros (lengthSound, Nmax);
hzDeviation = zeros (lengthSound, Nmax);

% Lvec = linspace(floor(Nfrac) * h, ceil(Nfrac)*h, lengthSound);

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

Dxpu = spdiags([-eu, eu], 0:1, M, M);
Dxpw = spdiags([-eu, eu], 0:1, Mw, Mw);

        
%% Initialise output
out = zeros(floor(lengthSound), 1);

nn = 2;
nChangeSave = [1];
Nsave = [Nstart];
%% main loop
for n = 1:lengthSound  

    c = cVec(n);
    E = Evec(n);
    L = Lvec(n);
%     kappaSq = E * I / (rho * A);
    kappaSq = kappaSqVec(n);

    % recalculate gridspacing and (fractional) number of intervals
    h = sqrt((c^2 * k^2 + 4 * sig1 * k + sqrt((c^2 * k^2 + 4 * sig1 * k)^2 + 16 * kappaSq * k^2))/2);
    Nfrac = L / h;
    N = floor (round (Nfrac * 10000)/10000);
    NsaveLoop(n) = Nfrac;
    lambdaSq = c^2 * k^2 / h^2; 
    muSq = kappaSq * k^2 / h^4; % should always be 1

    % Calculate alpha (fractional part of Nfrac)
    alf = Nfrac - N;
    alfSave(n) = alf;
    % Can only add/remove one point at a time
%     if abs(N - NPrev) > 1
%         error('Can only add/remove one grid point at a time')
%     end
    
    %% Check whether to add or remove points
    if N ~= NPrev
        nChangeSave(nn) = n + 1;
        Nsave(nn) = N;
        nn = nn+1;
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
        
        Dxpu = spdiags([-eu, eu], 0:1, M, M);
        Dxpw = spdiags([-eu, eu], 0:1, Mw, Mw);

        Dxxu = spdiags([eu -2*eu eu], -1:1, M,M);
        ew = ones(Mw, 1);
        Dxxw = spdiags([ew -2*ew ew], -1:1, Mw, Mw);
    end
    
%     Dxp = zeros(M + Mw);
%     Dxp(1:M, 1:M) = Dxpu;
%     Dxp(M+1:end, M+1:end) = Dxpw;

    Dxx = zeros(M + Mw);
    Dxx(1:M, 1:M) = Dxxu;
    Dxx(M+1:end, M+1:end) = Dxxw;

    %% Calculate (quadratic) interpolator and virtual grid points
    ip = [-(alf - 1) / (alf + 1), 1, (alf - 1) / (alf + 1)];

    if numFromBound == 1
%         Dxp(M, M) = Dxp(M, M) + ip(3);
%         Dxp(M, M+1) = 1/3;
%         Dxp(M+1, M) = 1/3;
%         Dxx = -Dxp' * Dxp;
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
    
    if sig0 == 0 && sig1 == 0 && ~maxEig
        if kappaSq == 0
            f = sort(real(1/(pi * k) * asin(c * k / 2 * sqrt(-eig(Dxx)))));
        elseif c == 0
            f = sort(real(1/(pi * k) * asin(sqrt(kappaSq) * k / 2 * sqrt(eig(Dxxxx)))));
        else
            f = sort(real(1/(pi * k) * asin(1/2 * sqrt(-eig(c^2*k^2 * Dxx - kappaSq * k^2 * Dxxxx)))));

        end
        modesSave(n, 1:(M+Mw)) = f(1:(M+Mw));

    else
        Q = [Amat\B, Amat\C; ...
                  eye(size(B)), zeros(size(B))];
        [f, sigma, zSav] = analyseQ(Q, k);
%         zSave(n, 1:(M+Mw)) = zSav;
        modesSave(n, 1:(M+Mw)) = f(1:(M+Mw));
        sigmaSave(n, 1:(M+Mw)) = sigma(1:(M+Mw));
    end
    %     subplot(211)
%     hold off;
%     plot(f(1:(M+Mw)))
%     drawnow;
%     hold on;
%     beta = (1:length(f))*pi; % bar 
    p = 1:(M+Mw);
    if kappaSq ~= 0
%         bCoeff = c^2 - 2 * sig0 * sig1;
%         aCoeff = kappaSq - sig1^2;
%         cCoeff = -(sig0 + c^2 * p.^2 * pi^2 / L^2 + kappaSq * pi^4 * p.^4 / L^4);
%         beta = sqrt((-bCoeff + sqrt(bCoeff^2 - 4 * aCoeff * cCoeff)) / (2 * aCoeff));

        % gives machine imprecision for lower values of kappaSq
        if c == 0
            beta = pi * p / L; % might be wrong..?
        else
            beta = sqrt(((sqrt(4 * pi^4 * kappaSq^2 * p.^4 / L^4 + c^4 + 4 * c^2 * pi^2 * kappaSq * p.^2 / L^2)) - c^2)/(2*kappaSq));
        end
    else
        beta = p * pi / L;
    end
    expectedF(n, 1:(M+Mw)) = 1/(pi*k) * asin(sqrt(lambdaSq * sin(beta * h / 2).^2 + 4 * muSq * sin(beta*h / 2).^4));
%     plot(expectedF);
    
%     subplot(212)
    centDeviation(n, 1:(M+Mw)) = 1200 * log2(f(1:(M+Mw))./expectedF(n, 1:(M+Mw))');
%     centDeviation(n, 1:(M+Mw)) = f(1:(M+Mw))-expectedF(n, 1:(M+Mw))';
%     plot (centDeviation(n, 1:(M+Mw)));
%     drawnow;
   

    NPrev = N;
    
end
%%
figure;
modesSave(modesSave == 0) = nan;
expectedF(expectedF == 0) = nan;
% plot(modesSave(1:n, :))

plot (Nvec, centDeviation, 'Linewidth', 1);
xlabel("$\mathcal{N}$ (fractional no. of intervals)", 'interpreter', 'latex')
ylabel("Deviation (in cents)", 'interpreter', 'latex')
xticks(Nsave);
set(gca, 'Fontsize', 18)
grid on;
hold on;

%%

dxdy = zeros(max(Nsave)-1, (Nend-Nstart));
diffs = zeros(max(Nsave)-1, (Nend-Nstart));
diffX = zeros(max(Nsave)-1, (Nend-Nstart));
diffY = zeros(max(Nsave)-1, (Nend-Nstart));
percentageMinOfX = zeros(max(Nsave), (Nend-Nstart));
minN = zeros(max(Nsave), (Nend-Nstart));
minNVal = zeros(max(Nsave), (Nend-Nstart));
for NN = 1:(Nend-Nstart)
    for i = 1:Nsave(NN)-1
        xLoc1 = find(centDeviation(nChangeSave(NN):nChangeSave(NN+1)-1, i)==min(centDeviation(nChangeSave(NN):nChangeSave(NN+1)-1, i))) ...
            + nChangeSave(NN) - 1;
        xLoc2 = find(centDeviation(nChangeSave(NN):nChangeSave(NN+1)-1, i+1)==min(centDeviation(nChangeSave(NN):nChangeSave(NN+1)-1, i+1))) ...
            + nChangeSave(NN) - 1;

        yLoc1 = centDeviation(xLoc1, i);
        yLoc2 = centDeviation(xLoc2, i+1);
        
        if i == 1
            lowestFreqMin(NN) = (xLoc1 - nChangeSave(NN)) / (nChangeSave(NN+1)-1 - nChangeSave(NN));
        end
        minN(i, NN) = NsaveLoop(xLoc1);
        minNVal(i, NN) = yLoc1;
        percentageMinOfX(i, NN) = (xLoc1 - nChangeSave(NN)) / (nChangeSave(NN+1)-1 - nChangeSave(NN));

        if i == Nsave(NN)-1
            highestFreqMin(NN) = (xLoc2 - nChangeSave(NN)) / (nChangeSave(NN+1)-1 - nChangeSave(NN));
            NfracSave(NN) = Nvec(xLoc2);
            percentageMinOfX(i+1, NN) = (xLoc2 - nChangeSave(NN)) / (nChangeSave(NN+1)-1 - nChangeSave(NN));
            minN(i+1, NN) = NsaveLoop(xLoc2);
            minNVal(i+1, NN) = yLoc2;
            minValLast(NN) = yLoc2;
        end
        
        dxdy(i, NN) = (yLoc2 - yLoc1) / (Nvec(xLoc2) - Nvec(xLoc1));
        diffX(i, NN) = Nvec(xLoc2) - Nvec(xLoc1);
        diffY(i, NN) = yLoc2 - yLoc1;
        diffs(i, NN) = sqrt(diffY(i, NN)^2 + diffX(i, NN)^2);
        plot([Nvec(xLoc1), Nvec(xLoc2)], [yLoc1, yLoc2], 'k');
        hold on;
    end
    dxdyMean(NN) = mean(dxdy(i-8:i, NN));
    lastDxDy(NN) = dxdy(i, NN);
    lastMinN(NN) = 0.5 * (minN(i, NN) + minN(i+1, NN));
    minNMean(NN) = mean(minN(i-8:i+1, NN));
    minSave(NN) = min(min(centDeviation(nChangeSave(NN):nChangeSave(NN+1)-1, :)));
end
% ylim([-120, max(max(centDeviation))]);

figure('Position', [489 542 560 315]);
% 
plot(real(expectedF(1:n, :)), '--', 'color', [1, 0, 0, 0.5], 'Linewidth', 1.5)
hold on
plot(modesSave(1:n, :), 'k', 'Linewidth', 1.5)
hold on;
xticks(nChangeSave);
xticklabels(string(Nstart:Nend));
% figure
title ("Modal Analysis $\mathcal{N} = " + Nstart + " \rightarrow" + Nend + "$", 'interpreter', 'latex');
xlabelsave = num2cell(Nstart:sign(Nend-Nstart):Nend);
set(gca, 'Linewidth', 2, 'Fontsize', 16, 'TickLabelInterpreter', 'latex', ...
    'Position', [0.0929 0.1100 0.8696 0.8086])
xLab = xlabel("$\mathcal{N}$", 'interpreter', 'latex');
ylabel("Frequency (Hz)", 'interpreter', 'latex')
xLab.Position(1) = 567;
xLab.Position(2) =-1000;
ylim([0, fs / 2])
xlim([0, lengthSound+1])
grid on

% minNVal(minNVal>=0) = nan;
% plot(minNVal)
% %%
% figure
% subplot(311)
% diffs(diffs<=0) = nan;
% plot(diffs)
% 
% subplot(312)
% diffX(diffX>=0) = nan;
% plot(diffX)
% 
% subplot(313)
% diffY(diffY>=0) = nan;
% plot(diffY)

%%
figure
dxdy(dxdy<=0) = nan;
plot(fliplr(dxdy .* minN(1:end-1, :)))

