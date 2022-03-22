% close all;
clear all;

drawThings = false;
drawSpeed = 1;
drawStart = 0;

numFromBoundX = 1;
numFromBoundY = 1;

gridSim = false; 
lockY = false;
dispCorr = false;
waveSpeed = false;

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]

NxStart = 15.0;
NxEnd = 20.0;
NyStart = 15.0;
NyEnd = 20.0;

    
if gridSim
    detail = 50;
    NxDiff = NxEnd-NxStart;
    NyDiff = NyEnd-NyStart;
    lengthSound = NxDiff * detail * NyDiff * detail; % detail x detail per Nx Ny combination
%     lengthSound = 100^2;
else
    lengthSound = 1000;   % Length of the simulation
end
rho = 7850;
H = 0.005;
if waveSpeed
    cVec = linspace(1000, 1000, lengthSound);
    Evec = linspace(sqrt(0), sqrt(0), lengthSound).^2;

else
    cVec = linspace(0, 0, lengthSound);
    Evec = linspace(sqrt(2e11), sqrt(2e11), lengthSound).^2;
end
nu = 0.3;
Dvar = Evec * H^3 / (12 * (1 - nu^2));
kappaSqVec = Dvar / (rho * H);

sig0 = 0;
sig1 = 0.000;

% h = 2 * sqrt(k * (sig1 + sqrt(sig1^2 + kappaSqVec(1)))); 
% hVec = 2 * sqrt(k * (sig1 + sqrt(sig1^2 + kappaSqVec))); 
h = sqrt(cVec(1)^2 * k^2 + 4 * sig1 * k + sqrt((cVec(1)^2 * k^2 + 4 * sig1 * k)^2 + 16 * kappaSqVec(1) * k^2));
hVec = sqrt(cVec.^2 * k^2 + 4 * sig1 * k + sqrt((cVec.^2 * k^2 + 4 * sig1 * k).^2 + 16 * kappaSqVec * k^2));

if gridSim
    Lx = 1./linspace(1/(NxStart * h), 1/(NxEnd * h), NxDiff * detail);
    Ly = 1./linspace(1/(NyStart * h), 1/(NyEnd * h), NyDiff * detail);
else
    Lx = 1./linspace(1/(NxStart * h), 1/(NxEnd * h), lengthSound);
    Ly = 1./linspace(1/(NyStart * h), 1/(NyEnd * h), lengthSound);

end
NxFrac = Lx(1)/h;
Nx = floor(NxFrac);
NxPrev = Nx;

if numFromBoundX == -1
    Mx = ceil(Nx * 0.5);
    Mwx = floor(Nx * 0.5);
else
    Mx = Nx-numFromBoundX;
    Mwx = numFromBoundX;
end

NyFrac = Ly(1)/h;
Ny = floor(NyFrac);
NyPrev = Ny;
if lockY
    Ly = ones(1, length(Ly)) * Ny * h;
end

if numFromBoundY == -1
    My = ceil(Ny * 0.5);
    Mwy = floor(Ny * 0.5);
else
    My = Ny-numFromBoundY;
    Mwy = numFromBoundY;
end

nCounter = 0;
percentCounter = 0;

if gridSim
    Lx = repmat(linspace(NxStart*h, (NxEnd - 1/(NxDiff * detail)) * h, NxDiff * detail), NyDiff * detail, 1);
    Ly = repmat(linspace(NyStart*h, (NyEnd - 1/(NyDiff * detail)) * h, NyDiff * detail)', 1, NxDiff * detail);
    Lx = reshape(Lx, lengthSound, 1);
    Ly = reshape(Ly, lengthSound, 1);
    LxIdx = 1;
    LyIdx = 1;
    maxCentDeviation = zeros(NyDiff * detail, NyDiff * detail);

end

% assuming h doesn't change
NxMax = ceil(Lx./h);
NyMax = ceil(Ly./h);
Nmax = max (NxMax .* NyMax);
modesSave = zeros (lengthSound, Nmax);
sigmaSave = zeros (lengthSound, Nmax);
centDeviation = zeros (lengthSound, Nmax);
centDeviationSquare = zeros (lengthSound, Nmax);
expectedModes = zeros (lengthSound, Nmax);

% alfXMat = repmat((0:sqrt(lengthSound):lengthSound-1)/lengthSound,sqrt(lengthSound),1);
% 
% alfYMat = alfXMat';
% alfXVec = reshape(alfXMat, lengthSound, 1);
% alfYVec = reshape(alfYMat, lengthSound, 1);

%% Simulation loop
nChangeSave = [1];
for n = 1:lengthSound

    h = sqrt(cVec(1)^2 * k^2 + 4 * sig1 * k + sqrt((cVec(1)^2 * k^2 + 4 * sig1 * k)^2 + 16 * kappaSqVec(1) * k^2));
    NxFrac = Lx(n)/h;
    Nx = floor(NxFrac);

    alfX = NxFrac - Nx;
    alfXSave(n) = alfX;
%     alfX = alfXVec(n);
    if Nx ~= NxPrev
        nChangeSave = [nChangeSave; n];

        if numFromBoundX == -1
            Mx = ceil(Nx * 0.5);
            Mwx = floor(Nx * 0.5); % CHANGING Mx AND Mwx HERE!!
        else
            Mx = Nx - numFromBoundX;
            Mwx = numFromBoundX;
        end
    end
    NxPrev = Nx;
    
    NyFrac = Ly(n)/h;
    Ny = floor(NyFrac);

    alfY = NyFrac - Ny;
%     alfY = alfYVec(n);
    alfYSave(n) = alfY;

    if Ny ~= NyPrev
        nChangeSave = [nChangeSave; n];
        if numFromBoundY == -1
            My = ceil(Ny * 0.5);
            Mwy = floor(Ny * 0.5); % CHANGING Mx AND Mwx HERE!!
        else
            My = Ny - numFromBoundY;
            Mwy = numFromBoundY;
        end
    end
    NyPrev = Ny;
    
%     alfX = 1;
    ipX = [-(alfX - 1) / (alfX + 1), 1, (alfX - 1) / (alfX + 1)];
    
    Dxx = zeros(Nx);
    euX = ones(Mx, 1);
    ewX = ones(Mwx, 1);
    Dxx(1:Mx, 1:Mx) = spdiags([euX -2*euX euX], -1:1, Mx, Mx);
    Dxx(Mx+1:end, Mx+1:end) = spdiags([ewX -2*ewX ewX], -1:1, Mwx, Mwx);
    
%     Dxx(1:Mx, 1:Mx) = toeplitz([-2, 1, zeros(1, Mx-2)]);
%     Dxx(Mx+1:end, Mx+1:end) = toeplitz([-2, 1, zeros(1, Mwx-2)]);
    
    if numFromBoundX == 1
        Dxx(Mx, Mx:(Mx+1)) = Dxx(Mx, Mx:(Mx+1)) + fliplr(ipX(2:3));
    else
        Dxx(Mx, Mx:(Mx+2)) = Dxx(Mx, Mx:(Mx+2)) + fliplr(ipX);
    end
    Dxx(Mx+1, (Mx-1):(Mx+1)) = Dxx(Mx+1, (Mx-1):(Mx+1)) + ipX;

    Dxx = Dxx / h^2;
    
    
    ipY = [-(alfY - 1) / (alfY + 1), 1, (alfY - 1) / (alfY + 1)];

    Dyy = zeros(Ny);
    euY = ones(My, 1);
    ewY = ones(Mwy, 1);
    Dyy(1:My, 1:My) = spdiags([euY -2*euY euY], -1:1, My, My);
    Dyy(My+1:end, My+1:end) = spdiags([ewY -2*ewY ewY], -1:1, Mwy, Mwy);
    
%     Dyy(1:My, 1:My) = toeplitz([-2, 1, zeros(1, My-2)]);
%     Dyy(My+1:end, My+1:end) = toeplitz([-2, 1, zeros(1, Mwx-2)]);
    
    if Ny > 2
        if numFromBoundY == 1
            Dyy(My, My:(My+1)) = Dyy(My, My:(My+1)) + fliplr(ipY(2:3));
        else
            Dyy(My, My:(My+2)) = Dyy(My, My:(My+2)) + fliplr(ipY);
        end
        Dyy(My+1, (My-1):(My+1)) = Dyy(My+1, (My-1):(My+1)) + ipY;
    end
    Dyy = Dyy / h^2;
    
    D = kron(speye(Nx), Dyy) + kron(Dxx, speye(Ny));
    DD = D * D;
    B = 2 * speye(Ny*Nx) + cVec(n)^2 * k^2 * D - kappaSqVec(n) * k^2 * DD + 2 * sig1 * k * D;
    Amat = speye(Ny*Nx) * (1 + sig0 * k);
    C = (-1 + sig0 * k) * speye(Ny*Nx) - 2 * sig1 * k * D;
        
    if dispCorr
        epsilon = 1e-6;
        sigDC = 10;
        if alfX == 0
            betaX = 0;
        else
            betaX = (1-alfX) / (alfX + epsilon);
        end
        if alfY == 0
            betaY = 0;
        else
            betaY = (1-alfY) / (alfY + epsilon);

        end
        JDCXEtaSave = sparse (zeros ((My + Mwy) * (Mx + Mwx)));
        
        for yLoc = 1:My+Mwy
            JDCX = sparse(zeros (My + Mwy, Mx + Mwx));     
            etaDCX = sparse(zeros (My + Mwy, Mx + Mwx));

            if yLoc ~= My && yLoc ~= My+1
                JDCX(yLoc, Mx) = 1/h^2;
                JDCX(yLoc, Mx+1) = -1/h^2;
                etaDCX(yLoc, Mx) = -1;
                etaDCX(yLoc, Mx+1) = 1;
            end
            JDCX = reshape(JDCX, (My + Mwy) * (Mx + Mwx), 1);
            etaDCX = reshape(etaDCX, (My + Mwy) * (Mx + Mwx), 1)';
            JDCXEtaSave = JDCXEtaSave + JDCX * etaDCX;
        end
        
        JDCYEtaSave = sparse (zeros ((My + Mwy) * (Mx + Mwx)));

        for xLoc = 1:Mx+Mwx
            JDCY = sparse(zeros (My + Mwy, Mx + Mwx));     
            etaDCY = sparse(zeros (My + Mwy, Mx + Mwx));

            if xLoc ~= Mx && xLoc ~= Mx+1
                JDCY(My, xLoc) = 1/h^2;
                JDCY(My+1, xLoc) = -1/h^2;
                etaDCY(My, xLoc) = -1;
                etaDCY(My+1, xLoc) = 1;
            end
            JDCY = reshape(JDCY, (My + Mwy) * (Mx + Mwx), 1);
            etaDCY = reshape(etaDCY, (My + Mwy) * (Mx + Mwx), 1)';
            JDCYEtaSave = JDCYEtaSave + JDCY * etaDCY;
        end

        Adispcorr =  -betaX * k^2 * (1 + sigDC / k) / 2 * JDCXEtaSave - ...
            betaY * k^2 * (1 + sigDC / k) / 2 * JDCYEtaSave;
        Cdispcorr = betaX * k^2 * (1 - sigDC / k) / 2 * JDCXEtaSave + ...
            betaY * k^2 * (1 - sigDC / k) / 2 * JDCYEtaSave;

        Amat = Amat + Adispcorr;
        C = C + Cdispcorr;
    end
    numExpectedModes = (Mx+Mwx) * (My+Mwy);
%     [f, sigma, ~] = analyseQ([Amat\B, Amat\C; ...
%               eye(size(B)), zeros(size(B))], k);
    f = sort(real(1/(2*pi*k) * acos(1/2 * eig(full(B)))));
    if numExpectedModes ~= length(f)
        disp(length(f) - numExpectedModes)
    end
    numExpectedModes = min(numExpectedModes, length(f));
    f = f(1:numExpectedModes);
    modesSave(n, 1:numExpectedModes) = f(1:numExpectedModes);
%     sigmaSave(n, 1:numExpectedModes) = sigma(1:numExpectedModes);
    %     subplot(211)
%     hold off;
%     plot(f(1:(M+Mw)))
%     hold on;
%     beta = (1:length(f))*pi; % bar 
    px = 1:(Mx+Mwx);
    py = 1:(My+Mwy);
    
    betaX = px * pi / Lx(n);
    betaY = py * pi / Ly(n);
    SxMat = repmat(sin(betaX * h / 2).^2, My+Mwy, 1);
    SyMat = repmat(sin(betaY' * h / 2).^2, 1, Mx+Mwx);
    
    SxSquare = sin(betaX * h / 2).^2;
    SySquare = sin(betaY * h / 2).^2;
    if waveSpeed
        lambda = cVec(n) * k / h;
        expectedModes = 1/(pi*k) * asin(lambda * sqrt((SxMat + SyMat)));
    else
        mu = sqrt(kappaSqVec(n)) * k / h^2;
        expectedModes = 1/(pi*k) * asin(2 * mu * (SxMat + SyMat));
    end
    expectedModes = real(sort(reshape(expectedModes, 1, (Mx+Mwx) * (My+Mwy))));
    expectedModesSquare = real(1/(pi*k) * asin(2 * mu * (SxSquare + SySquare)));
    indexSave = [];
    for testIdx = 1:length(SxSquare)
        for testIdx2 = 1:length(expectedModes)
            if expectedModesSquare(testIdx) == expectedModes(testIdx2)
                indexSave = [indexSave; testIdx2];
                indexSaved = true;
                break;
            end   
        end
        if indexSaved
            continue;
        end
    end
    expectedF(n, 1:numExpectedModes) = expectedModes(1:numExpectedModes);
    expectedFSquare(n, 1:length(SxSquare)) = expectedModesSquare;
    modesSaveSquare(n, 1:length(SxSquare)) = f(indexSave);

    centDeviation(n, 1:numExpectedModes) = 1200 * log2(f./expectedModes');
    centDeviationSquare(n, 1:length(SxSquare)) = 1200 * log2(f(indexSave)./expectedModes(indexSave)');
    if gridSim
        idx = n - 1 + (n == 1);
        if Lx(n) > Lx(idx)
            LxIdx = LxIdx + 1;
        end
        if Ly(n) < Ly(idx)
            LyIdx = 1;
        end

        maxCentDeviation(LyIdx, LxIdx) = min(1200 * log2(f./expectedModes'));
        LyIdx = LyIdx + 1;
    end
    if drawThings && mod (n,drawSpeed) == 0
        subplot(211);
        hold off;
        plot(f)
        hold on;
        plot(expectedModes);

        subplot(212)
        plot (centDeviation(n, 1:numExpectedModes));
        drawnow;
    end
    if n > lengthSound * nCounter / 100
        percentCounter = percentCounter + 1;
        nCounter = nCounter + 1;
        disp((percentCounter) + "% done")
    end
    
end
%%
figure
% disp("plotting One Step Analysis")
% plotOneStepAnalysis
% return;
if gridSim
    surf(maxCentDeviation, 'Linestyle', 'none')
end
%%
figure
nChangeSave = [nChangeSave; n+1];
modesSave(modesSave == 0) = nan;
expectedF(expectedF == 0) = nan;
for i = 1:length(nChangeSave)-1
    plot(nChangeSave(i):nChangeSave(i+1)-1, modesSave(nChangeSave(i):nChangeSave(i+1)-1, :), 'k')
    hold on;
%     plot(nChangeSave(i):nChangeSave(i+1)-1, expectedF(nChangeSave(i):nChangeSave(i+1)-1, :), 'r')
    hold on;

end
%%
figure
plot(centDeviation)
%%
% for i = 1:sqrt(lengthSound)
%     plotRange = (i-1) * sqrt(lengthSound) + 1: i * sqrt(lengthSound);
%     plot(centDeviation(plotRange, :));
%     % find mode that contains minimum
%     [~, b] = find(centDeviation(plotRange, :) == min(min(centDeviation(plotRange, :))));
% 
%     % find minimum using interpolation (for higher quality)
%     detail = 10;
%     interpVect = interp1(1:sqrt(lengthSound),centDeviation(plotRange, b),1:1/detail:sqrt(lengthSound),'cubic');
%     idx = find(interpVect == min(interpVect));
%     minSave(i) = (idx / detail) / sqrt(lengthSound);
%     drawnow;
% end
% plot(minSave)