% close all;
clear all;

drawThings = false;
drawSpeed = 1;
drawStart = 0;

numFromBoundX = 1;
numFromBoundY = 1;

gridSim = false; 
lockY = false;
dispCorr = true;
waveSpeed = false;

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
if gridSim
    lengthSound = 100^2;
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

NxStart = 5.0;
NxEnd = 6.0;
NyStart = 5.0;
NyEnd = 5.0;

% h = 2 * sqrt(k * (sig1 + sqrt(sig1^2 + kappaSqVec(1)))); 
% hVec = 2 * sqrt(k * (sig1 + sqrt(sig1^2 + kappaSqVec))); 
h = sqrt(cVec(1)^2 * k^2 + 4 * sig1 * k + sqrt((cVec(1)^2 * k^2 + 4 * sig1 * k)^2 + 16 * kappaSqVec(1) * k^2));
hVec = sqrt(cVec.^2 * k^2 + 4 * sig1 * k + sqrt((cVec.^2 * k^2 + 4 * sig1 * k).^2 + 16 * kappaSqVec * k^2));

Lx = linspace(NxStart * h, NxEnd * h, lengthSound);
Ly = linspace(NyStart * h, NyEnd * h, lengthSound);

NxMax = ceil(Lx./hVec);
NyMax = ceil(Ly./hVec);
Nmax = max (NxMax .* NyMax);


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

modesSave = zeros (lengthSound, Nmax);
sigmaSave = zeros (lengthSound, Nmax);
centDeviation = zeros (lengthSound, Nmax);
expectedModes = zeros (lengthSound, Nmax);

nCounter = 0;
percentCounter = 0;

if gridSim
    Lx = repmat(linspace(NxStart*h, (NxEnd - 1/sqrt(lengthSound)) * h, sqrt(lengthSound)), sqrt(lengthSound), 1);
    Ly = repmat(linspace(NyStart*h, (NyEnd - 1/sqrt(lengthSound)) * h, sqrt(lengthSound))', 1, sqrt(lengthSound));
    Lx = reshape(Lx, lengthSound, 1);
    Ly = reshape(Ly, lengthSound, 1);
    LxIdx = 1;
    LyIdx = 1;
    maxCentDeviation = zeros(sqrt(lengthSound), sqrt(lengthSound));

end
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
        sigDC = 0.01;
        eps = 1e-15;
        betaXDC = (1 - alfX) / (alfX + eps);
        betaYDC = (1 - alfY) / (alfY + eps);
        JDCX = sparse(zeros (My + Mwy, Mx + Mwx));
        JDCX(:, Mx) = 	1/h^2;
        JDCX(:, Mx+1) = -1/h^2;
        JDCX(My:My+1, Mx:Mx+1) = zeros(2);

        JDCY = sparse(zeros (My + Mwy, Mx + Mwx));
        JDCY(My, :) = 1/h^2;
        JDCY(My+1, :) = -1/h^2;
        JDCY(My:My+1, Mx:Mx+1) = zeros(2);

        
        etaDCX = sparse(zeros (My + Mwy, Mx + Mwx));
        etaDCX(:, Mx) = -1;
        etaDCX(:, Mx+1) = 1;
        etaDCX(My:My+1, Mx:Mx+1) = zeros(2);

        etaDCY = sparse(zeros (My + Mwy, Mx + Mwx));
        etaDCY(My, :) = -1;
        etaDCY(My+1, :) = 1;
        etaDCY(My:My+1, Mx:Mx+1) = zeros(2);
        
        JDCX = reshape(JDCX, (My + Mwy) * (Mx + Mwx), 1);
        JDCY = reshape(JDCY, (My + Mwy) * (Mx + Mwx), 1);
        etaDCX = reshape(etaDCX, (My + Mwy) * (Mx + Mwx), 1)';
        etaDCY = reshape(etaDCY, (My + Mwy) * (Mx + Mwx), 1)';
            
        Adispcorr =  -betaXDC * k^2 * (1 + sigDC / k) / 2 * JDCX * etaDCX;% - ...
%             betaYDC * k^2 * (1 + sigDC / k) / 2 * JDCY * etaDCY;
        Cdispcorr = betaXDC * k^2 * (1 - sigDC / k) / 2 * JDCX * etaDCX + ...
%             betaYDC * k^2 * (1 - sigDC / k) / 2 * JDCY * etaDCY;
        Amat = Amat + Adispcorr;
        C = C + Cdispcorr;
    end
    numExpectedModes = (Mx+Mwx) * (My+Mwy);
    [f, sigma, ~] = analyseQ([Amat\B, Amat\C; ...
              eye(size(B)), zeros(size(B))], k);
    f = f(1:numExpectedModes);
    modesSave(n, 1:numExpectedModes) = f(1:numExpectedModes);
    sigmaSave(n, 1:numExpectedModes) = sigma(1:numExpectedModes);
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
    if waveSpeed
        lambda = cVec(n) * k / h;
        expectedF = 1/(pi*k) * asin(lambda * sqrt((SxMat + SyMat)));
    else
        mu = sqrt(kappaSqVec(n)) * k / h^2;
        expectedF = 1/(pi*k) * asin(2 * mu * (SxMat + SyMat));
    end
    expectedF = real(sort(reshape(expectedF, 1, (Mx+Mwx) * (My+Mwy))));
    expectedModes(n, 1:numExpectedModes) = expectedF;
    centDeviation(n, 1:numExpectedModes) = 1200 * log2(f./expectedF');
    if gridSim
        idx = n - 1 + (n == 1);
        if Lx(n) > Lx(idx)
            LxIdx = LxIdx + 1;
        end
        if Ly(n) < Ly(idx)
            LyIdx = 1;
        end

        maxCentDeviation(LyIdx, LxIdx) = min(1200 * log2(f./expectedF'));
        LyIdx = LyIdx + 1;
    end
    if drawThings && mod (n,drawSpeed) == 0
        subplot(211);
        hold off;
        plot(f)
        hold on;
        plot(expectedF);

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
if gridSim
    surf(maxCentDeviation, 'Linestyle', 'none')
end
%%
figure
modesSave(modesSave == 0) = nan;
for i = 1:length(nChangeSave)-1
    plot(nChangeSave(i):nChangeSave(i+1)-1, modesSave(nChangeSave(i):nChangeSave(i+1)-1, :), 'k')
    hold on;
%     plot(nChangeSave(i):nChangeSave(i+1)-1, expectedModes(nChangeSave(i):nChangeSave(i+1)-1, :), 'r')
%     hold on;

end
%%
figure
plotOneStepAnalysis
%%
% figure
% plot(centDeviation)
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