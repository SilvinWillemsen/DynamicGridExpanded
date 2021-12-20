% close all;
clear all;

drawThings = false;
energyCalc = true;
plotPropagation = false;
if plotPropagation 
    figure('Position', [440 585 815 213])
end
drawSpeed = 10;

dispCorr = false;
plotLim = 1e-4;
numFromBoundX = -1;
numFromBoundY = -1;
%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs * 2;   % Length of the simulation (1 second) [samples]             
drawStart = 1.5*fs;

rho = 7850;
H = 0.005;
Evec = linspace(sqrt(2e12), sqrt(2e12), lengthSound).^2;
nu = 0.3;
Dvar = Evec * H^3 / (12 * (1 - nu^2));
kappaSqVec = Dvar / (rho * H);

sig0 = 0;
sig1 = 0.000;

Lx = linspace(0.3, 1, lengthSound);

% Lx = 0.5;
Ly = linspace(0.4, 0.4, lengthSound);

h = 2 * sqrt(k * (sig1 + sqrt(sig1^2 + kappaSqVec(1)))); 

Nx = floor(Lx(1)/h);
NxPrev = Nx;

if numFromBoundX == -1
    Mx = ceil(Nx * 0.5);
    Mwx = floor(Nx * 0.5);
else
    Mx = Nx-numFromBoundX;
    Mwx = numFromBoundX;
end


Ny = floor(Ly(1)/h);
NyPrev = Ny;

if numFromBoundY == -1
    My = ceil(Ny * 0.5);
    Mwy = floor(Ny * 0.5);
else
    My = Ny-numFromBoundY;
    Mwy = numFromBoundY;
end

%% Initialise state vectors (one more grid point than the number of intervals)

% don't subtract from as there is one point for overlap
qNext = zeros(Ny * Nx, 1);
q = zeros(Ny * Nx, 1);

%% Initial conditions (raised cosine)
% halfWidth = floor(min(Nx, Ny) / 10);
halfWidth = 1;
width = 2 * halfWidth + 1;
xInLoc = 1;
yInLoc = 1;

% Implemented to be "number-of-point from right / bottom boundary
xOutLoc = 1;
yOutLoc = 1; 

% Set initial velocity to zero
qPrev = q;

out = zeros(lengthSound, 1);

nCounter = 0;
percentCounter = 0;
%% Simulation loop
for n = 1:lengthSound
    if sig0 == 0 && sig1 == 0
        if n == 1
            q((xInLoc-1) * Ny + yInLoc) = 1;
        end
    elseif mod(n, fs/2) == 1
        q((xInLoc-1) * Ny + yInLoc) = 1;
    end
    
    h = 2 * sqrt(k * (sig1 + sqrt(sig1^2 + kappaSqVec(n)))); 
    NxFrac = Lx(n)/h;
    Nx = floor(NxFrac);

    alfX = NxFrac - Nx;

    if Nx ~= NxPrev
        if Nx > NxPrev
            
            qNextTmp = reshape(qNext, Ny, NxPrev);
            qTmp = reshape(q, Ny, NxPrev);
            qPrevTmp = reshape(qPrev, Ny, NxPrev);

            if numFromBoundX == 1
                qBordersNextX = [qNextTmp(:, Mx-1:Mx), qNextTmp(:, Mx+1)]; % unnecessary
                qBordersX = [qTmp(:, Mx-1:Mx), qTmp(:, Mx+1)];
                qBordersPrevX = [qPrevTmp(:, Mx-1:Mx), qPrevTmp(:, Mx+1)];
            else
                qBordersNextX = [qNextTmp(:, Mx-1:Mx), qNextTmp(:, Mx+1:Mx+2)]; % unnecessary
                qBordersX = [qTmp(:, Mx-1:Mx), qTmp(:, Mx+1:Mx+2)];
                qBordersPrevX = [qPrevTmp(:, Mx-1:Mx), qPrevTmp(:, Mx+1:Mx+2)];
            end
            
            cubicIpX = [alfX * (alfX + 1) / -((alfX + 2) * (alfX + 3)); ...
                2 * alfX / (alfX + 2); ...
                2 / (alfX + 2); ...
                2 * alfX / -((alfX + 3) * (alfX + 2))]';

            if numFromBoundX == -1
                if mod(Nx, 2) == 1 
                    qNext = reshape([qNextTmp(:, 1:Mx), qBordersNextX * cubicIpX', qNextTmp(:, Mx+1:end)], Ny*Nx, 1);
                    q = reshape([qTmp(:, 1:Mx), qBordersX * cubicIpX', qTmp(:, Mx+1:end)], Ny*Nx, 1);
                    qPrev = reshape([qPrevTmp(:, 1:Mx), qBordersPrevX * cubicIpX', qPrevTmp(:, Mx+1:end)], Ny*Nx, 1);
                else
                    qNext = reshape([qNextTmp(:, 1:Mx), qBordersNextX * fliplr(cubicIpX)', qNextTmp(:, Mx+1:end)], Ny*Nx, 1);
                    q = reshape([qTmp(:, 1:Mx), qBordersX * fliplr(cubicIpX)', qTmp(:, Mx+1:end)], Ny*Nx, 1);
                    qPrev = reshape([qPrevTmp(:, 1:Mx), qBordersPrevX * fliplr(cubicIpX)', qPrevTmp(:, Mx+1:end)], Ny*Nx, 1);

                end
            elseif numFromBoundX == 1
                qNext = reshape([qNextTmp(:, 1:Mx), qBordersNextX * cubicIpX(1:3)', qNextTmp(:, Mx+1:end)], Ny*Nx, 1);
                q = reshape([qTmp(:, 1:Mx), qBordersX * cubicIpX(1:3)', qTmp(:, Mx+1:end)], Ny*Nx, 1);
                qPrev = reshape([qPrevTmp(:, 1:Mx), qBordersPrevX * cubicIpX(1:3)', qPrevTmp(:, Mx+1:end)], Ny*Nx, 1);
            else
                qNext = reshape([qNextTmp(:, 1:Mx), qBordersNextX * cubicIpX', qNextTmp(:, Mx+1:end)], Ny*Nx, 1);
                q = reshape([qTmp(:, 1:Mx), qBordersX * cubicIpX', qTmp(:, Mx+1:end)], Ny*Nx, 1);
                qPrev = reshape([qPrevTmp(:, 1:Mx), qBordersPrevX * cubicIpX', qPrevTmp(:, Mx+1:end)], Ny*Nx, 1);
              
            end
           disp("added column at n = " + num2str(n))

        else
            qNextTmp = reshape(qNext, Ny, NxPrev);
            qTmp = reshape(q, Ny, NxPrev);
            qPrevTmp = reshape(qPrev, Ny, NxPrev);
            
            qBordersDiffX = qTmp(:, Mx+1) - qTmp(:, Mx);
            
            subplot(2,1,1)
            plot(qBordersDiffX)
            ylim([-plotLim, plotLim])
            drawnow;
            if numFromBoundY == -1
                if mod(Nx, 2) == 0
                    qNext = reshape([qNextTmp(:, 1:Mx-1), qNextTmp(:, Mx+1:end)], Ny*Nx, 1);
                    q = reshape([qTmp(:, 1:Mx-1), qTmp(:, Mx+1:end)], Ny*Nx, 1);
                    qPrev = reshape([qPrevTmp(:, 1:Mx-1), qPrevTmp(:, Mx+1:end)], Ny*Nx, 1);
                else
                    qNext = reshape([qNextTmp(:, 1:Mx), qNextTmp(:, Mx+2:end)], Ny*Nx, 1);
                    q = reshape([qTmp(:, 1:Mx), qTmp(:, Mx+2:end)], Ny*Nx, 1);
                    qPrev = reshape([qPrevTmp(:, 1:Mx), qPrevTmp(:, Mx+2:end)], Ny*Nx, 1);
                end
            else
                qNext = reshape([qNextTmp(:, 1:Mx-1), qNextTmp(:, Mx+1:end)], Ny*Nx, 1);
                q = reshape([qTmp(:, 1:Mx-1), qTmp(:, Mx+1:end)], Ny*Nx, 1);
                qPrev = reshape([qPrevTmp(:, 1:Mx-1), qPrevTmp(:, Mx+1:end)], Ny*Nx, 1);
            end
            disp("removed column at n = " + num2str(n))

        end
        
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

    if Ny ~= NyPrev
        if Ny > NyPrev
            
            qNextTmp = reshape(qNext, NyPrev, Nx);
            qTmp = reshape(q, NyPrev, Nx);
            qPrevTmp = reshape(qPrev, NyPrev, Nx);

            if numFromBoundY == 1
                qBordersNextY = [qNextTmp(My-1:My, :); qNextTmp(My+1, :)]; % unnecessary
                qBordersY = [qTmp(My-1:My, :); qTmp(My+1, :)];
                qBordersPrevY = [qPrevTmp(My-1:My, :); qPrevTmp(My+1, :)];
            else
                qBordersNextY = [qNextTmp(My-1:My, :); qNextTmp(My+1:My+2, :)]; % unnecessary
                qBordersY = [qTmp(My-1:My, :); qTmp(My+1:My+2, :)];
                qBordersPrevY = [qPrevTmp(My-1:My, :); qPrevTmp(My+1:My+2, :)];
            end
            
            cubicIpY = [alfY * (alfY + 1) / -((alfY + 2) * (alfY + 3)); ...
                2 * alfY / (alfY + 2); ...
                2 / (alfY + 2); ...
                2 * alfY / -((alfY + 3) * (alfY + 2))]';

            if numFromBoundY == -1
                if mod(Ny, 2) == 1 
                    qNext = reshape([qNextTmp(1:My, :); cubicIpY * qBordersNextY; qNextTmp(My+1:end, :)], Ny*Nx, 1);
                    q = reshape([qTmp(1:My, :); cubicIpY * qBordersY; qTmp(My+1:end, :)], Ny*Nx, 1);
                    qPrev = reshape([qPrevTmp(1:My, :); cubicIpY * qBordersPrevY; qPrevTmp(My+1:end, :)], Ny*Nx, 1);
                else
                    qNext = reshape([qNextTmp(1:My, :); fliplr(cubicIpY) * qBordersNextY; qNextTmp(My+1:end, :)], Ny*Nx, 1);
                    q = reshape([qTmp(1:My, :); fliplr(cubicIpY) * qBordersY; qTmp(My+1:end, :)], Ny*Nx, 1);
                    qPrev = reshape([qPrevTmp(1:My, :);  fliplr(cubicIpY) * qBordersPrevY; qPrevTmp(My+1:end, :)], Ny*Nx, 1);

                end
            elseif numFromBoundY == 1
                qNext = reshape([qNextTmp(1:My, :); cubicIpY(1:3) * qBordersNextY; qNextTmp(My+1:end, :)], Ny*Nx, 1);
                q = reshape([qTmp(1:My, :); cubicIpY(1:3) * qBordersY; qTmp(My+1:end, :)], Ny*Nx, 1);
                qPrev = reshape([qPrevTmp(1:My, :); cubicIpY(1:3) * qBordersPrevY; qPrevTmp(My+1:end, :)], Ny*Nx, 1);
            else
                qNext = reshape([qNextTmp(1:My, :); cubicIpY * qBordersNextY; qNextTmp(My+1:end, :)], Ny*Nx, 1);
                q = reshape([qTmp(1:My, :); cubicIpY * qBordersY; qTmp(My+1:end, :)], Ny*Nx, 1);
                qPrev = reshape([qPrevTmp(1:My, :); cubicIpY * qBordersPrevY; qPrevTmp(My+1:end, :)], Ny*Nx, 1);
              
            end
           	disp("added row at n = " + num2str(n))

        else
            qNextTmp = reshape(qNext, NyPrev, Nx);
            qTmp = reshape(q, NyPrev, Nx);
            qPrevTmp = reshape(qPrev, NyPrev, Nx);
            
            qBordersDiffY = qTmp(My+1, :) - qTmp(My, :);
            subplot(2,1,2)
            plot(qBordersDiffY);
            ylim([-plotLim, plotLim])

            drawnow;
            if numFromBoundY == -1
                if mod(Ny, 2) == 0
                    qNext = reshape([qNextTmp(1:My-1, :); qNextTmp(My+1:end, :)], Ny*Nx, 1);
                    q = reshape([qTmp(1:My-1, :); qTmp(My+1:end, :)], Ny*Nx, 1);
                    qPrev = reshape([qPrevTmp(1:My-1, :); qPrevTmp(My+1:end, :)], Ny*Nx, 1);
                else
                    qNext = reshape([qNextTmp(1:My, :); qNextTmp(My+2:end, :)], Ny*Nx, 1);
                    q = reshape([qTmp(1:My, :); qTmp(My+2:end, :)], Ny*Nx, 1);
                    qPrev = reshape([qPrevTmp(1:My, :); qPrevTmp(My+2:end, :)], Ny*Nx, 1);
                end
            else
                qNext = reshape([qNextTmp(1:My-1, :); qNextTmp(My+1:end, :)], Ny*Nx, 1);
                q = reshape([qTmp(1:My-1, :); qTmp(My+1:end, :)], Ny*Nx, 1);
                qPrev = reshape([qPrevTmp(1:My-1, :); qPrevTmp(My+1:end, :)], Ny*Nx, 1);
            end
            disp("removed row at n = " + num2str(n))

        end

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
    
    if numFromBoundY == 1
        Dyy(My, My:(My+1)) = Dyy(My, My:(My+1)) + fliplr(ipY(2:3));
    else
        Dyy(My, My:(My+2)) = Dyy(My, My:(My+2)) + fliplr(ipY);
    end
    Dyy(My+1, (My-1):(My+1)) = Dyy(My+1, (My-1):(My+1)) + ipY;

    Dyy = Dyy / h^2;
    
    D = kron(speye(Nx), Dyy) + kron(Dxx, speye(Ny));
    DD = D * D;
    B = 2 * speye(Ny*Nx) - kappaSqVec(n) * k^2 * DD + 2 * sig1 * k * D;
    Amat = speye(Ny*Nx) * (1 + sig0 * k);
    C = (-1 + sig0 * k) * speye(Ny*Nx) - 2 * sig1 * k * D;
    
    %% Update equation
    qNext = Amat \ (B * q + C * qPrev);
    
    %% Displacement correction
    if dispCorr
        epsilon = 0; % Calculation is still defined for epsilon = 0 
        sigDC = 0.1; % Damping coefficient
        
        qNextTmp = reshape(qNext, Ny, Nx);
        qTmp = reshape(q, Ny, Nx);
        qPrevTmp = reshape(qPrev, Ny, Nx);
            
        % exclude connections
        startIndices = [1, Mx+2, My+2, 1];
        endIndices = [My-1, Nx, Ny, Mx-1];
        
        %% Calculate correction force
        
        % horizontal states (along vertical inner boundaries)
        etaPrevX1 = qPrevTmp(1:My-1, Mx+1) - qPrevTmp(1:My-1, Mx);
        etaPrevX2 = qPrevTmp(My+2:Ny, Mx+1) - qPrevTmp(My+2:Ny, Mx);
        etaNextX1 = qNextTmp(1:My-1, Mx+1) - qNextTmp(1:My-1, Mx);
        etaNextX2 = qNextTmp(My+2:Ny, Mx+1) - qNextTmp(My+2:Ny, Mx);

        % vertical states (along horizontal inner boundaries)
        etaPrevY1 = qPrevTmp(My+1, 1:Mx-1) - qPrevTmp(My, 1:Mx-1);
        etaPrevY2 = qPrevTmp(My+1, Mx+2:Nx) - qPrevTmp(My, Mx+2:Nx);
        etaNextY1 = qNextTmp(My+1, 1:Mx-1) - qNextTmp(My, 1:Mx-1);
        etaNextY2 = qNextTmp(My+1, Mx+2:Nx) - qNextTmp(My, Mx+2:Nx);

        rForce = (1 - sigDC / k) / (1 + sigDC / k);
        oOPX = (h^2 * (1 + sigDC / k) * (1-alfX)) / (2 * h^2 * (alfX + epsilon) + 2 * k^2 * (1 + sigDC / k) * (1-alfX));
        oOPY = (h^2 * (1 + sigDC / k) * (1-alfY)) / (2 * h^2 * (alfY + epsilon) + 2 * k^2 * (1 + sigDC / k) * (1-alfY));
        
        % Here, uNext and wNext are the 'intermediate' states of u and w (schemes without connection forces)
        FX1 = (etaNextX1 + rForce * etaPrevX1) * oOPX;
        FX2 = (etaNextX2 + rForce * etaPrevX2) * oOPX;
        FY1 = (etaNextY1 + rForce * etaPrevY1) * oOPY;
        FY2 = (etaNextY2 + rForce * etaPrevY2) * oOPY;
        
        %% Add displacement correction to inner boundaries
        qNextTmp(1:My-1, Mx) = qNextTmp(1:My-1, Mx) + k^2/h^2 * FX1;
        qNextTmp(1:My-1, Mx+1) = qNextTmp(1:My-1, Mx+1) - k^2/h^2 * FX1;

        qNextTmp(My+2:end, Mx) = qNextTmp(My+2:end, Mx) + k^2/h^2 * FX2;
        qNextTmp(My+2:end, Mx+1) = qNextTmp(My+2:end, Mx+1) - k^2/h^2 * FX2;
        
        qNextTmp(My, 1:Mx-1) = qNextTmp(My, 1:Mx-1) + k^2/h^2 * FY1;
        qNextTmp(My+1, 1:Mx-1) = qNextTmp(My+1, 1:Mx-1) - k^2/h^2 * FY1;
        
        qNextTmp(My, Mx+2:Nx) = qNextTmp(My, Mx+2:Nx) + k^2/h^2 * FY2;
        qNextTmp(My+1, Mx+2:Nx) = qNextTmp(My+1, Mx+2:Nx) - k^2/h^2 * FY2;
        
        betaX = (1-alfX) / (alfX + epsilon);
        betaY = (1-alfY) / (alfY + epsilon);
        pxPlus = betaX * (1 + sigDC / k) * k^2 / (2 * h^2);
        pyPlus = betaY * (1 + sigDC / k) * k^2 / (2 * h^2);
        pxMin = betaX * (1 - sigDC / k) * k^2 / (2 * h^2);
        pyMin = betaY * (1 - sigDC / k) * k^2 / (2 * h^2);
        etaPrevX1 = qPrevTmp(My, Mx+1) - qPrevTmp(My, Mx); % etaX1
        etaPrevY1 = qPrevTmp(My+1, Mx) - qPrevTmp(My, Mx); % etaY1
        etaPrevY2 = qPrevTmp(My+1, Mx+1) - qPrevTmp(My, Mx+1); % etaY2
        etaPrevX2 = qPrevTmp(My+1, Mx+1) - qPrevTmp(My+1, Mx); % etaX2

        pMat = [1 + pxPlus + pyPlus, -pxPlus, -pyPlus, 0;
                -pxPlus, 1 + pxPlus + pyPlus, 0, -pyPlus;
                -pyPlus, 0, 1 + pxPlus + pyPlus, -pxPlus;
                0, -pyPlus, -pxPlus, 1 + pxPlus + pyPlus];
        v2D = [qNextTmp(My, Mx)     + pxMin * etaPrevX1 + pyMin * etaPrevY1;
               qNextTmp(My, Mx+1)   - pxMin * etaPrevX1 + pyMin * etaPrevY2;
               qNextTmp(My+1, Mx)   + pxMin * etaPrevX1 - pyMin * etaPrevY1;
               qNextTmp(My+1, Mx+1) - pxMin * etaPrevX2 - pyMin * etaPrevY2];
        solut2D = pMat \ v2D;
        qNextTmp(My, Mx) = solut2D(1);
        qNextTmp(My, Mx+1) = solut2D(2);
        qNextTmp(My+1, Mx) = solut2D(3);
        qNextTmp(My+1, Mx+1) = solut2D(4);
        qNext = reshape(qNextTmp, Ny*Nx, 1);

    end

    
    out(n) = q((xInLoc-1 + 3) * Ny + (yInLoc + 3));

    
    if drawThings && mod(n, drawSpeed) == 0 && n > drawStart
        reshapedQ = reshape(q, Ny, Nx);
        qWithBoundaries = zeros(Ny + 2, Nx + 2);
        qWithBoundaries(2:end-1, 2:end-1) = reshapedQ;
        imagesc(qWithBoundaries);
        set(gca, 'Position', [0, 0, 1, 1], 'Clim', [-0.5, 0.5]);
        colormap gray;

        drawnow;
        
    end
    if n > lengthSound * nCounter / 100
        percentCounter = percentCounter + 1;
        nCounter = nCounter + 1;
        disp((percentCounter) + "% done")
    end
    
    % Update system states
    qPrev = q;
    q = qNext;
    
end
spectrogram (out,512,64,512, 44100, 'yaxis');
