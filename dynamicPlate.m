close all;
clear all;

drawThings = false;
energyCalc = true;
plotPropagation = false;
if plotPropagation 
    figure('Position', [440 585 815 213])
end
drawSpeed = 1;
drawStart = 0;

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs * 5;   % Length of the simulation (1 second) [samples]             

rho = 7850;
H = 0.05;
E = 2e11;
nu = 0.3;
Dvar = E * H^3 / (12 * (1 - nu^2));
kappaSq = Dvar / (rho * H);

sig0 = 1;
sig1 = 0.05;

Lx = linspace(1, 5, lengthSound);

% Lx = 0.5;
Ly = linspace(0.5, 2, lengthSound);

h = 2 * sqrt(k * (sig1 + sqrt(sig1^2 + kappaSq))); 

Nx = floor(Lx(1)/h);
NxPrev = Nx;

Mx = ceil(Nx * 0.5);
Mwx = floor(Nx * 0.5);

Ny = floor(Ly(1)/h);
NyPrev = Ny;

My = ceil(Ny * 0.5);
Mwy = floor(Ny * 0.5);

h = min(Lx(1)/Nx, Ly(1)/Ny);  % Recalculation of grid spacing based on integer N

%% Initialise state vectors (one more grid point than the number of intervals)

% don't subtract from Nx as there is one point for overlap
qNext = zeros(Ny * Nx, 1);
q = zeros(Ny * Nx, 1);

%% Initial conditions (raised cosine)
% halfWidth = floor(min(Nx, Ny) / 10);
halfWidth = 1;
width = 2 * halfWidth + 1;
xInLoc = 5;
yInLoc = 5;

% Implemented to be "number-of-point from right / bottom boundary
xOutLoc = 5;
yOutLoc = 5; 

% Set initial velocity to zero
qPrev = q;

out = zeros(lengthSound, 1);

nCounter = 0;
percentCounter = 0;
%% Simulation loop
for n = 1:lengthSound
    if mod(n, fs/2) == 1
        q((xInLoc-1) * Ny + yInLoc) = 1;
    end
    
    NxFrac = Lx(n)/h;
    Nx = floor(NxFrac);
    alfX = NxFrac - Nx;

    if Nx ~= NxPrev
        if Nx > NxPrev
            
            qNextTmp = reshape(qNext, Ny, NxPrev);
            qTmp = reshape(q, Ny, NxPrev);
            qPrevTmp = reshape(qPrev, Ny, NxPrev);

            qBordersNextX = [qNextTmp(:, Mx-1:Mx), qNextTmp(:, Mx+1:Mx+2)]; % unnecessary
            qBordersX = [qTmp(:, Mx-1:Mx), qTmp(:, Mx+1:Mx+2)];
            qBordersPrevX = [qPrevTmp(:, Mx-1:Mx), qPrevTmp(:, Mx+1:Mx+2)];

            
            cubicIpX = [alfX * (alfX + 1) / -((alfX + 2) * (alfX + 3)); ...
                2 * alfX / (alfX + 2); ...
                2 / (alfX + 2); ...
                2 * alfX / -((alfX + 3) * (alfX + 2))]';

            if mod(Nx, 2) == 1 
                qNext = reshape([qNextTmp(:, 1:Mx), qBordersNextX * cubicIpX', qNextTmp(:, Mx+1:end)], Ny*Nx, 1);
                q = reshape([qTmp(:, 1:Mx), qBordersX * cubicIpX', qTmp(:, Mx+1:end)], Ny*Nx, 1);
                qPrev = reshape([qPrevTmp(:, 1:Mx), qBordersPrevX * cubicIpX', qPrevTmp(:, Mx+1:end)], Ny*Nx, 1);
            else
                qNext = reshape([qNextTmp(:, 1:Mx), qBordersNextX * fliplr(cubicIpX)', qNextTmp(:, Mx+1:end)], Ny*Nx, 1);
                q = reshape([qTmp(:, 1:Mx), qBordersX * fliplr(cubicIpX)', qTmp(:, Mx+1:end)], Ny*Nx, 1);
                qPrev = reshape([qPrevTmp(:, 1:Mx), qBordersPrevX * fliplr(cubicIpX)', qPrevTmp(:, Mx+1:end)], Ny*Nx, 1);

            end
   
        else
            %%%
        end
        disp("added column at n = " + num2str(n))
        Mx = ceil(Nx * 0.5);
        Mwx = floor(Nx * 0.5); % CHANGING Mx AND Mwx HERE!!
    end
    NxPrev = Nx;
    
    
    NyFrac = Ly(n)/h;
    Ny = floor(NyFrac);
    alfY = NyFrac - Ny;
%     if (alfX > 0.99 && ~drawThings) || (alfY > 0.99 && ~drawThings)
%         drawThings = true;
%     end
%     if (alfX > 0.01 && alfX < 0.99 && drawThings) && (alfY > 0.01 && alfY < 0.99 && drawThings) 
%         drawThings = false;
%     end
%     if alfY > 0.99 && ~drawThings
%         drawThings = true;
%     end
%     if alfY > 0.01 && alfY < 0.99 && drawThings 
%         drawThings = false;
%     end

%     if alfY > 0.95 && ~drawThings
%         drawThings = true;
%     end
%     if alfY > 0.05 && alfY < 0.95 && drawThings
%         drawThings = false;
%     end

    if Ny ~= NyPrev
        if Ny > NyPrev
            
            qNextTmp = reshape(qNext, NyPrev, Nx);
            qTmp = reshape(q, NyPrev, Nx);
            qPrevTmp = reshape(qPrev, NyPrev, Nx);

            qBordersNextY = [qNextTmp(My-1:My, :); qNextTmp(My+1:My+2, :)]; % unnecessary
            qBordersY = [qTmp(My-1:My, :); qTmp(My+1:My+2, :)];
            qBordersPrevY = [qPrevTmp(My-1:My, :); qPrevTmp(My+1:My+2, :)];

            
            cubicIpY = [alfY * (alfY + 1) / -((alfY + 2) * (alfY + 3)); ...
                2 * alfY / (alfY + 2); ...
                2 / (alfY + 2); ...
                2 * alfY / -((alfY + 3) * (alfY + 2))]';

            if mod(Ny, 2) == 1 
                qNext = reshape([qNextTmp(1:My, :); cubicIpY * qBordersNextY; qNextTmp(My+1:end, :)], Ny*Nx, 1);
                q = reshape([qTmp(1:My, :); cubicIpY * qBordersY; qTmp(My+1:end, :)], Ny*Nx, 1);
                qPrev = reshape([qPrevTmp(1:My, :); cubicIpY * qBordersPrevY; qPrevTmp(My+1:end, :)], Ny*Nx, 1);
            else
                qNext = reshape([qNextTmp(1:My, :); fliplr(cubicIpY) * qBordersNextY; qNextTmp(My+1:end, :)], Ny*Nx, 1);
                q = reshape([qTmp(1:My, :); fliplr(cubicIpY) * qBordersY; qTmp(My+1:end, :)], Ny*Nx, 1);
                qPrev = reshape([qPrevTmp(1:My, :); fliplr(cubicIpY) * qBordersPrevY; qPrevTmp(My+1:end, :)], Ny*Nx, 1);

            end
   
        else
            %%%
        end
        disp("added row at n = " + num2str(n))

        My = ceil(Ny * 0.5);
        Mwy = floor(Ny * 0.5); % CHANGING Mx AND Mwx HERE!!
    end
    NyPrev = Ny;
    
%     alfX = 1;
    ipX = [(alfX - 1) / (alfX + 1), 1, -(alfX - 1) / (alfX + 1)];

    Dxx = zeros(Nx);
    Dxx(1:Mx, 1:Mx) = toeplitz([-2, 1, zeros(1, Mx-2)]);
    Dxx(Mx+1:end, Mx+1:end) = toeplitz([-2, 1, zeros(1, Mwx-2)]);
    
    Dxx(Mx, Mx:(Mx+2)) = Dxx(Mx, Mx:(Mx+2)) + ipX;
    Dxx(Mx+1, (Mx-1):(Mx+1)) = Dxx(Mx+1, (Mx-1):(Mx+1)) + fliplr(ipX);

    Dxx = Dxx / h^2;
    
    ipY = [(alfY - 1) / (alfY + 1), 1, -(alfY - 1) / (alfY + 1)];

    Dyy = zeros(Ny);
    Dyy(1:My, 1:My) = toeplitz([-2, 1, zeros(1, My-2)]);
    Dyy(My+1:end, My+1:end) = toeplitz([-2, 1, zeros(1, Mwy-2)]);

%     Dyy(My:(My+2), My) = Dyy(My:(My+2), My) + ipY';
%     Dyy((My-1):(My+1), My+1) = Dyy((My-1):(My+1), My+1) + fliplr(ipY)';
    Dyy(My, My:(My+2)) = Dyy(My, My:(My+2)) + ipY;
    Dyy(My+1, (My-1):(My+1)) = Dyy(My+1, (My-1):(My+1)) + fliplr(ipY);

    Dyy = Dyy / h^2;
    
    D = kron(speye(Nx), Dyy) + kron(Dxx, speye(Ny));
    DD = D * D;
    B = 2 * speye(Ny*Nx) - kappaSq * k^2 * DD + 2 * sig1 * k * D;
    Amat = speye(Ny*Nx) * (1 + sig0 * k);
    C = (-1 + sig0 * k) * speye(Ny*Nx) - 2 * sig1 * k * D;
    %% Update equation
    qNext = Amat \ (B * q + C * qPrev);
    
    out(n) = q((Nx - xOutLoc-1) * Ny + (Ny - yOutLoc));

    
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

