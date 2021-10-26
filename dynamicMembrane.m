close all;
clear all;

drawThings = false;
energyCalc = true;
plotPropagation = false;
if plotPropagation 
    figure('Position', [440 585 815 213])
end
drawSpeed = 5;
drawStart = 0;

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs * 1;   % Length of the simulation (1 second) [samples]             

c = 100;

Lx = linspace(0.5, 2, lengthSound);

% Lx = 0.5;
Ly = 0.05;                 % Length in y direction [m]

h = sqrt(2) * c * k; 
Nx = floor(Lx(1)/h);
NxPrev = Nx;

Mx = ceil(Nx * 0.5);
Mwx = floor(Nx * 0.5);
Ny = floor(Ly/h);       % Number of intervals between grid points
h = min(Lx(1)/Nx, Ly/Ny);  % Recalculation of grid spacing based on integer N
% changing c seems to require max(Lx/Nx, Ly/Ny) to work smoothly..
lambdaSq = c^2 * k^2 / h^2;

%% Initialise state vectors (one more grid point than the number of intervals)
uNext = zeros(Ny+1, Mx+1); 
u = zeros(Ny+1, Mx+1); 

wNext = zeros(Ny+1, Mwx+1); 
w = zeros(Ny+1, Mwx+1); 

% don't subtract from Nx as there is one point for overlap
qNext = zeros((Ny-1) * (Nx), 1);
q = zeros((Ny-1) * (Nx), 1);

%% Initial conditions (raised cosine)
% halfWidth = floor(min(Nx, Ny) / 10);
halfWidth = 1;
width = 2 * halfWidth + 1;
xLoc = 10;
yLoc = floor((Ny+1) / 2);
xRange = xLoc-halfWidth : xLoc+halfWidth;
yRange = yLoc-halfWidth : yLoc+halfWidth;

rcMatHalf = zeros(Ny+1, Mx+1);
rcMatHalf(yRange, xRange) = 1 * hann(width) * hann(width)';
u = rcMatHalf; % initialise current state  \

rcMat = zeros(Ny-1, Nx);
rcMat(yRange-1, xRange-1) = 1 * hann(width) * hann(width)';

q = reshape(rcMat, (Ny-1)*(Nx), 1);
outlocQ = find(q == max(q));

% Set initial velocity to zero
uPrev = u;
wPrev = w;
qPrev = q;

out = zeros(lengthSound, 1);


nCounter = 0;
percentCounter = 0;
%% Simulation loop
for n = 1:lengthSound
    
    NxFrac = Lx(n)/h;
    Nx = floor(NxFrac);
    alf = NxFrac - Nx;
    if Nx ~= NxPrev
        if Nx > NxPrev
            
            bordersNext = [uNext(:, Mx:Mx+1), wNext(:, 1:2)]; % unnecessary
            borders = [u(:, Mx:Mx+1), w(:, 1:2)];
            bordersPrev = [uPrev(:, Mx:Mx+1), wPrev(:, 1:2)];

            
            qNextTmp = reshape(qNext, Ny-1, NxPrev);
            qTmp = reshape(q, Ny-1, NxPrev);
            qPrevTmp = reshape(qPrev, Ny-1, NxPrev);

            qBordersNext = [qNextTmp(:, Mx-1:Mx), qNextTmp(:, Mx+1:Mx+2)]; % unnecessary
            qBorders = [qTmp(:, Mx-1:Mx), qTmp(:, Mx+1:Mx+2)];
            qBordersPrev = [qPrevTmp(:, Mx-1:Mx), qPrevTmp(:, Mx+1:Mx+2)];

            
            cubicIp = [alf * (alf + 1) / -((alf + 2) * (alf + 3)); ...
                2 * alf / (alf + 2); ...
                2 / (alf + 2); ...
                2 * alf / -((alf + 3) * (alf + 2))]';

            if mod(Nx, 2) == 1 
                uNext = [uNext, bordersNext * cubicIp'];
                u = [u, borders * cubicIp'];
                uPrev = [uPrev, bordersPrev * cubicIp'];
%                 disp("Added to u with alf = " + num2str(alf))
                qNext = reshape([qNextTmp(:, 1:Mx), qBordersNext * cubicIp', qNextTmp(:, Mx+1:end)], (Ny-1) * Nx, 1);
                q = reshape([qTmp(:, 1:Mx), qBorders * cubicIp', qTmp(:, Mx+1:end)], (Ny-1) * Nx, 1);
                qPrev = reshape([qPrevTmp(:, 1:Mx), qBordersPrev * cubicIp', qPrevTmp(:, Mx+1:end)], (Ny-1) * Nx, 1);


            else
                wNext = [bordersNext * fliplr(cubicIp)', wNext];
                w = [borders * fliplr(cubicIp)', w];
                wPrev = [bordersPrev * fliplr(cubicIp)', wPrev];
%                 disp("Added to w with alf = " + num2str(alf))
                qNext = reshape([qNextTmp(:, 1:Mx), qBordersNext * fliplr(cubicIp)', qNextTmp(:, Mx+1:end)], (Ny-1) * Nx, 1);
                q = reshape([qTmp(:, 1:Mx), qBorders * fliplr(cubicIp)', qTmp(:, Mx+1:end)], (Ny-1) * Nx, 1);
                qPrev = reshape([qPrevTmp(:, 1:Mx), qBordersPrev * fliplr(cubicIp)', qPrevTmp(:, Mx+1:end)], (Ny-1) * Nx, 1);

            end
            

            
            

        else
            %%%
        end
        Mx = ceil(Nx * 0.5);
        Mwx = floor(Nx * 0.5); % CHANGING Mx AND Mwx HERE!!
    end
    NxPrev = Nx;
    
    ip = [(alf - 1) / (alf + 1), 1, -(alf - 1) / (alf + 1)];

    Dxx = zeros(Nx);
    Dxx(1:Mx, 1:Mx) = toeplitz([-2, 1, zeros(1, Mx-2)]);
    Dxx(Mx+1:end, Mx+1:end) = toeplitz([-2, 1, zeros(1, Mwx-2)]);
    
    Dxx(Mx, Mx:(Mx+2)) = Dxx(Mx, Mx:(Mx+2)) + ip;
    Dxx(Mx+1, (Mx-1):(Mx+1)) = Dxx(Mx+1, (Mx-1):(Mx+1)) + fliplr(ip);

    Dxx = Dxx / h^2;
    
    Dyy = toeplitz([-2, 1, zeros(1, Ny-3)]);
    Dyy = Dyy / h^2;
    
    D = kron(speye(Nx), Dyy) + kron(Dxx, speye(Ny-1));

    B = 2 * speye((Ny-1) * Nx) + c^2 * k^2 * D;
    
    borders = [u(2:Ny, Mx:Mx+1), w(2:Ny, 1:2)];
    uInterp = borders * [0, ip]';
    wInterp = borders * fliplr([0, ip])';
    
    %% Update equation
    uNext(2:Ny, 2:Mx) = 2 * u(2:Ny, 2:Mx) - uPrev(2:Ny, 2:Mx) ...
        + lambdaSq * (u(3:Ny+1, 2:Mx) + u(1:Ny-1, 2:Mx) + u(2:Ny, 3:Mx+1) + u(2:Ny, 1:Mx-1) - 4 * u(2:Ny, 2:Mx)); 
       
    uNext(2:Ny, Mx+1) = 2 * u(2:Ny, Mx+1) - uPrev(2:Ny, Mx+1) ...
        + lambdaSq * (u(3:Ny+1, Mx+1) + u(1:Ny-1, Mx+1) + uInterp + u(2:Ny, Mx) - 4 * u(2:Ny, Mx+1)); 

    
    wNext(2:Ny, 2:Mwx) = 2 * w(2:Ny, 2:Mwx) - wPrev(2:Ny, 2:Mwx) ...
        + lambdaSq * (w(3:Ny+1, 2:Mwx) + w(1:Ny-1, 2:Mwx) + w(2:Ny, 3:Mwx+1) + w(2:Ny, 1:Mwx-1) - 4 * w(2:Ny, 2:Mwx)); 

    wNext(2:Ny, 1) = 2 * w(2:Ny, 1) - wPrev(2:Ny, 1) ...
        + lambdaSq * (w(3:Ny+1, 1) + w(1:Ny-1, 1) + w(2:Ny, 2) + wInterp - 4 * w(2:Ny, 1)); 

    qNext = B * q - qPrev;
    
    out(n) = u(yLoc, xLoc);
    outQ(n) = q(outlocQ);
    
    if drawThings && mod(n, drawSpeed) == 0 && n > drawStart
        subplot(321)
        imagesc(u);
        set(gca, 'Position', [0, 0.66, 0.5, 0.33], 'CLim', [-0.05, 0.05]);
        colormap gray;
        subplot(322)
        imagesc(w);
        set(gca, 'Position', [0.5, 0.66, 0.5, 0.33], 'CLim', [-0.05, 0.05]);
        colormap gray;
        
        subplot(3,2, [3,4])
        reshapedQ = reshape(q, Ny-1, Nx);
        qWithBoundaries = zeros(Ny + 1, Nx + 2);
        qWithBoundaries(2:end-1, 2:end-1) = reshapedQ;
        imagesc(qWithBoundaries);
        set(gca, 'Position', [0, 0.33, 1, 0.33], 'CLim', [-0.05, 0.05]);
        colormap gray;

        subplot(3,2, [5,6])
        imagesc([u, w] - qWithBoundaries);
        set(gca, 'Position', [0, 0, 1, 0.33], 'CLim', [-0.05, 0.05]);
        colormap gray;
        drawnow;
        
    end
    if n > lengthSound * nCounter / 100
        percentCounter = percentCounter + 1;
        nCounter = nCounter + 1;
        disp((percentCounter) + "% done")
    end
        % Update system states
    uPrev = u;
    u = uNext;
    wPrev = w;
    w = wNext;

    qPrev = q;
    q = qNext;
end

