fs = 44100;
k = 1/fs;
c = 100;

Lx = linspace(0.25, 0.5, 1000);
Ly = 0.01;                 % Length in y direction [m]

h = sqrt(2) * c * k; 
NxMax = floor(max(Lx)/h);
NxMin = floor(min(Lx)/h);

freqsSave = zeros(length(Lx), NxMax);
    
for i = 1:length(Lx)
    h = sqrt(2) * c * k; 
    NxFrac = (Lx(i)/h);
    Nx = floor(NxFrac);
    NxPrev = Nx;
    alf = NxFrac - Nx;
%     alf = 0; % test
    Mx = ceil(Nx * 0.5);
    Mwx = floor(Nx * 0.5);
    Ny = floor(Ly/h);       % Number of intervals between grid points
    h = max(Lx(i)/Nx, Ly/Ny);  % Recalculation of grid spacing based on integer N

    lambdaSq = c^2 * k^2 / h^2;


    % matrices
    DxxU = toeplitz([-2, 1, zeros(1, Mx-2)]);
    DxxW = toeplitz([-2, 1, zeros(1, Mwx-2)]);
    Dxx = zeros(Nx);
    Dxx(1:Mx, 1:Mx) = DxxU;
    Dxx(Mx+1:end, Mx+1:end) = DxxW;
    ip = [(alf - 1) / (alf + 1), 1, -(alf - 1) / (alf + 1)];
    Dxx(Mx, Mx:Mx+2) = Dxx(Mx, Mx:Mx+2) + ip;
    Dxx(Mx+1, Mx-1:Mx+1) = Dxx(Mx+1, Mx-1:Mx+1) + fliplr(ip);

    
    Dyy = toeplitz([-2, 1, zeros(1, Ny-3)]);
    
    D = 1/h^2 * (kron(speye(Nx), Dyy) + kron(Dxx, speye(Ny-1)));
    B = 2 * speye(Nx*(Ny-1)) + c^2 * k^2 * D;

    freqs = 1/(pi * k) * asin (c * k / 2 * sqrt(-eig(full(D))));
    %     freqsUse = sort(freqs);
    %     freqsUse = freqsUse(freqsUse < 22050);
    freqsSave(i, 1:length(freqs)) = sort(freqs);
    hold off;
    scatter(1:length(freqs), sort(freqs));
    hold on;
    scatter(NxMin*(Ny-1), freqsSave(i, NxMin*(Ny-1)), 80, 'r', 'filled');
    eps = Lx(i)/Ly;
    NxMat = repmat(1:Nx, Ny, 1);
    NyMat = repmat((1:Ny)', 1, Nx);
    modes = c / (2 * sqrt(Lx(i) * Ly)) * sqrt(NxMat.^2/eps + eps * NyMat.^2);
    modes = sort(reshape(modes, Nx * Ny, 1));
    plot(modes);
    drawnow; 
    
   disp((floor(i/length(Lx) * 1000) / 10) + "% done");
end
    hold off

for i = 1:NxMin*(Ny-1)
    plot(freqsSave(:, i));
    hold on;
    drawnow;
%     pause(0.3)
end

