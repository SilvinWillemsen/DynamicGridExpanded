%{
    An implementation of a dynamic grid for finite-difference schemes 
    accompanying the following publication:
        S. Willemsen, S. Bilbao, M. Ducceschi and S. Serafin, "Dynamic 
        Grids for Finite-Difference Schemes in Musical Instrument 
        Simulations", Proceedings of the 23rd Int. Conf. on Digital Audio 
        Effects (DAFx), 2021.
%}

clear all;
close all;
clc;

%% Draw settings
drawThings = false;
drawStart = 0;
drawSpeed = 100;

fs = 44100;             % Sample rate (in Hz)
k = 1/fs;               % Time step (in s)
lengthSound = fs * 10;   % Length of the simulation (in samples)

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
%% Initialise wave speeds for the entire simulation
L = 1;                  % System length (does not really matter)
rho = 7850;
r = 0.0005;
A = pi * r^2; 
T0 = 300;
E = 2e11;
sig0 = 1;

halfWidth = 5;
initState = hann(2*halfWidth + 1);
T = T0 + initState' * initState;
c = sqrt(T / (rho * A));
h = c * k;
Nfrac = L / h;
N = floor(Nfrac);
NPrev = N;
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

% Excite left system using u_1 = 1 (excites all modes)
startLoc = M/2 - halfWidth;
endLoc = M/2 + halfWidth;
u(startLoc:endLoc) = initState; 
uPrev = u;

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

    % Retrieve new c
    T = T0 + u' * u;
    c = sqrt(T / (rho * A));
    
    % recalculate gridspacing and (fractional) number of intervals
    h = c * k;
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
                    
                end
                
            % Otherwise, add grid points to u using interpolation
            elseif numFromBound == 1
                uNext = [uNext; 0];
                u = [u; cubicIp(1:3) * [u(end-1:end); w(1)]];
                uPrev = [uPrev; cubicIp(1:3) * [uPrev(end-1:end); wPrev(1)]];

            else
                uNext = [uNext; 0];
                u = [u; cubicIp * [u(end-1:end); w(1:2)]];
                uPrev = [uPrev; cubicIp * [uPrev(end-1:end); wPrev(1:2)]];

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
                    
                else 
                    % if N is odd, remove from right system
                    wNext = wNext(2:end);
                    w = w(2:end);
                    wPrev = wPrev(2:end);
                    
                end
                
            % Otherwise, remove grid point from u
            else
                uNext = uNext(1:end-1);
                u = u(1:end-1);
                uPrev = uPrev(1:end-1);
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
    
    %% Calculate (quadratic) interpolator and virtual grid points
    ip = [-(alf - 1) / (alf + 1), 1, (alf - 1) / (alf + 1)];

    % virtualGridPoints = [u_{M+1}^n; w_{-1}^n]
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

    %% Displacement correction
    if dispCorr
        epsilon = 0; % Calculation is still defined for epsilon = 0 
        sig0 = 1; % Damping coefficient

        %% Calculate correction force
        etaPrev = wPrev(1) - uPrev(end);
        rForce = (1 - sig0 / k) / (1 + sig0 / k);
        oOP = (h * (1 + sig0 / k) * (1-alf)) / (2 * h * (alf + epsilon) + 2 * k^2 * (1 + sig0 / k) * (1-alf));
        
        % Here, uNext and wNext are the 'intermediate' states of u and w (schemes without connection forces)
        F = (wNext(1) - uNext(end) + rForce * etaPrev) * oOP;
        
        %% Add displacement correction to inner boundaries
        uNext(end) = uNext(end) + k^2/h * F;
        wNext(1) = wNext(1) - k^2/h * F;

    end

    %% save output
    out(n) = u(1);

    %% draw stuff
    if n > drawStart && drawThings && mod(n, drawSpeed) == 0

        % Grid point locations
        hLocsLeft = (0:(length(u)))*h;
        hLocsRight = (fliplr(L - (0:(length(w)))*h));
        
        hold off;
        
        % Plot left system (with left outer boundary)
        plot(hLocsLeft, [0;u], 'LineWidth' ,2, 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
        hold on;
        
        % Plot right system (with right outer boundary)
        plot(hLocsRight, [w; 0], 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
        
        % Settings
        xlim([0, L])
%         ylim([-1, 1])
        grid on;
        xlabel("$x$ (m)", 'interpreter', 'latex')
        ylabel("displacement (m)")
        title("$\alpha = " + num2str(round(alf*100)/100) + "$", 'interpreter', 'latex')
        legend(["$u$", "$w$"], 'interpreter', 'latex')
        set(gca, "Fontsize", 16, 'Linewidth', 2)
        
        drawnow;
    end
    
    %% Update states
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;
    
    NPrev = N;
    
end

figure;
spectrogram (out,512,64,512, 44100, 'yaxis');