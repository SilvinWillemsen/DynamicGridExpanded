syms alf Dxx A
N = 20;
numFromBound = 1;
if numFromBound == -1
    M = ceil (0.5 * N);
    Mw = floor (0.5 * N);
else
    M = N-1;
    Mw = 1;
end
Dxp = sym(zeros(M+Mw));
eu = ones(M+Mw, 1);
% Dxp(1:end, 1:end) = spdiags([-eu, eu], 0:1, M+Mw, M+Mw);
% Dxp(M, M) = Dxp(M, M) + A;
% Dxp(M, M+2) = - ((alf - 1) / (alf+1));


Dxx = sym(zeros(M+Mw));

Dxx(1:M, 1:M) = spdiags([eu -2*eu, eu], -1:1, M, M);
Dxx(M+1:end, M+1:end) = spdiags([eu -2*eu, eu], -1:1, Mw, Mw);

%% Calculate (quadratic) interpolator and virtual grid points
ip = [-A, 1, A];
testMat = sym(zeros(4));
testMat(2,2:4) = fliplr(ip);
testMat(3,1:3) = ip;
testMat * testMat;
if numFromBound == 1
    Dxx(M, M:(M+1)) = Dxx(M, M:(M+1)) + fliplr(ip(2:end));
    Dxx(M+1, (M-1):(M+1)) = Dxx(M+1, (M-1):(M+1)) + ip;
else
    Dxx(M, M:(  M+2)) = Dxx(M, M:(M+2)) + fliplr(ip);
    Dxx(M+1, (M-1):(M+1)) = Dxx(M+1, (M-1):(M+1)) + ip;
    %%%
%     Dxx(M+2, M) = ip(1);
%     Dxx(M-1, M+1) = ip(1);
    %%%


end
Dxxxx = Dxx * Dxx;

kappaSq = 1;
k= 1/44100;
h = sqrt(2 * sqrt(kappaSq) * k);
B = 2 * eye(size(Dxx)) + Dxx
Q = [B, -eye(size(Dxx));eye(size(Dxx)), zeros(size(Dxx))]
answer = eig(Q);
% shouldBeOnes = abs(subs(answer, 0:0.01:1))
%%
close all
for alfSubs = 0:0.01:1
    eigs = subs(answer, alfSubs);
    plot(exp(1j*2*pi*(0:0.01:1)))
    hold on
    scatter(real(eigs), imag(eigs))
    hold off;
    xlim([-1.5, 1.5])
    ylim([-1.5, 1.5])
    drawnow;
end