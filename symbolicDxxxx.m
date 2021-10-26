syms alf Dxx
N = 30;
numFromBound = -1;
M = ceil (0.5 * N);
Mw = floor (0.5 * N);

Dxx = sym(zeros(M+Mw));

Dxx(1:M, 1:M) = toeplitz([-2, 1, zeros(1, M-2)]);
Dxx(M+1:end, M+1:end) = toeplitz([-2, 1, zeros(1, Mw-2)]);

%% Calculate (quadratic) interpolator and virtual grid points
ip = [-(alf - 1) / (alf + 1), 1, (alf - 1) / (alf + 1)];
testMat = sym(zeros(4));
testMat(2,2:4) = fliplr(ip);
testMat(3,1:3) = ip;
testMat * testMat
if numFromBound == 1
    Dxx(M, M:(M+1)) = Dxx(M, M:(M+1)) + fliplr(ip(2:end));
    Dxx(M+1, (M-1):(M+1)) = Dxx(M+1, (M-1):(M+1)) + ip;
else
    Dxx(M, M:(M+2)) = Dxx(M, M:(M+2)) + fliplr(ip);
    Dxx(M+1, (M-1):(M+1)) = Dxx(M+1, (M-1):(M+1)) + ip;

end
Dxxxx = Dxx * Dxx;