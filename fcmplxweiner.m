function F = fcmplxweiner(T_SigNoise, T_DesSigNoise, Nfilt)

Ninput = length(T_SigNoise);
Nd = Ninput + Nfilt - 1;

d = T_DesSigNoise;
d = [d, zeros(1, Nd-length(d))];

X = convmtx(T_SigNoise.', Nfilt);
R = ctranspose(X) * X;

P = real(R)
Q = imag(R)

g = ctranspose(X) * d.';

top = [P Q.'];
bottom = [Q P];
NE = [top;bottom];

g = g.';
g = ([g, zeros(1, length(NE) - length(g))]).';

f = NE\g;

Nf = length(f) / 2;
F = zeros(1, Nf);
for i = 1:Nf
    F(1,i) = f(i,1) + (f(Nf+i, 1)*sqrt(-1));
end

%Convolve the signal with the filter now

Y = conv(F, T_SigNoise);
ghf = (ctranspose(g(1:Nf, 1)))*F.';
delta = (ctranspose(d)).'*d.';
epsilon = 1.0 - ghf/delta

return;



