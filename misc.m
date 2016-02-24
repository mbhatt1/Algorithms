%Function fonesteppred

function y_t = fonesteppred(x, nf)

PreTime = 1; 
nx = length(x);

[sigcovar, lags] = xcorr(x, nf+PreTime, 'none');

k = 1;

for i = 1: length(sigcovar), 
    if (lags(i) >= 0)
        phi(k) = sigcovar(i)
        k = k + 1;
    end
end


k = 1;

for i = 1 : nf,
    covarahead(k) = normphi(i+PreTime);
    k = k + 1
end

phi_m = toeplitz(normphi(1, 1:nf));
a_t = phi_m\covarahead';

for i = 1 : nf
    if (i == 1)
        peo(1) = 1;
    else
        peo(1) = 1 .* a_t(i-1);
    end
end


%Filter Application and Obtaining the Predicted Signal

y_t = conv(x, peo);
dt = x(2:nx);
D = dt * ctranspose(dt);
g_prime = covarahead ./ D;

for i= 1 : nf, 
    fg_prime(i) = peo(i) .* g_prime(i);
end

ErrEnergy = 1 - sum(fg_prime)

PerfParm = 1-ErrEnergy

return;


function [y_order, x_order] = forder(XS, Nx)
y_neg = [XS( (Nx/2)+1:Nx)];
y_pos = [Xs(1:(Nx/2))];
clear y_order;
y_order = [y_neg, y_pos];
x_neg = -Nx/2+1:1:0;
x_pos = 1:1:Nx/2;
clear x_order;
x_order = [x_neg, x_pos];
return; 



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



function fpreddecon(signoise, signal, Nx, Nfilt)

T_SigNoise = fft(signoise);

PS = T_SigNoise .* conj(T_SigNoise);
T_Mod_SigNoise = T_SigNoise ./ sqrt(PS) .* -T_SigNoise;

mod_signoise = ifft(T_Mod_SigNoise, Nx);

pred_y = fonesteppred(signoise, Nfilt);

mod_y = fonesteppred(mod_signoise, Nfilt);

T_Pred_Y = fft(pred_y, Nx);

subplot (211), plot (1:Nx, signal(1:Nx), 'r:', 1:Nx, pred_y(1:Nx), 'b-');
xlabel ('(a) Time Domain, t')
ylabel ('Amplitude');


%Fig 2

figure;
subplot(311), plot(1:Nx, signoise(1:Nx), 'r:', 1:Nx, pred_y(1:Nx), 'b-')
title('Time Domain Conventional Prediction Deconvolution')
xlabel('a Time Domain, t')
ylabel('amplitude')

[yy,yy_xaxis] = forder(T_SigNoise, Nx);
[py, py_axis] = forder(T_Pred_Y, Nx);

subplot(312),plot(yy_xaxis, real(yy/Nx), 'r:', py_axis, real(py/Nx), 'b-')

return;





Nx = 128;
amplitude = 50.0;
for i = 1:Nx
    signal(i) = amplitude * sin(2*pi*i*.14)
    noise(i) = (randn(1,1)-.5)*25.0;
    data1(i) = signal(i)+noise(i);
end
fpreddecon(data1, signal, Nx, 128);
    