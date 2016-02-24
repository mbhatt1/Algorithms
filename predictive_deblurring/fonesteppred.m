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
for i = 1: length(phi),
    normphi(k) = phi(i)/phi(1);
    k  = k + 1;
end 

k = 1
for i = 1 : nf,
    covarahead(k) = normphi(i+PreTime);
    k = k + 1;
end

phi_m = toeplitz(normphi(1, 1:nf));
a_t = phi_m\covarahead';

for i = 1 : nf
    if (i == 1)
        peo(1) = 1;
    else
        peo(i) = 1 .* a_t(i-1);
    end
end


%Filter Application and Obtaining the Predicted Signal

y_t = conv(x, peo);
dt = x(2:nx);
D = dt * ctranspose(dt);
g_prime = covarahead ./ D;

for i= 1 : nf, 
    y1 = peo(i)
    fg_prime(i) = y1 .* g_prime(i);
end

ErrEnergy = 1 - sum(fg_prime)

PerfParm = 1-ErrEnergy

return;


