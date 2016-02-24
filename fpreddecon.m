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





