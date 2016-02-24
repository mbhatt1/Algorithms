for i = 1 : 128
    signal(i) = 50 * sin(2*22/7*.4* i)
    noise(i) = (randn(1,1)- .5)*.25
    data(i) = signal(i) + noise(i)
end
f = fpreddecon(data, signal, 128, 128)
