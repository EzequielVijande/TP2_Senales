function [x,s,sp,e,Gp,fs] = TP2_ej2(n,snr,audioPath)
% Abro audio, ploteo y escucho
[s,fs] = audioread(audioPath);
% Sumo ruido blanco
variance = 10^(-snr/10);
% Calculo varianza del ruido blanco
x = awgn(s,snr);
% Calculo de predicotres lineales optimos y ganancias
muestras = round(fs*20*10^(-3)); % cantidad de muestras a analizar en la cual la funcion es ergodica y estacionaria
cantH = floor(size(s)/muestras);
cantH = cantH(1);
coef = zeros(n,cantH);
N = cantH*muestras;
sp = zeros(1,N);
for i = 0.0:1.0:cantH-1
    xAux = x(i*muestras + 1:i*muestras + muestras);
    % calculo el estimador de la funcion autocorrelacion
    acf = zeros(1,muestras-1);
    for j = 0.0:1.0:muestras-1
        factor = (muestras-j);
        cumsum = 0;
        for k = 0.0:1.0:(muestras-j-1)
            cumsum = cumsum + xAux(k+1)*xAux(k+j+1);
        end
        acf(j+1) = cumsum/factor;
    end
    acf(1) = acf(1)-variance^2; % resto autocorrelacion ocasionada por el ruido
    matrix = zeros(n,n);
    vector = zeros(n,1);
    for j = 1.0:1.0:n % armo filas de la matriz y armo vector
        vector(j,1) = acf(j+1);
        for k = 1.0:1.0:n % armo columnas de la matriz
            matrix(j,k)= acf(abs(j-k)+1);
        end
    end
    h = matrix\vector; % equivalente a inv(matrix)*vector
    % agrego vector h a coef
    for l = 1.0:1.0:n
        coef(l,i+1) = h(l);
    end 
    ht = [0;h];
    spAux = conv(xAux,ht,'same');
    p = 1;
    for m = (i*muestras+1):1.0:(i*muestras + muestras)
        sp(m) = spAux(p);
        p = p + 1;
    end
end
sp = transpose(sp);
x = x(1:N);
s = s((1+ceil(n/2)):N);
s = [s;zeros(ceil(n/2),1)];
% Calculo error
e = zeros(1,N);
for i = 1.0:1.0:N % señal cumsum
    e(i) = s(i)-sp(i); 
end
% Calculo ganancia
sCumsum = 0;
eCumsum = 0;
for i = 1.0:1.0:N % señal cumsum
    sCumsum = sCumsum + s(i)^2;
    eCumsum = eCumsum + e(i)^2; 
end
Gp = sCumsum/eCumsum;
end