function [ Rxx ] = NpCorrelationEstimate( x ,len)
%NpCorrelationEstimate devulve los valores de la autocorrelacion de la
%entrada
%   La funcion estima los valores de la autocorrelacion mediante el
% estimador no polarizado (1/N-k). sum{X(k)X(n+K)}, donde la sumatoria
% va desde el primer elemetno de X hasta el N-k esimo. len es el numero
% de valores ue se desean de la autocorrelacion.
N= size(x);
N = N(2);
Rxx =zeros(1,len+1);
for k = 1:(len+1)
    
    for n=1:N-k
        Rxx(k) = Rxx(k)+ x(n)*x(n+k-1);
    end
    Rxx(k) = Rxx(k) / (N-k);
end

