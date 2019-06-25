function [ Rxx ] = PolCorrelationEstimate( x,len )
%PolCorrelationEstimate devulve los valores de la autocorrelacion de la
%entrada
%   La funcion estima los valores de la autocorrelacion mediante el
% estimador polarizado (1/N). sum{X(k)X(n+K)}, donde la sumatoria
% va desde el primer elemento de X hasta el N-k esimo. len es el numero
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
end

