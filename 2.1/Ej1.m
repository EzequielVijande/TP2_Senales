x = load('Archivo_4.mat');
x= x.x;
len_x = size(x);
len_x = len_x(2);
len = 128;
k = [0:len,0:len];
Rxx_n = NpCorrelationEstimate(x,len);
rxx_n = Rxx_n ./ Rxx_n(1);
Rxx_p = PolCorrelationEstimate( x,len );
rxx_p = Rxx_p ./ Rxx_p(1);
y = [rxx_n,rxx_p];
figure(1);
hold on;
stem(k,y);
title('rxx(k) estimadas a partir de la funcion muestra')
grid on;
hold off;
Coefs = PartialCorrelation( rxx_p,len );
figure(2);
hold on;
stem(Coefs);
grid on;
title('Coeficientes de autocorrelacion parcial')
hold off;
%Calculo los parametros para un modelo AR de orden 2 a partir
%de los estimadores.
phi_1 = rxx_p(2)*( 1-rxx_p(3) );
phi_1 = phi_1 / ( 1 - (rxx_p(2)^2) );
phi_2 = Coefs(2);
Var_x = 1.0/ (1-phi_1*rxx_p(2)-phi_2*rxx_p(3));
rxx_teorico = zeros(1,len);
rxx_teorico(1) = 1;
rxx_teorico(2) = phi_1 / (1-phi_2);
for i = 3:len
    rxx_teorico(i) = phi_1*rxx_teorico(i-1) + phi_2*rxx_teorico(i-2);
end
figure(3)
subplot(2,1,1);
stem(rxx_p);
title('Autocorrelacion normalizada estimada');

subplot(2,1,2);
stem(rxx_teorico);
title('Autocorrelacion del modelo AR de orden 2');
%Espectrograma utilizando periodograma.
X_period = (abs( fft(x) ).^2) ./ ( len_x );
%Espectograma transformand la rxx
Sxx = abs( fft([fliplr(Rxx_n(2:end)),Rxx_n]) );
len_S = size(Sxx);
len_S = len_S(2);
figure(4)
hold on;
f = linspace(0,1,len_x);
plot(f,X_period)
f = linspace(0,1,len_S);
plot(f,Sxx)
xlabel('f(Hz)');
ylabel('Sxx')
title('Densidad espectral de potencia estimada.')
xlim([0 0.5]);