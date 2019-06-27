x = load('Archivo_4.mat');
x= x.x;
len_x = size(x);
len_x = len_x(2);
len = 128;
k = 0:len-1;
Rxx_n = NpCorrelationEstimate(x,len);
rxx_n = Rxx_n ./ Rxx_n(1);
Rxx_p = PolCorrelationEstimate( x,len );
rxx_p = Rxx_p ./ Rxx_p(1);
%Grafica de las autocorrelaciones estimadas.
figure(1)
hold on;
subplot(2,1,1);
stem(k,Rxx_n);
ylabel('Rxx(n)');
xlabel('n');
grid on;
title('Autocorrelacion estimada con estimador no polarizado');

subplot(2,1,2);
stem(k,Rxx_p);
ylabel('Rxx(n)')
xlabel('n')
title('Autocorrelacion estimada con estimador polarizado');
grid on;
hold off;

figure(2);
hold on;
subplot(2,1,1);
stem(k,rxx_n);

ylabel('rxx(n)')
xlabel('n')
title('Autocorrelacion normalizada, estimador no polarizado');
grid on;

subplot(2,1,2);
stem(k,rxx_p);

ylabel('rxx(n)')
xlabel('n')
title('Autocorrelacion estimada,estimador polarizado');
grid on;
hold off;


Coefs = PartialCorrelation( rxx_p,len-1 );
figure(3);
hold on;
stem(Coefs);
grid on;
title('Coeficientes de autocorrelacion parcial')
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
figure(4)
subplot(2,1,1);
stem(rxx_p);
grid on;
title('Autocorrelacion normalizada estimada');

subplot(2,1,2);
stem(rxx_teorico);
grid on;
title('Autocorrelacion del modelo AR de orden 2');
%Espectrograma utilizando periodograma.
X_period = (abs( fft(x) ).^2) ./ ( len_x );
%Espectograma transformand la rxx
Sxx = abs( fft([fliplr(Rxx_n(2:end)),Rxx_n]) );
len_S = size(Sxx);
len_S = len_S(2);
figure(5)
hold on;
f = linspace(0,1,len_x);
plot(f,X_period)
f = linspace(0,1,len_S);
plot(f,Sxx)
grid on;
xlabel('f(Hz)');
ylabel('Sxx')
title('Densidad espectral de potencia estimada.')
xlim([0 0.5]);