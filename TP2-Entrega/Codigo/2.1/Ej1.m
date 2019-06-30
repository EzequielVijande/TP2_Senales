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
title('Autocorrelacion normalizada,estimador polarizado');
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
Var_x = (1-phi_2)/( (1+phi_2)*(1-phi_1-phi_2)*(1+phi_1-phi_2) );
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
title('Autocorrelacion normalizada del modelo AR de orden 2');

figure(5)
subplot(2,1,1);
stem(Rxx_p);
grid on;
title('Autocorrelacion estimada');

subplot(2,1,2);
stem(rxx_teorico.*Var_x);
grid on;
title('Autocorrelacion del modelo AR de orden 2');

%Espectrograma utilizando promediacion de periodogramas.
n=32;
chunk_size = len_x/n;
X_period = zeros(1,chunk_size);
for i = 1:n
    aux = x(i:i+chunk_size-1);
    X_period = X_period + (abs( fft(aux) ).^2) ./ ( n*chunk_size);
end
%Espectograma transformand la rxx
Sxx = abs( fft([fliplr(Rxx_n(2:end)),Rxx_n]) );
len_S = size(Sxx);
len_S = len_S(2);

figure(6)
hold on;
f = linspace(0,1,chunk_size);
plot(f,X_period,'DisplayName','Promediacion Periodograma')
f = linspace(0,1,len_S);
plot(f,Sxx,'DisplayName','Transfromada de la correlacion estimada')
grid on;
xlabel('f(Hz)');
ylabel('Sxx')
title('Densidad espectral de potencia estimada.')
xlim([0 0.5]);

%Espectrograma teorico de modelo AR orden 2
f = linspace(0,0.5,len_x);
arg1 = complex(0, -(2.*pi.*f) );
arg2 = complex(0, -(4.*pi.*f) );
Sxx_teorico = 1./( (abs(1-phi_1.*exp(arg1)-phi_2.*exp(arg2))).^2);
plot(f,Sxx_teorico,'DisplayName','Espectrograma teorico')
legend