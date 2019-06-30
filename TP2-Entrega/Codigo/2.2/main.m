snr = 30;
N = 20;
audioPath = 'whereIam8Khz.wav';
GpArray = zeros(1,N);
for n = 1.0:1.0:N
    [~,~,~,~,Gp,fs] = TP2_ej2(n,snr,audioPath);
    GpArray(n) = Gp;
end
fprintf('fs: %i Hz\n',fs);
figure('Name','Ganancias de Predicción');
plot(1.0:1.0:N,GpArray);
xlabel('N (orden del predictor)');
ylabel('Gp (veces)');
[GpMax,orden] = max(GpArray); % encuentro orden optimo
[x,s,sp,e,Gp,fs] = TP2_ej2(orden,snr,audioPath);
t = 0:seconds(1/fs):seconds(size(s,1)/fs);
t = t(1:size(s,1));
figure('Name','Señal + Ruido vs Señal original');
hold on;
plot(t,x);
plot(t,s);
hold off;
figure('Name','Señal + Ruido vs Señal predecida');
hold on;
plot(t,x);
plot(t,sp);
hold off;
xlabel('Time (s)');
ylabel('Audio Signal');
figure('Name','Señal Original vs Señal predecida');
hold on;
plot(t,sp);
plot(t,s);
hold off;
xlabel('Time (s)');
ylabel('Audio Signal');
figure('Name','Señal Error');
plot(t,e);
xlabel('Time (s)');
ylabel('Audio Signal');
disp('Press enter to listen to original audio');
pause;
sound(s,fs);
disp('Press enter to continue after the audio ends');
pause;
disp('Press enter to listen to original audio + AWGN');
pause;
sound(x,fs);
disp('Press enter to continue after the audio ends');
pause;
disp('Press enter to listen to estimated audio');
pause;
sound(sp,fs);
disp('Press enter to continue after the audio ends');
pause;


