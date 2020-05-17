%% Z
clear, clc, close all;

data = load("C:\Users\milja\source\repos\DOSR\DOSR\data\ekg0.mat");
x = data.ekg;
N = length(x);
fs = data.fs;
F=0:1/fs:(length(x)-1)/fs;
%figure(1);
%plot(F(1:20000),x(1:20000));
%title('Loaded ECG data')
%ylabel('x[t]')
%xlabel('t[s]')

fosa=0:fs/N:(fs/2);
X=fft(x,N)/length(x)*2;
Xamp=abs(X(1:(N/2+1)));
figure (2);
plot(fosa,Xamp);
xlabel('f[Hz]');
ylabel('|X(jf)|[Db]');
%{
fcutlow = 2;
fcuthigh = 400;
Wp = [fcutlow fcuthigh]/500;                                 % Passband Frequency (Normalised)
Ws = [fcutlow-1 fcuthigh+1]/500;                             % Stopband Frequency (Normalised)
Rp =   2;                                                   % Passband Ripple (dB)
Rs =  40;                                                   % Stopband Ripple (dB)
[n,Ws]  = cheb1ord(Wp,Ws,Rp,Rs);                            % Filter Order
[z,p,k] = cheby1(n,Rs,Ws);                                  % Filter Design, Sepcify Bandpass
[sos,g] = zp2sos(z,p,k);                                    % Convert To Second-Order-Section For Stability
figure(3);
freqz(sos, 2^16, fs)                                        % Filter Bode Plot
signal_Filtered = filtfilt(sos, g, x);


%{
Rp=2;
Rs=40;
Wp = [45 55]/500;
Ws = [1 500]/500;
[n,Wn] = buttord(Wp,Ws,Rp,Rs,'s');

[b,a]=butter(n,Wn,'s');
h=freqs(b,a,length(fosa));

figure(3);
plot(fosa,20*log10(abs(h)));
xlabel('Frequency, Hz');
ylabel('Magnitude, dB');
[bd,ad] = bilinear(b,a,fs);    
%}

%{
Wp = [100 200]/500;
Ws = [50 250]/500;
Rp = 3;
Rs = 40;
[n,Wn] = buttord(Wp,Ws,Rp,Rs)
[b,a] = butter(n,Wn);
sos = zp2sos(b,a);

freqz(sos,128,1000)
title(sprintf('n = %d Butterworth Bandpass Filter',n))

[bd,ad] = bilinear(b,a,fs);    
%}
disp(signal_Filtered)
%y = filter(bd,ad,signal_Filtered);
figure (4);
plot(F, signal_Filtered);
xlabel('Frequency, Hz');
ylabel('Magnitude, dB');
disp('end');
X1=fft(signal_Filtered,N)/length(signal_Filtered)*2;
Xamp=abs(X1(1:(N/2+1)));
figure (5);
plot(fosa,Xamp);
xlabel('f[Hz]');
ylabel('|X(jf)|[Db]');

%[b,a]=butter(n,Wn,'s');
%w = logspace(-1,1);
%h=freqs(b,a,w);
%plot(1,20*log10(abs(h)));

%}

%% butter


clear all; close all; clc; 


data = load("C:\Users\milja\source\repos\DOSR\DOSR\data\ekg0.mat");
x = data.ekg;
Fs = data.fs;
t=0:1/Fs:(length(x)-1)/Fs;

N = 4*2^nextpow2(length(x));
f1 = 0:Fs/N:Fs/2;
X = fft(x,N)/length(x); 
X1 = abs(X(1:N/2+1));
X1(2:N/2+1) = 2*X1(2:N/2+1); 
Xphase = unwrap(angle(X(1:N/2+1))); 

figure
subplot(2, 2, 1:2)
    plot(t, x); 
    xlabel('t[s]'); ylabel('x(t)');
    title('Signal x(t)'); grid on;
subplot(2, 2, 3)
    plot(f1, X1); xlim([0 Fs/2]); 
    xlabel('f[Hz]'); ylabel('|X(jf)|');
    title('Amplitudska karakteristika'); grid on;
subplot(2,2,4)
    plot(f1, Xphase); xlim([0 Fs/2]); ylim([-400 10]);
    xlabel('f[Hz]'); ylabel('arg{X(jf)}');
    title('Fazna karakteristika'); grid on
    
Wp = [0.001 499.999]*2*pi;
Ws = [49.99 50.01]*2*pi;
Rp = 2; Rs = 40;
[n, Wn] = buttord(Wp, Ws, Rp, Rs, 's');
[b, a] = butter(n, Wn, 's');

[h, w] = freqs(b, a, N/2+1);

[bz, az] = bilinear(b, a, Fs);

[hz, fz] = freqz(bz, az, N/2+1, Fs); 

figure
    plot(w/(2*pi), 20*log10(abs(h)), 'k-', 'Linewidth', 1.5); hold on; 
    plot(fz, 20*log10(abs(hz)), 'r', 'Linewidth', 1.5);
    xlabel('f [Hz]');
    title('Amplitudska karakteristika filtra'); grid on;
    legend('analogni', 'digitalni');
    
y = filter(bz, az, x); 
y = x-y;
Y = fft(y,N)/length(y); 
Y1 = abs(Y(1:N/2+1));
Y1(2:N/2+1) = 2*Y1(2:N/2+1); 


figure
    subplot(2, 2, [1, 3])
        plot(f1, X1, 'Linewidth', 1.5); hold on
        plot(f1, Y1, 'Linewidth', 1.5); 
        xlabel('f[Hz]'); ylabel('|X(jf)|, |Y(jf)|');
        legend('ulazni', 'izlazni signal')
        title('Amplitudske karakteristike ulaznog i izlaznog signala'); grid on;
    subplot(2, 2, 2)
        plot(t, x);
        xlabel('t[s]'); ylabel('x(t)');
        title('Ulazni signal'); grid on;
    subplot(2, 2, 4)
        plot(t, y, 'r');
        xlabel('t[s]'); ylabel('y(t)');
        title('Izlazni signal'); grid on;


%%  cheby1


clear all; close all; clc; 


data = load("C:\Users\milja\source\repos\DOSR\DOSR\data\ekg0.mat");
x = data.ekg;
Fs = data.fs;
t=0:1/Fs:(length(x)-1)/Fs;

N = 4*2^nextpow2(length(x));
f1 = 0:Fs/N:Fs/2;
X = fft(x,N)/length(x); 
X1 = abs(X(1:N/2+1));
X1(2:N/2+1) = 2*X1(2:N/2+1); 
Xphase = unwrap(angle(X(1:N/2+1))); 

figure
subplot(2, 2, 1:2)
    plot(t, x); 
    xlabel('t[s]'); ylabel('x(t)');
    title('Signal x(t)'); grid on;
subplot(2, 2, 3)
    plot(f1, X1); xlim([0 Fs/2]); 
    xlabel('f[Hz]'); ylabel('|X(jf)|');
    title('Amplitudska karakteristika'); grid on;
subplot(2,2,4)
    plot(f1, Xphase); xlim([0 Fs/2]); ylim([-400 10]);
    xlabel('f[Hz]'); ylabel('arg{X(jf)}');
    title('Fazna karakteristika'); grid on
    
Wp = [49.5 50.5]*2*pi;
Ws = [49.999 50.001]*2*pi;
Rp = 2; Rs = 40;
[n, Wn] = cheb1ord(Wp, Ws, Rp, Rs, 's');
[b, a] = cheby1(n, Rp, Wp,'s');

[h, w] = freqs(b, a, N/2+1);

[bz, az] = bilinear(b, a, Fs);

[hz, fz] = freqz(bz, az, N/2+1, Fs); 

figure
    plot(w/(2*pi), 20*log10(abs(h)), 'k-', 'Linewidth', 1.5); hold on; 
    plot(fz, 20*log10(abs(hz)), 'r', 'Linewidth', 1.5);
    xlabel('f [Hz]');
    title('Amplitudska karakteristika filtra'); grid on;
    legend('analogni', 'digitalni');
    
y = filter(bz, az, x); 
y = x-y;
Y = fft(y,N)/length(y); 
Y1 = abs(Y(1:N/2+1));
Y1(2:N/2+1) = 2*Y1(2:N/2+1); 


figure
    subplot(2, 2, [1, 3])
        plot(f1, X1, 'Linewidth', 1.5); hold on
        plot(f1, Y1, 'Linewidth', 1.5); 
        xlabel('f[Hz]'); ylabel('|X(jf)|, |Y(jf)|');
        legend('ulazni', 'izlazni signal')
        title('Amplitudske karakteristike ulaznog i izlaznog signala'); grid on;
    subplot(2, 2, 2)
        plot(t, x);
        xlabel('t[s]'); ylabel('x(t)');
        title('Ulazni signal'); grid on;
    subplot(2, 2, 4)
        plot(t, y, 'r');
        xlabel('t[s]'); ylabel('y(t)');
        title('Izlazni signal'); grid on;

%% experiment

clear all; close all; clc; 


data = load("C:\Users\milja\source\repos\DOSR\DOSR\data\ekg0.mat");
x = data.ekg;
Fs = data.fs;
t=0:1/Fs:(length(x)-1)/Fs;

N = 4*2^nextpow2(length(x));
f1 = 0:Fs/N:Fs/2;
X = fft(x,N)/length(x); 
X1 = abs(X(1:N/2+1));
X1(2:N/2+1) = 2*X1(2:N/2+1); 
Xphase = unwrap(angle(X(1:N/2+1))); 

figure
subplot(2, 2, 1:2)
    plot(t, x); 
    xlabel('t[s]'); ylabel('x(t)');
    title('Signal x(t)'); grid on;
subplot(2, 2, 3)
    plot(f1, X1); xlim([0 Fs/2]); 
    xlabel('f[Hz]'); ylabel('|X(jf)|');
    title('Amplitudska karakteristika'); grid on;
subplot(2,2,4)
    plot(f1, Xphase); xlim([0 Fs/2]); ylim([-400 10]);
    xlabel('f[Hz]'); ylabel('arg{X(jf)}');
    title('Fazna karakteristika'); grid on
    
Fs = 1000;  % Sampling Frequency

Fpass1 = 47;          % First Passband Frequency
Fstop1 = 49;          % First Stopband Frequency
Fstop2 = 51;          % Second Stopband Frequency
Fpass2 = 53;          % Second Passband Frequency
Apass1 = 2;           % First Passband Ripple (dB)
Astop  = 40;          % Stopband Attenuation (dB)
Apass2 = 1;           % Second Passband Ripple (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);

figure
    plot(w/(2*pi), 20*log10(abs(h)), 'k-', 'Linewidth', 1.5); hold on; 
    plot(fz, 20*log10(abs(hz)), 'r', 'Linewidth', 1.5);
    xlabel('f [Hz]');
    title('Amplitudska karakteristika filtra'); grid on;
    legend('analogni', 'digitalni');
    
y = filter(bz, az, x); 
y = x-y;
Y = fft(y,N)/length(y); 
Y1 = abs(Y(1:N/2+1));
Y1(2:N/2+1) = 2*Y1(2:N/2+1); 


figure
    subplot(2, 2, [1, 3])
        plot(f1, X1, 'Linewidth', 1.5); hold on
        plot(f1, Y1, 'Linewidth', 1.5); 
        xlabel('f[Hz]'); ylabel('|X(jf)|, |Y(jf)|');
        legend('ulazni', 'izlazni signal')
        title('Amplitudske karakteristike ulaznog i izlaznog signala'); grid on;
    subplot(2, 2, 2)
        plot(t, x);
        xlabel('t[s]'); ylabel('x(t)');
        title('Ulazni signal'); grid on;
    subplot(2, 2, 4)
        plot(t, y, 'r');
        xlabel('t[s]'); ylabel('y(t)');
        title('Izlazni signal'); grid on;
