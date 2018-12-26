clc;
clear;
%% -----------Q(a)------------
% A = 339V; f = 50Hz; fs = 500Hz
% x(t) = 339*exp(j*100*pi*t)
% FT:x(w) = 339*2*pi*delta(w-100*pi) = 678*pi*delta(w-100*pi)
% x(n) = 339*exp(j*2*pi*1/10*n) = 339*exp(j*pi/5*n)
%DTFT:x(w) = 678*pi*delta(w-pi/5)

%% -----------Q(b)------------
n = 1:64;
x = 339*exp(1i*2*pi/10*n);
x_r = 339*cos(pi/5*n);
x_i = 339*sin(pi/5*n);
figure(1);
subplot(2,1,1);
stem(n,x_r);
title('The real part of the signal');
xlabel('n/2*10^{-4}s');
ylabel('magnitude/V');
grid on;
subplot(2,1,2);
stem(n,x_i);
title('The imaginary part of the signal');
xlabel('n/2*10^{-4}s');
ylabel('magnitude/V');
grid on;

%% -----------Q(c)------------
k = 0:0.1:500-0.1;
X = zeros(1,length(k));
for count1 = 1:length(k)
    for count2 = 1:length(n)
        X(count1) = X(count1)+x(count2)*exp(-1*1i*2*pi*(count2-1)*k(count1)*10/length(k));
    end
end
figure(2);
plot(k,abs(X)/max(abs(X)));
title('The DFT of the signal');
xlabel('w/Hz');
ylabel('magnitude');
grid on;
% We can easily find when the w = 50Hz,DFT has the max coefficient.

%% ------------Q(d)-------------
X_m = fft(x,64);
k_2 = 0:500/64:500-1/64;
figure(3);
plot(k,abs(X)/max(abs(X)),k_2,abs(X_m)/max(abs(X_m)),'or');
title('The FFT(64) & DFT of the signal');
xlabel('w/Hz');
ylabel('magnitude');
legend('DFT','FFT');
grid on;
% From the plot, we can find there exists an error between the frequenct of
% the maximum value of FFT and 50, whitch is 50-46.88=3.12Hz.

%% ------------Q(e)-------------
X_m_2 = fft(x,512);
figure(4);
k_2=0:500/512:500-1/512;
plot(k,abs(X)/max(abs(X)),k_2,abs(X_m_2)/max(abs(X_m_2)),'or');
title('The FFT(512) & DFT of the signal');
xlabel('w/Hz');
ylabel('magnitude');
legend('DFT','FFT');
grid on;
%There still exists an error but it is smaller now. The error at this point
%is 50-49.8=0.2Hz, whitch means that when the L gets bigger, the frequency
%is more accurate.

%% ------------Q(f)-------------
fre = zeros(16,1);
ma = zeros(16,1);
count_3 = 0;
for L = 64:64:1024
    count_3 = count_3+1;
    X_m_3 = fft(x,L);
    [ma(count_3),fre(count_3)] = max(abs(X_m_3));
    fre(count_3) = (fre(count_3)-1)*500/L;
end
error = abs(fre-50);
L_p = 64:64:1024;
figure(5);
plot(L_p,error);
title('The frequency error');
xlabel('L');
ylabel('frequency/Hz');
grid on;
% Easy to find the bigger L is, the smaller error is.

%% ----------------Q(g)-------------
count_5 = 0;
figure(6);
for L = [64 128 256 512]
    count_5 = count_5+1;
    count_6 = 0;
    for SNR = -15:45
        count_6 = count_6+1;
        row(count_6) = SNR;
        sigma = sqrt(339^2/(10^(SNR/10)));
        sum = 0;
        for count_7 = 1:200
            x_n = x+sigma/sqrt(2)*(randn(1,64)+1i*randn(1,64));
            X_m_4 = fft(x_n,L);
            [ma_2(count_7),fre_2(count_7)] = max(abs(X_m_4));
            fre_2(count_7) = (fre_2(count_7)-1)*500/L;
            sum = sum+(fre_2(count_7)-50)^2;
        end
        MSE(count_6) = sum/200;
    end
    subplot(2,2,count_5);
    plot(row,10*log10(MSE));%use dB to see the difference more clearly
    xlabel('SNR/dB');
    ylabel('MSE/dB');
    title(['L is ',num2str(L)]);
    grid on;
end
% Noise will have a bad influence on the frequency estimation. Espacially
% when the space number L is small. Respactively, when the L is big enough,
% the noise will have much smaller influence on the estimation.













