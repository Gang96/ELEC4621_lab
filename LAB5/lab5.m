clear;
clc;
%% ----------------------Q(1)------------------
%Power Spetral Density(PSD) is the DTFT of the auto-correlation function.
%So to get the PSD of x[n], we need first get the auto-correlation function
%of x[n]. Then use DTFT with R[n] we will get the SPD, which is the job we
%will done in the next several questions.

%% ----------------------Q(2)------------------
n = 1:64;
R = zeros(1,64);
SNR = 100; %assume SNR(sigma to noise ratio) is 100
sigma = sqrt(339^2/(10^(SNR/10)));
x = 339*exp(1i*2*pi/10*n)+sigma/sqrt(2)*(randn(1,64)+1i*randn(1,64));
for i = 1:64
    for count = 1:(65-i)
    R(i) = R(i)+conj(x(count))*x(count+i-1);
    end
    R(i) = R(i)/(65-i);
end
figure(1)
plot(n,abs(R));
xlabel('n');
ylabel('Magnitude');
title('Unbiased Auto-correlation');
grid on;

%% ----------------------Q(3)------------------------
R_2 = zeros(1,64);
for i = 1:64
    for count = 1:(65-i)
    R_2(i) = R_2(i)+conj(x(count))*x(count+i-1);
    end
    R_2(i) = R_2(i)/64;
end
figure(2)
plot(n,abs(R_2));
xlabel('n');
ylabel('Magnitude');
title('Biased Auto-correlation');
grid on;

%% ---------------------Q(4)--------------------------
n2 = -pi:2*pi/300:pi;
PSD = zeros(1,length(n2));
for w = 1:length(n2)
    for k = 1:length(R)
    PSD(w) = PSD(w)+R(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

PSD_2 = zeros(1,length(n2));
for w = 1:length(n2)
    for k = 1:length(R_2)
    PSD_2(w) = PSD_2(w)+R_2(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

figure(3)
subplot(2,1,1);
plot(n2,abs(PSD));
xlabel('w');
ylabel('PSD');
title('PSD of Unbiased Autocorrelation');
grid on;

subplot(2,1,2);
plot(n2,abs(PSD_2));
xlabel('w');
ylabel('PSD');
title('PSD of Biased Autocorrelation');
grid on;

%From figure(4), we can easily find that the unbiased PSD is more curly
%than the one of the biased autocorrelation. In other words, the frequency
%components of biased PSD is much leaner than the unbiased PSD, this is
%because in biased PSD, the amplitude of other frequency except 50Hz has
%been weakened. But in the unbasid PSD, the autocorrelation function has
%the same amplitude. So it will get more curly.
%% -------------------Q(5)--------------------
% We will get N-M+1 blocks
M=32;
blocks=64-M+1;
x_block=zeros(blocks,M);
R_1_block=zeros(blocks,2*M-1);
R_2_block=zeros(blocks,2*M-1);
for i=1:blocks
     x_block(i,:)=x(i:i+M-1);
     for i2 = 1:(2*M-1)
         for count = 1:(2*M-i2)
             R_1_block(i,i2) = R_1_block(i,i2)+conj(x(count))*x(count+i2-1);
         end
         R_1_block(i,i2) = R_1_block(i,i2)/(2*M-i2);
     end
     for i2 = 1:(2*M-1)
         for count = 1:(2*M-i2)
             R_2_block(i,i2) = R_2_block(i,i2)+conj(x(count))*x(count+i2-1);
         end
         R_2_block(i,i2) = R_2_block(i,i2)/(2*M);
     end
end
R_1_aver=mean(R_1_block,1);
R_2_aver=mean(R_2_block,1);

PSD_b_1 = zeros(1,length(n2));
for w = 1:length(n2)
    for k = 1:length(R_1_aver)
    PSD_b_1(w) = PSD_b_1(w)+R_1_aver(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

PSD_b_2 = zeros(1,length(n2));
for w = 1:length(n2)
    for k = 1:length(R_2_aver)
    PSD_b_2(w) = PSD_b_2(w)+R_2_aver(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

figure(4)
subplot(2,1,1);
plot(n2,abs(PSD_b_1));
xlabel('w');
ylabel('PSD');
title('PSD of Unbiased Autocorrelation (BLOCK METHOD)');
grid on;

subplot(2,1,2);
plot(n2,abs(PSD_b_2));
xlabel('w');
ylabel('PSD');
title('PSD of Biased Autocorrelation (BLOCK METHOD)');
grid on;

%% ------------------Q(7)---------------
A=1;
SNR_2=-20:2:30;
N=128;
var=A^2./(10.^(SNR_2./10)); 
p=A^2./var;
CRB_f=6*500^2/(2*pi)^2./p./N./(N^2-1);
CRB_A=var./N;
figure(5)
subplot(2,1,1);
plot(SNR_2,CRB_f);
title('CRB of the frequency versus SNR');
xlabel('SNR');
ylabel('CRB');
subplot(2,1,2);
plot(SNR_2,CRB_A);
title('CRB of the amplitude versus SNR');
xlabel('SNR');
ylabel('CRB');

%% --------------------------Q(8)--------------------
n_2 = 1:128;
SNR_3 = -20:2:30;
f_estimate = zeros(1,1000);
A_estimate = zeros(1,1000);
f_MSE = zeros(1,length(SNR_3));
A_MSE = zeros(1,length(SNR_3));
for count_2 = 1:length(SNR_3)
    for count_3 = 1:1000
        sigma = sqrt(1^2/(10^(SNR_3(count_2)/10)));
        f = 45+10*rand(1,1);
        x_2 = exp(1i*2*pi/500*f*n_2)+sigma/sqrt(2)*(randn(1,128)+1i*randn(1,128));
        x_w = fft(x_2,2048);
        [a,b]=max(abs(x_w).^2);
        b = b/2048*500;
        f_estimate(count_3)=b;
        A_estimate(count_3)=max(abs(x_w))/length(x_2);
        f_MSE(count_2)=f_MSE(count_2)+(f_estimate(count_3)-50).^2;
        A_MSE(count_2)=A_MSE(count_2)+(A_estimate(count_3)-1).^2;
    end
    f_MSE(count_2) = f_MSE(count_2)/1000;
    A_MSE(count_2) = A_MSE(count_2)/1000;
end

figure(6);
subplot(2,1,1);
semilogy(SNR_3,CRB_f,SNR_3,f_MSE);
title('frequency');
xlabel('SNR');
ylabel('MSE&CRB');
legend('CRB','MSE');
subplot(2,1,2);
semilogy(SNR_3,CRB_A,SNR_3,A_MSE);
title('amplitude');
xlabel('SNR');
ylabel('MSE&CRB');
ylim([10^(-6),10^(1)]);
legend('CRB','MSE');

%Briefly speaking the MSE of both frequency and amplitude will get smaller
%with the smaller noise.











