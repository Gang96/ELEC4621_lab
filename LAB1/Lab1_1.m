clc;
clear;

%-------------------------Q(a)-------------------------
%The minimum sampling frequency should be twice as the one of the 
%original signal, so the fmin = 2*f = 4400Hz

%-------------------------Q(b)-------------------------
%initialize variables
f = 2200;
fmin = 2*f;
t = (0:1e-6:5/2200);
x = cos(2*pi*f*t); %continuous signal
n2 = (0:1/(1.6*fmin):5/2200);
n3 = (0:1/(0.6*fmin):5/2200);
x2 = cos(2*pi*f*n2); %discrete signal for Q(c)
x3 = cos(2*pi*f*n3); %discrete signal for Q(d)

%plot the original signal
figure(1);
subplot(3,1,1);
plot(t,x);
grid on;
ylim([-1.5 1.5]);

%--------------------------Q(c)-----------------------

%plot the wave for Q(c)
subplot(3,1,2);
stem(n2,x2);
grid on;
ylim([-1.5 1.5]);

%1.6*2 = 3.2 samples we get for each period
%The minimum frequency to fit the points is 2200 Hz

%--------------------------Q(d)-----------------------

%plot the wave for Q(d)
subplot(3,1,3);
stem(n3,x3);
grid on;
ylim([-1.5 1.5]);

%0.6*2 = 1.2 samples we get for each period
%The minimum frequency to fit the points is 2200/5 = 440 Hz

%------------------------------Q(e)---------------------------
%original signal's spectra
n4 = -15000:15000;
h = zeros(length(n4),1);
for w = 1:length(n4)
    for k = 1:length(x)
    h(w) = h(w)+x(k)*exp(-1*1i*2*pi*n4(w)*t(k));
    end
end
figure(2);
subplot(3,1,1);
plot(n4,abs(h));

%original Q(c)'s signal's spectra
hc = zeros(length(n4),1);
for w = 1:length(n4)
    for k = 1:length(x2)
    hc(w) = hc(w)+x2(k)*exp(-1*1i*2*pi*n4(w)*n2(k));
    end
end
subplot(3,1,2);
plot(n4,abs(hc));

%original Q(d)'s signal's spectra
hd = zeros(length(n4),1);
for w = 1:length(n4)
    for k = 1:length(x3)
    hd(w) = hd(w)+x3(k)*exp(-1*1i*2*pi*n4(w)*n3(k));
    end
end
subplot(3,1,3);
plot(n4,abs(hd));







