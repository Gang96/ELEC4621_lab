clc;
clear;
%% -----------Q(a)-------------

%variables initialize
r1 = exp(-1/8);
r2 = 0.925;
r3 = r2;
R1 = r1;
R2 = r2*exp(1i*0.45*pi);
R3 = r2*exp(-1i*0.45*pi);
R4 = r3*exp(1i*0.8*pi);
R5 = r3*exp(-1i*0.8*pi);

% from the roots to calculate the coefficients of a polynomial
v = [R1;R2;R3;R4;R5];
h = poly(v);

% plot the impulse response
n = 1:6;
figure(1);
stem(n,h);

%% -----------Q(b)-------------
% use the matlab code in lab1
n2 = [-pi:2*pi/100:pi];
h2 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(h)
    h2(w) = h2(w)+h(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

%plot the magnitude and phase responses
figure(2);
subplot(2,1,1)
plot(n2,abs(h2));
subplot(2,1,2)
plot(n2,angle(h2));

%% -----------Q(c)&Q(d)-------------
% --------r1------------
r1 = exp(-1/8);
r2 = 0.925;
r3 = r2;
r1 = 1/r1;
R1 = r1;
R2 = r2*exp(1i*0.45*pi);
R3 = r2*exp(-1i*0.45*pi);
R4 = r3*exp(1i*0.8*pi);
R5 = r3*exp(-1i*0.8*pi);
v = [R1;R2;R3;R4;R5];
h_1 = poly(v);

n = 1:6;
figure(3);
subplot(3,1,1);
stem(n,h_1);

n2 = [-pi:2*pi/100:pi];
h2_1 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(h_1)
    h2_1(w) = h2_1(w)+h_1(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

%plot the magnitude and phase responses
subplot(3,1,2)
plot(n2,abs(h2_1));
subplot(3,1,3)
plot(n2,angle(h2_1));

% --------r2------------
r1 = exp(-1/8);
r2 = 0.925;
r3 = r2;
r2 = 1/r2;
R1 = r1;
R2 = r2*exp(1i*0.45*pi);
R3 = r2*exp(-1i*0.45*pi);
R4 = r3*exp(1i*0.8*pi);
R5 = r3*exp(-1i*0.8*pi);
v = [R1;R2;R3;R4;R5];
h_2 = poly(v);

n = 1:6;
figure(4);
subplot(3,1,1);
stem(n,h_2);

n2 = [-pi:2*pi/100:pi];
h2_2 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(h_1)
    h2_2(w) = h2_2(w)+h_2(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

%plot the magnitude and phase responses
subplot(3,1,2)
plot(n2,abs(h2_2));
subplot(3,1,3)
plot(n2,angle(h2_2));

% --------r3------------
r1 = exp(-1/8);
r2 = 0.925;
r3 = r2;
r3 = 1/r3;
R1 = r1;
R2 = r2*exp(1i*0.45*pi);
R3 = r2*exp(-1i*0.45*pi);
R4 = r3*exp(1i*0.8*pi);
R5 = r3*exp(-1i*0.8*pi);
v = [R1;R2;R3;R4;R5];
h_3 = poly(v);

n = 1:6;
figure(5);
subplot(3,1,1);
stem(n,h_3);

n2 = [-pi:2*pi/100:pi];
h2_3 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(h_1)
    h2_3(w) = h2_3(w)+h_3(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

%plot the magnitude and phase responses
subplot(3,1,2)
plot(n2,abs(h2_3));
subplot(3,1,3)
plot(n2,angle(h2_3));

% --------r1&r2------------
r1 = exp(-1/8);
r2 = 0.925;
r3 = r2;
r1 = 1/r1;
r2 = 1/r2;
R1 = r1;
R2 = r2*exp(1i*0.45*pi);
R3 = r2*exp(-1i*0.45*pi);
R4 = r3*exp(1i*0.8*pi);
R5 = r3*exp(-1i*0.8*pi);
v = [R1;R2;R3;R4;R5];
h_4 = poly(v);

n = 1:6;
figure(6);
subplot(3,1,1);
stem(n,h_4);

n2 = [-pi:2*pi/100:pi];
h2_4 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(h_1)
    h2_4(w) = h2_4(w)+h_4(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

%plot the magnitude and phase responses
subplot(3,1,2)
plot(n2,abs(h2_4));
subplot(3,1,3)
plot(n2,angle(h2_4));

% --------r2&r3------------
r1 = exp(-1/8);
r2 = 0.925;
r3 = r2;
r2 = 1/r2;
r3 = 1/r3;
R1 = r1;
R2 = r2*exp(1i*0.45*pi);
R3 = r2*exp(-1i*0.45*pi);
R4 = r3*exp(1i*0.8*pi);
R5 = r3*exp(-1i*0.8*pi);
v = [R1;R2;R3;R4;R5];
h_5 = poly(v);

n = 1:6;
figure(7);
subplot(3,1,1);
stem(n,h_5);

n2 = [-pi:2*pi/100:pi];
h2_5 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(h_1)
    h2_5(w) = h2_5(w)+h_5(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

%plot the magnitude and phase responses
subplot(3,1,2)
plot(n2,abs(h2_5));
subplot(3,1,3)
plot(n2,angle(h2_5));

% --------r1&r3------------
r1 = exp(-1/8);
r2 = 0.925;
r3 = r2;
r1 = 1/r1;
r3 = 1/r3;
R1 = r1;
R2 = r2*exp(1i*0.45*pi);
R3 = r2*exp(-1i*0.45*pi);
R4 = r3*exp(1i*0.8*pi);
R5 = r3*exp(-1i*0.8*pi);
v = [R1;R2;R3;R4;R5];
h_6 = poly(v);

n = 1:6;
figure(8);
subplot(3,1,1);
stem(n,h_6);

n2 = [-pi:2*pi/100:pi];
h2_6 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(h_1)
    h2_6(w) = h2_6(w)+h_6(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

%plot the magnitude and phase responses
subplot(3,1,2)
plot(n2,abs(h2_6));
subplot(3,1,3)
plot(n2,angle(h2_6));

% --------r1&r2&r3------------
r1 = exp(-1/8);
r2 = 0.925;
r3 = r2;
r1 = 1/r1;
r2 = 1/r2;
r3 = 1/r3;
R1 = r1;
R2 = r2*exp(1i*0.45*pi);
R3 = r2*exp(-1i*0.45*pi);
R4 = r3*exp(1i*0.8*pi);
R5 = r3*exp(-1i*0.8*pi);
v = [R1;R2;R3;R4;R5];
h_7 = poly(v);

n = 1:6;
figure(9);
subplot(3,1,1);
stem(n,h_7);

n2 = [-pi:2*pi/100:pi];
h2_7 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(h_1)
    h2_7(w) = h2_7(w)+h_7(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

%plot the magnitude and phase responses
subplot(3,1,2)
plot(n2,abs(h2_7));
subplot(3,1,3)
plot(n2,angle(h2_7));

%Question(c): Only the original filter has the minimum phase. And we can 
%know which response belongs to which system by simply analysing the phase
%response.

%Question(d): The original system has the minimum sum of the noms of
%the response i.e.sum(abs(h(n))), which means it has the minimum power.
%% -----------Q(e)-------------
%----------rounded to integer---------
h_8 = zeros(6,1);
h_8(1) = round(h(1),0);
h_8(2) = round(h(2),0);
h_8(3) = round(h(3),0);
h_8(4) = round(h(4),0);
h_8(5) = round(h(5),0);
h_8(6) = round(h(6),0);

n = 1:6;
figure(10);
subplot(3,1,1);
stem(n,h_8);

n2 = [-pi:2*pi/100:pi];
h2_8 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(h_1)
    h2_8(w) = h2_8(w)+h_8(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

%plot the magnitude and phase responses
subplot(3,1,2)
plot(n2,abs(h2_8));
subplot(3,1,3)
plot(n2,angle(h2_8));

%----------rounded to 1 number of decimal places----------
h_9 = zeros(6,1);
h_9(1) = round(h(1),1);
h_9(2) = round(h(2),1);
h_9(3) = round(h(3),1);
h_9(4) = round(h(4),1);
h_9(5) = round(h(5),1);
h_9(6) = round(h(6),1);

n = 1:6;
figure(11);
subplot(3,1,1);
stem(n,h_9);

n2 = [-pi:2*pi/100:pi];
h2_9 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(h_1)
    h2_9(w) = h2_9(w)+h_9(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

%plot the magnitude and phase responses
subplot(3,1,2)
plot(n2,abs(h2_9));
subplot(3,1,3)
plot(n2,angle(h2_9));

%----------rounded to 2 number of decimal places----------
h_10 = zeros(6,1);
h_10(1) = round(h(1),2);
h_10(2) = round(h(2),2);
h_10(3) = round(h(3),2);
h_10(4) = round(h(4),2);
h_10(5) = round(h(5),2);
h_10(6) = round(h(6),2);

n = 1:6;
figure(12);
subplot(3,1,1);
stem(n,h_10);

n2 = [-pi:2*pi/100:pi];
h2_10 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(h_1)
    h2_10(w) = h2_10(w)+h_10(k)*exp(-1*1i*n2(w)*(k-1));
    end
end

%plot the magnitude and phase responses
subplot(3,1,2)
plot(n2,abs(h2_10));
subplot(3,1,3)
plot(n2,angle(h2_10));

%Not much difference between the systems with the different numbers of
%decimal.
A = [1];
figure(13);
subplot(2,2,1);
zplane(h,A);
subplot(2,2,2);
zplane(h_8',A);
subplot(2,2,3);
zplane(h_9',A);
subplot(2,2,4);
zplane(h_10',A);






