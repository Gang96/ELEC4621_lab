clc;
clear;
%% -----------Q(a)-------------
% From question we know:
% H(z)=1+0.74417*z^(-1)+0.52604*z^(-2)+0.625*z^(-3)-0.1296*z(-4)
% The Nyquist gain is the h(w) when w = pi. So it is 0.02727.
%% -----------Q(b)-------------
h1 = [1,0.74417,0.52604,0.625,-0.1296];
zero1 = roots(h1);
%So, the zeros are zi = [-0.9873 0.0345+1i*0.8677 0.0345-1i*0.8677 0.1741]

%% -----------Q(c)-------------
% Calculate to get the 7 other filters' transfer functions
zero2 = zero1;
zero2(1) = 1/zero1(1);
h2 = poly(zero2);

zero3 = zero1;
zero3(4) = 1/zero1(4);
h3 = poly(zero3);

zero4 = zero1;
zero4(2) = 1/zero1(2);
zero4(3) = 1/zero1(3);
h4 = poly(zero4);

zero5 = zero1;
zero5(1) = 1/zero1(1);
zero5(4) = 1/zero1(4);
h5 = poly(zero5);

zero6 = zero1;
zero6(1) = 1/zero1(1);
zero6(2) = 1/zero1(2);
zero6(3) = 1/zero1(3);
h6 = poly(zero6);

zero7 = zero1;
zero7(4) = 1/zero1(4);
zero7(2) = 1/zero1(2);
zero7(3) = 1/zero1(3);
h7 = poly(zero7);

zero8 = zero1;
zero8(1) = 1/zero1(1);
zero8(4) = 1/zero1(4);
zero8(2) = 1/zero1(2);
zero8(3) = 1/zero1(3);
h8 = poly(zero8);

% Get the magnitude responses for the 8 different filters
n = -pi:2*pi/100:pi;
h1_2 = zeros(length(n),1);
for w = 1:length(n)
    for k = 1:length(h1)
   h1_2(w) = h1_2(w)+h1(k)*exp(-1*1i*n(w)*(k-1));
    end
end

h2_2 = zeros(length(n),1);
for w = 1:length(n)
    for k = 1:length(h2)
   h2_2(w) = h2_2(w)+h2(k)*exp(-1*1i*n(w)*(k-1));
    end
end

h3_2 = zeros(length(n),1);
for w = 1:length(n)
    for k = 1:length(h3)
   h3_2(w) = h3_2(w)+h3(k)*exp(-1*1i*n(w)*(k-1));
    end
end

h4_2 = zeros(length(n),1);
for w = 1:length(n)
    for k = 1:length(h4)
   h4_2(w) = h4_2(w)+h4(k)*exp(-1*1i*n(w)*(k-1));
    end
end

h5_2 = zeros(length(n),1);
for w = 1:length(n)
    for k = 1:length(h5)
   h5_2(w) = h5_2(w)+h5(k)*exp(-1*1i*n(w)*(k-1));
    end
end

h6_2 = zeros(length(n),1);
for w = 1:length(n)
    for k = 1:length(h6)
   h6_2(w) = h6_2(w)+h6(k)*exp(-1*1i*n(w)*(k-1));
    end
end

h7_2 = zeros(length(n),1);
for w = 1:length(n)
    for k = 1:length(h7)
   h7_2(w) = h7_2(w)+h7(k)*exp(-1*1i*n(w)*(k-1));
    end
end

h8_2 = zeros(length(n),1);
for w = 1:length(n)
    for k = 1:length(h8)
   h8_2(w) = h8_2(w)+h8(k)*exp(-1*1i*n(w)*(k-1));
    end
end

figure(1);
plot(n,abs(h1_2)/max(abs(h1_2)),'k-o',n,abs(h2_2)/max(abs(h2_2)),'-x',n,abs(h3_2)/max(abs(h3_2)),'-+',n,abs(h4_2)/max(abs(h4_2)),'-*',n,abs(h5_2)/max(abs(h5_2)),'-s',n,abs(h6_2)/max(abs(h6_2)),'-d',n,abs(h7_2)/max(abs(h7_2)),'-v',n,abs(h8_2)/max(abs(h8_2)),'-p');
grid on;
xlabel('n');
ylabel('magnitude');
legend('h1','h2','h3','h4','h5','h6','h7','h8');
% Notice that they are the same

%% -----------Q(d)-------------
figure(2);
plot(n,angle(h1_2)/pi*360,'k-o',n,angle(h2_2)/pi*360,'-x',n,angle(h3_2)/pi*360,'-+',n,angle(h4_2)/pi*360,'-*',n,angle(h5_2)/pi*360,'-s',n,angle(h6_2)/pi*360,'-d',n,angle(h7_2)/pi*360,'-v',n,angle(h8_2)/pi*360,'-p');
grid on;
xlabel('n');
ylabel('phase');
legend('h1','h2','h3','h4','h5','h6','h7','h8');
% h1 is the minimum phase with the smoothest phase response.

%% ------------Q(e)-------------
figure(3);
m = n(1:100);
plot(m,diff(angle(h1_2)/pi*360)/(pi/100),'k-o',m,diff(angle(h2_2)/pi*360)/(pi/100),'-x',m,diff(angle(h3_2)/pi*360)/(pi/100),'-+',m,diff(angle(h4_2)/pi*360)/(pi/100),'-*',m,diff(angle(h5_2)/pi*360)/(pi/100),'-s',m,diff(angle(h6_2)/pi*360)/(pi/100),'-d',m,diff(angle(h7_2)/pi*360)/(pi/100),'-v',m,diff(angle(h8_2)/pi*360)/(pi/100),'-p');
grid on;
xlabel('n');
ylabel('phase');
legend('h1','h2','h3','h4','h5','h6','h7','h8');
% The minimum phase system h1 has the smallest group delay.


