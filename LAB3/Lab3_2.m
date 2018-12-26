clear;
clc;
%% ---------------------------Q(a)-----------------------------
% From the poles give we know the filter's transfer function is
% H(z) =
% z^4/((z+0.55)*(z+0.65)*(z-0.98*exp(1i*0.67*pi)*(z-0.98*exp(-1i*0.67*pi)))

%% ---------------------------Q(b)-----------------------------
% All-pole filter so ai = 0. After simplify the brackets we get
% b3 = -1.2; b4 = -0.3575; b1 = -0.98; b2 = -0.9604;

%% ---------------------------Q(d)-----------------------------
pole = [-0.55; -0.65; 0.98*exp(1i*0.67*pi); 0.98*exp(-1i*0.67*pi)];
a = poly(pole);
X = zeros(1000,1);
X(1) = 1;
Y = filter(1,a,X);

G = sum(abs(Y));

pole1 = [-0.55; -0.65];
a1 = poly(pole1);
Y1 = filter(1,a1,X);

G1 = sum(abs(Y1));
% After knowing the BIBO gain, we can get the representation of these accumulator blocks
% Block 1: (4.4);
% Block 2: (7.1);
% Block 3: (7.1);

%% -------------------------Q(e)-------------------------
pole2(1) = pole(3);
pole2(2) = pole(4);
a2 = poly(pole2);
Y2 = filter(1,a2,X);

P1 = sum(Y.^2);
P2 = sum(Y2.^2);
P3 = 1;

% Noise power are:
E1 = 1/3072*P1;
E2 = 1/48*P2;
E3 = 1/48*P3;

%% -------------------------Q(f)-------------------------
w = -pi:pi/100:pi;
z = exp(1i*w);
T1 = z.^4./((z+0.55).*(z+0.65).*(z-0.98*exp(1i*0.67*pi)).*(z-0.98*exp(-1i*0.67*pi)));
T1_2 = (z.^(-1)).^4./((z.^(-1)+0.55).*(z.^(-1)+0.65).*(z.^(-1)-0.98*exp(1i*0.67*pi)).*(z.^(-1)-0.98*exp(-1i*0.67*pi)));
T2 = z.^2./((z-0.98*exp(1i*0.67*pi)).*(z-0.98*exp(-1i*0.67*pi)));
T2_2 = (z.^(-1)).^2./(((z.^(-1))-0.98*exp(1i*0.67*pi)).*((z.^(-1))-0.98*exp(-1i*0.67*pi)));

N1 = 1/3072.*T1.*T1_2;
N2 = 1/48.*T2.*T2_2;
N3 = E3;
N = N1+N2+N3;
figure(4);
plot(w/pi, abs(N));
grid on;







