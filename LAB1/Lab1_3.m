clc;
clear;

%---------------------------Q(a)-----------------------------
%represent the two filters, suppose the sampling time to be 10^-6s
a = 1/2;
Ts = 1e-6;
H = tf([1 -1*a],1,Ts,'Variable','z^-1');
G = tf(1,[1 -1*a],Ts,'Variable','z^-1');

%get the impulse responses with a = 1/2
figure(1);
subplot(2,1,1);
impulse(H);
subplot(2,1,2);
impulse(G);

%get the impulse responses with a = -1/2
a = -1/2;
H = tf([1 -1*a],1,Ts,'Variable','z^-1');
G = tf(1,[1 -1*a],Ts,'Variable','z^-1');
figure(2);
subplot(2,1,1);
impulse(H);
subplot(2,1,2);
impulse(G);

%---------------------------Q(b)-----------------------------
%present the input signal
f = 2200;
n = (0:Ts:5/2200);
x = cos(2*pi*f*n);

%To simplify the calculation, the transfer function can
%be represented by difference equations.
%first, try a = 1/2
y = zeros(length(n),1);
w = zeros(length(n),1);
a = 1/2;
for k = 1:length(n)
    if k == 1
        y(k) = x(k);
    else
        y(k) = x(k)-a*x(k-1);
    end
end
for k = 1:length(n)
    if k == 1
        w(k) = y(k);
    else
        w(k) = y(k)+a*w(k-1);
    end
end

%draw the output 
figure(3);
subplot(3,1,1);
plot(n,x);
subplot(3,1,2);
plot(n,y);
subplot(3,1,3);
plot(n,w);

%then, try a = -1/2
y2 = zeros(length(n),1);
w2 = zeros(length(n),1);
a = -1/2;
for k = 1:length(n)
    if k == 1
        y2(k) = x(k);
    else
        y2(k) = x(k)-a*x(k-1);
    end
end
for k = 1:length(n)
    if k == 1
        w2(k) = y2(k);
    else
        w2(k) = y2(k)+a*w2(k-1);
    end
end

%draw the output when a = 1/2
figure(4);
subplot(3,1,1);
plot(n,x);
subplot(3,1,2);
plot(n,y2);
subplot(3,1,3);
plot(n,w2);

%---------------------------Q(c)-----------------------------
n2 = [-pi:pi/100:pi];
h2 = zeros(length(n2),1);
for w2 = 1:length(n2)
    for k = 1:length(y)
    h2(w2) = h2(w2)+y(k)*exp(-1*1i*n2(w2)*(k-1));
    end
end
figure(5);
plot(n2,abs(h2));

g2 = zeros(length(n2),1);
for w2 = 1:length(n2)
    for k = 1:length(w)
    g2(w2) = g2(w2)+w(k)*exp(-1*1i*n2(w2)*(k-1));
    end
end
figure(6);
plot(n2,abs(g2));
















