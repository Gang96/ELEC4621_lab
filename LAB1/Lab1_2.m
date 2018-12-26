clc;
clear;
%----------------------------Q(a)--------------------------
%initialize variables

x = [1 2 3 4 5 6 7 6 5 4 3 2 1]';
h = [-0.0625 0.25 0.625 0.25 -0.0625]';
%in this part, use matrix to present sequences
h2 = [0 0 0 0 0 0 0 0 0 0 0 0 ...
    -0.0625 0.25 0.625 0.25 -0.0625 ...
    0 0 0 0 0 0 0 0 0 0 0 0]';
%set several zeros to avoid error
y = zeros(17,1);

%represent the convolution equation
for k = 1:(length(x)+4)
    for n = 1:length(x)
        y(k) = y(k)+x(n)*h2(k-n+13);
    end
end

%use convolution function to compare the results
yc = conv(x,h);

%--------------------------Q(b)--------------------------
y3 = zeros(17,1);
for k = 1:(length(x)+4)
    y3(k) = (h2(k+12:-1:k))'*x;
end
y3 = y3';

%--------------------------Q(c)--------------------------
n2 = [-pi:2*pi/100:pi];
h3 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(h)
    h3(w) = h3(w)+h(k)*exp(-1*1i*n2(w)*(k-1));
    end
end
figure(1);
subplot(2,1,1)
plot(n2,abs(h3));
subplot(2,1,2)
plot(n2,angle(h3));

%-------------------------Q(d)-----------------------------
%initialize variables

x = [1 2 3 4 5 6 7 6 5 4 3 2 1]';
hd = [-0.1 0.6 0.6 -0.1]';
%in this part, use matrix to present sequences
hd2 = [0 0 0 0 0 0 0 0 0 0 0 0 ...
    -0.1 0.6 0.6 -0.1 ...
    0 0 0 0 0 0 0 0 0 0 0 0]';
%set several zeros to avoid error
yd = zeros(16,1);

%represent the convolution equation
for k = 1:(length(x)+3)
    for n = 1:length(x)
        yd(k) = yd(k)+x(n)*hd2(k-n+13);
    end
end

n2 = [-pi:2*pi/100:pi];
hd3 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(hd)
    hd3(w) = hd3(w)+hd(k)*exp(-1*1i*n2(w)*(k-1));
    end
end
figure(2);
subplot(2,1,1);
plot(n2,abs(hd3));
subplot(2,1,2);
plot(n2,angle(hd3));

%-------------------------Q(e)-----------------------------
%initialize variables

x = [1 2 3 4 5 6 7 6 5 4 3 2 1]';
he = [0.2 0.5 0.2 0.1]';
%in this part, use matrix to present sequences
he2 = [0 0 0 0 0 0 0 0 0 0 0 0 ...
    0.2 0.5 0.2 0.1 ...
    0 0 0 0 0 0 0 0 0 0 0 0]';
%set several zeros to avoid error
ye = zeros(16,1);

%represent the convolution equation
for k = 1:(length(x)+3)
    for n = 1:length(x)
        ye(k) = ye(k)+x(n)*he2(k-n+13);
    end
end

n2 = [-pi:2*pi/100:pi];
he3 = zeros(length(n2),1);
for w = 1:length(n2)
    for k = 1:length(he)
    he3(w) = he3(w)+he(k)*exp(-1*1i*n2(w)*(k-1));
    end
end
figure(3);
subplot(2,1,1);
plot(n2,abs(he3));
subplot(2,1,2);
plot(n2,angle(he3));







