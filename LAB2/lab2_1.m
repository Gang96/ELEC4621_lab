clc;
clear;
%% -----------Q(a)-------------
%a = 2cos(w0)
%b = 1

%% -----------Q(b)-------------

% variables initialize
w = pi/4;
a = 2*cos(w);
b = 1;
N = (1:20)';
y = zeros(20,1);

% output calculation
y(1) = 1; %in matlab, y(0) is meaningless, so the series begin with 1
for n = 2:length(y)
    if n == 2
        y(n) = a*y(n-1);
    else
        y(n) = a*y(n-1)-b*y(n-2);
    end
end

%plot the output
figure(1);
stem(N,y);

%% -----------Q(c)-------------
y(1) = sqrt(1/2);
for n = 2:length(y)
    if n == 2
        y(n) = a*y(n-1);
    else
        y(n) = a*y(n-1)-b*y(n-2);
    end
end
figure(2);
stem(N,y);

% x[0] = sqrt(1/2);
% x[-1] = 0;

