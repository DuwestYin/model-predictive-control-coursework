clear; clc;close all;
num = 1;
den = [0.5, 1];
sys = tf(num,den);
t1 = 0:0.01:3;
[y, t1] = step(sys,t1);
figure,plot(t1,y);
n = length(t1);
n1 = abs(y-y(n))/y(n) > 0.05;
n2 = sum(n1);
N = 10;
T = n2*0.01/N;
delta = linspace(2,n2,N);
delta = floor(delta);
S = y(delta);
figure, plot(S);
I = eye(N);
Mss = zeros(1,N);
Mss(end) = 1;
Mss = [I;Mss];
Mss(1,:) = [];
Yk_1 = zeros(N,1);
Yk = zeros(N,1);
yk(1) = 0; yk(2) = 0;
t = 0:T:8;
t2 = 0:T:16;
u = sin(t);
u2 = sin(t2);
for i = 3:length(u);
    Yk = Mss*Yk_1 + S*(u(i-1)-u(i-2));
    yk(i) = Yk(1);
    Yk_1 = Yk;
end
figure, plot(t,yk,'r'); hold on;
[yy,t] = lsim(sys,u,t);
plot(t,yy,'b--');plot(t,u,'k');

y1k = 0;
k=0;

for L = N+1:length(u)
    for i = 1:N-1
        temp(i) = S(i)*(u(k+L-i)-u2(k+L-i-1));
    end
    temp(2) = sum(temp);
    y1k(L) = temp(2) + S(N)*u(k+L-N);
end
plot(t,y1k,'y')
