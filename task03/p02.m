%带约束的三容液位
clear; clc; close all;
Ts = 70; N = 15;
du_max = 5;   %控制增量约束
du_min = -5;
h1 = 36; h2 = 18; 
p = 15; m = 10; 
dt = 0.1;
k = 0:dt:2000;
X0 = [0; 0; 0; 20; 0];
[t, Y1] = ode45('OdeFunc_1', k, X0);
X0 = [0; 0; 0; 0; 50];
[t, Y2] = ode45('OdeFunc_1', k, X0);
y11 = Y1(:,1);
y21 = Y1(:,2);
y12 = Y2(:,1);
y22 = Y2(:,2);
for i = 1:15;
    S11(i) = y11(i*Ts/dt);
    S21(i) = y21(i*Ts/dt);
    S12(i) = y12(i*Ts/dt);
    S22(i) = y22(i*Ts/dt);
end
S11 = S11/20; S21 = S21/20;
S12 = S21/50; S22 = S22/50;
C = [eye(2),zeros(2,2*(N-1))];
Cc = eye(2); Cm = eye(2);
Mss = zeros(2*N, 2*N);
v=ones(2*(N-1),1);
Mss=diag(v,2);
Mss(end, end) = 1; Mss(end - 1, end - 1) = 1;
M = Mss;
for i = 1:N
    S(2*i-1:2*i, 1:2) = [S11(i),S12(i);S21(i),S22(i)];
end
Su = zeros(2*p, 2*m);
for i = 1:m
    Su(i*2-1:end,i*2-1:i*2) = S(1:2*p-2*(i-1),1:2);
end
Ty = eye(2*p);
Tu = 0.35*eye(2*m);
Kmpc = [eye(2),zeros(2,2*(m-1))]*pinv(Su'*(Ty'*Ty)*Su + Tu'*Tu)*Su'*(Ty'*Ty);
Y_ = zeros(2*N,1);
dU = zeros(2*m,1);
u = 0;
du = zeros(2,1);
ym = zeros(2,1);
Kf = zeros(2*N,2);
v = 0.2*eye(2);
Kf = repmat(v,15,1);
Rk = zeros(2*N,1);
Rk(1:2:2*N-1) = 36;
Rk(2:2:2*N) = 18;
Y = zeros(1,5);
y1 = zeros(50,0);
y2 = y1;
H = Su'*(Ty'*Ty)*Su + Tu'*Tu;   %QP问题参数矩阵
T = eye(2*m);
A = [T; -T];
b = zeros(4*m,1);
b(1:2*m) = du_max;
b(2*m+1:end) = -du_min;
for k = 1:50
    ym(1) = Y(end,1);
    ym(2) = Y(end,2);
    h3 = Y(end,3);
    y1(k) = ym(1);
    y2(k) = ym(2);
    Y_ = Mss*Y_ + S*du + Kf*(ym - C*Y_);
    Ep = Rk - M*Y_;
    f = 2*Su'*(Ty'*Ty)*Ep;
    f = -real(f);
    x = quadprog(H,f,A,b);   %解QP得到一串输入
    du(1) = x(1);            %取第一组输入（多输入被排成了一列）
    du(2) = x(2);
    u = u + du;
    X0 = [ym(1), ym(2), h3, u(1), u(2)];
    [t, Y] = ode45('OdeFunc_1', 0:0.01:70, X0);
end
p1 = plot(y1,'k');
p1.LineWidth = 2;
hold on;
p2 = plot(y2,'k--');
p2.LineWidth = 2;
title('有约束三容液位控制');
xlabel('控制步数(k)（Ts = 70s）');
ylabel('液位(h)/cm');
legend('h1','h2');
grid;

