%三容液位
clear; clc; close all;
Ts = 70; N = 15;    %采样时间和采样点数
h1 = 36; h2 = 18;   %液位设定值
p = 15; m = 10;     %预测周期和控制周期
dt = 0.1;           %仿真步长
k = 0:dt:2000;      %仿真时间
X0 = [0; 0; 0; 20; 0];  %微分方程初始条件
[t, Y1] = ode45('OdeFunc_1', k, X0);   %第一个输入
X0 = [0; 0; 0; 0; 50];
[t, Y2] = ode45('OdeFunc_1', k, X0);   %第二个输入
y11 = Y1(:,1);     %在不同输入下的输出值
y21 = Y1(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%作图
figure,
p3 = plot(y11,'k');
p3.LineWidth = 2;
hold on;
p4 = plot(y21,'k--');
p4.LineWidth = 2;
title('输入泵1的阶跃响应');
xlabel('时间t(s)');
ylabel('液位(h)/cm');
legend('h1','h2');
grid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y12 = Y2(:,1);
y22 = Y2(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%作图
figure,
p5 = plot(y12,'k');
p5.LineWidth = 2;
hold on;
p6 = plot(y22,'k--');
p6.LineWidth = 2;
title('输入泵2的阶跃响应');
xlabel('时间t(s)');
ylabel('液位(h)/cm');
legend('h1','h2');
grid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:15;
    S11(i) = y11(i*Ts/dt);   %对输出采样
    S21(i) = y21(i*Ts/dt);
    S12(i) = y12(i*Ts/dt);
    S22(i) = y22(i*Ts/dt);
end
S11 = S11/20; S21 = S21/20;     %构造S向量
S12 = S21/50; S22 = S22/50;
C = [eye(2),zeros(2,2*(N-1))];  %构造C矩阵
Cc = eye(2); Cm = eye(2);
Mss = zeros(2*N, 2*N);          %构造Mss矩阵
v=ones(2*(N-1),1);
Mss=diag(v,2);
Mss(end, end) = 1; Mss(end - 1, end - 1) = 1;
M = Mss;
for i = 1:N
    S(2*i-1:2*i, 1:2) = [S11(i),S12(i);S21(i),S22(i)];  %构造S矩阵
end
Su = zeros(2*p, 2*m);
for i = 1:m
    Su(i*2-1:end,i*2-1:i*2) = S(1:2*p-2*(i-1),1:2);     %构造Su矩阵
end
Ty = eye(2*p);    %加权系数矩阵
Tu = 0.35*eye(2*m);
%Kmpc增益
Kmpc = [eye(2),zeros(2,2*(m-1))]*pinv(Su'*(Ty'*Ty)*Su + Tu'*Tu)*Su'*(Ty'*Ty);
Y_ = zeros(2*N,1);   %观测的状态Y
u = 0;    %控制输入
du = zeros(2,1);     %控制增量
ym = zeros(2,1);     %当前测量值
Kf = zeros(2*N,2);   %观测器反馈系数
v = 0.2*eye(2);
Kf = repmat(v,15,1);   
Rk = zeros(2*N,1);   %设定值
Rk(1:2:2*N-1) = 36;     
Rk(2:2:2*N) = 18;
Y = zeros(1,5);      %保存微分方程解的最后一组值
y1 = zeros(50,0);    %用于显示仿真结果
y2 = y1;
for k = 1:50
    ym(1) = Y(end,1); %测量当前状态
    ym(2) = Y(end,2);
    h3 = Y(end,3);
    y1(k) = ym(1);    %保存以供接下来显示
    y2(k) = ym(2);
    Y_ = Mss*Y_ + S*du + Kf*(ym - C*Y_);  %观测器
    Ep = Rk - M*Y_;   %误差
    du = Kmpc*Ep;     %计算控制增量
    u = u + du;
    X0 = [ym(1), ym(2), h3, u(1), u(2)];  %更新初始条件
    [t, Y] = ode45('OdeFunc_1', 0:0.01:70, X0);  %计算微分方程
end
figure,
p1 = plot(y1,'k');
p1.LineWidth = 2;
hold on;
p2 = plot(y2,'k--');
p2.LineWidth = 2;
title('无约束三容液位控制');
xlabel('控制步数(k)（Ts = 70s）');
ylabel('液位(h)/cm');
legend('h1','h2');
grid;



