%����Һλ
clear; clc; close all;
Ts = 70; N = 15;    %����ʱ��Ͳ�������
h1 = 36; h2 = 18;   %Һλ�趨ֵ
p = 15; m = 10;     %Ԥ�����ںͿ�������
dt = 0.1;           %���沽��
k = 0:dt:2000;      %����ʱ��
X0 = [0; 0; 0; 20; 0];  %΢�ַ��̳�ʼ����
[t, Y1] = ode45('OdeFunc_1', k, X0);   %��һ������
X0 = [0; 0; 0; 0; 50];
[t, Y2] = ode45('OdeFunc_1', k, X0);   %�ڶ�������
y11 = Y1(:,1);     %�ڲ�ͬ�����µ����ֵ
y21 = Y1(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ͼ
figure,
p3 = plot(y11,'k');
p3.LineWidth = 2;
hold on;
p4 = plot(y21,'k--');
p4.LineWidth = 2;
title('�����1�Ľ�Ծ��Ӧ');
xlabel('ʱ��t(s)');
ylabel('Һλ(h)/cm');
legend('h1','h2');
grid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y12 = Y2(:,1);
y22 = Y2(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ͼ
figure,
p5 = plot(y12,'k');
p5.LineWidth = 2;
hold on;
p6 = plot(y22,'k--');
p6.LineWidth = 2;
title('�����2�Ľ�Ծ��Ӧ');
xlabel('ʱ��t(s)');
ylabel('Һλ(h)/cm');
legend('h1','h2');
grid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:15;
    S11(i) = y11(i*Ts/dt);   %���������
    S21(i) = y21(i*Ts/dt);
    S12(i) = y12(i*Ts/dt);
    S22(i) = y22(i*Ts/dt);
end
S11 = S11/20; S21 = S21/20;     %����S����
S12 = S21/50; S22 = S22/50;
C = [eye(2),zeros(2,2*(N-1))];  %����C����
Cc = eye(2); Cm = eye(2);
Mss = zeros(2*N, 2*N);          %����Mss����
v=ones(2*(N-1),1);
Mss=diag(v,2);
Mss(end, end) = 1; Mss(end - 1, end - 1) = 1;
M = Mss;
for i = 1:N
    S(2*i-1:2*i, 1:2) = [S11(i),S12(i);S21(i),S22(i)];  %����S����
end
Su = zeros(2*p, 2*m);
for i = 1:m
    Su(i*2-1:end,i*2-1:i*2) = S(1:2*p-2*(i-1),1:2);     %����Su����
end
Ty = eye(2*p);    %��Ȩϵ������
Tu = 0.35*eye(2*m);
%Kmpc����
Kmpc = [eye(2),zeros(2,2*(m-1))]*pinv(Su'*(Ty'*Ty)*Su + Tu'*Tu)*Su'*(Ty'*Ty);
Y_ = zeros(2*N,1);   %�۲��״̬Y
u = 0;    %��������
du = zeros(2,1);     %��������
ym = zeros(2,1);     %��ǰ����ֵ
Kf = zeros(2*N,2);   %�۲�������ϵ��
v = 0.2*eye(2);
Kf = repmat(v,15,1);   
Rk = zeros(2*N,1);   %�趨ֵ
Rk(1:2:2*N-1) = 36;     
Rk(2:2:2*N) = 18;
Y = zeros(1,5);      %����΢�ַ��̽�����һ��ֵ
y1 = zeros(50,0);    %������ʾ������
y2 = y1;
for k = 1:50
    ym(1) = Y(end,1); %������ǰ״̬
    ym(2) = Y(end,2);
    h3 = Y(end,3);
    y1(k) = ym(1);    %�����Թ���������ʾ
    y2(k) = ym(2);
    Y_ = Mss*Y_ + S*du + Kf*(ym - C*Y_);  %�۲���
    Ep = Rk - M*Y_;   %���
    du = Kmpc*Ep;     %�����������
    u = u + du;
    X0 = [ym(1), ym(2), h3, u(1), u(2)];  %���³�ʼ����
    [t, Y] = ode45('OdeFunc_1', 0:0.01:70, X0);  %����΢�ַ���
end
figure,
p1 = plot(y1,'k');
p1.LineWidth = 2;
hold on;
p2 = plot(y2,'k--');
p2.LineWidth = 2;
title('��Լ������Һλ����');
xlabel('���Ʋ���(k)��Ts = 70s��');
ylabel('Һλ(h)/cm');
legend('h1','h2');
grid;



