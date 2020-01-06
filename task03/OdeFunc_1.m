function dx = OdeFunc_1(t, x)
% x : h1 h2 h3 q1 q1
dx = zeros(5,1);
az1 = 0.4637; az2 = 0.6751; az3 = 0.4680;  %系统结构参数
g = 980; S = 154; Sn = 0.5;
Q13 = az1*Sn*sign(x(1)-x(3))*sqrt(2*g*abs(x(1)-x(3)));  %中间变量
Q32 = az3*Sn*sign(x(3)-x(2))*sqrt(2*g*abs(x(3)-x(2)));
Q20 = az2*Sn*sqrt(2*g*x(2));
dx(1) = (1/S)*(x(4) - Q13);       %微分方程
dx(2) = (1/S)*(x(5) + Q32 - Q20);
dx(3) = (1/S)*(Q13 - Q32);
%q1,q2的微分默认为0
end