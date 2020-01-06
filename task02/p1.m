clear;clc;close all;
T = 0.2; N = 8; m =3; p =5;
Suu = [0.3297, 0.5507, 0.6988, 0.7981, 0.8647, 0.9093, 0.9392, 0.9592]';
Sdd = ones(N,1);
KI = ones(N,1);
I = eye(N);
Ty = eye(5);
diag1 = 0.2*ones(1,m);
Tu = diag(diag1);
I = eye(N);
Mss = zeros(1,N);
Mss(end) = 1;
Mss = [I;Mss];
Mss(1,:) = [];
C = zeros(1,N);
C(1) = 1;
M = []; 
for i = 1:p
    M = [M;C*Mss^i];
end
Su = zeros(p,m);
for i = 1:m
    Su(i:end,i) = Suu(1:p-i+1);
end
Sd = Sdd(1:p);
Kmpc = [1, 0 , 0]*pinv(Su'*(Ty'*Ty)*Su + Tu'*Tu)*Su'*(Ty'*Ty);
a = (I-KI*C)*Mss;
b = (I-KI*C)*Suu;
c = (I-KI*C)*Sdd;