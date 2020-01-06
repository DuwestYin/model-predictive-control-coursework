function Kf = get_Kf(Mss, C)
A = Mss';
B = C';
Qc = ctrb(A,B);
P = Qc(:,1:30);
P_ = pinv(P);
A_ = P_*A*P
end