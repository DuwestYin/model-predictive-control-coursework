num = 1;
den = [2, 1];
sys = tf(num,den);
sys.InputDelay = 0.5;
step(sys);