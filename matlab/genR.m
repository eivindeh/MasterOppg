N = 6;
M = 4;
t1 = [zeros(N-1,1),eye(N-1);zeros(1,N)];
t1(6,6) = 0;
t1(1,6) = 0;
A1 = t1+t1';

A2 = diag([1,0,1,0,1,0]);
A3 = diag([0,1,0,1,0,1]);
A4 = zeros(N);
R = [A1,A2,A4,A4;A2,A1,A3,A4;A4,A3,A1,A2;A4,A4,A2,A1];
