clear all; clc;
v2_N = [-0.8393 0.4494 -0.3044]';
v1_N = [-0.1517 -0.9669 0.2050]';
v1_B = [0.8273 0.5541 -0.0920]';
v2_B = [-0.8285 0.5522 -0.0955]';
s1 = v1_B+v1_N;
s1_tilde = [0 -s1(3) s1(2); s1(3) 0 -s1(1); -s1(2) s1(1) 0]; 
s2 = v2_B+v2_N;
s2_tilde = [0 -s2(3) s2(2); s2(3) 0 -s2(1); -s2(2) s2(1) 0];
S = [s1_tilde; s2_tilde]
d = [v1_B-v1_N; v2_B-v2_N]
W = eye(6);
%d = d(:)''
qBar = inv(transpose(S)*W*S)*transpose(S)*W*d
BbarN = Gibbs2C(qBar)