clear all; clc
w1 =1; w2 = 1;
lambda = w1+w2;

v2_N = [-0.8393 0.4494 -0.3044]';
v1_N = [-0.1517 -0.9669 0.2050]';
v1_B = [0.8273 0.5541 -0.0920]';
v2_B = [-0.8285 0.5522 -0.0955]';

B = w1*v1_B*v1_N'+w2*v2_B*v2_N';
sigma = trace(B);
S = B + B';
Z = [ B(2,3) - B(3,2) ; ...
      B(3,1) - B(1,3) ; ...
      B(1,2) - B(2,1) ];
K = [ sigma Z' ; ...
      Z (S-sigma*eye(3))];
f = @(s) det(K-s*eye(4));
fd = @(s) f(s)*trace(inv(K-s*eye(4))*(-eye(4)));
lambda_optimize = QuestComp(f,fd,lambda,10)
q = inv((lambda_optimize+sigma)*eye(3)-S)*Z
qnorm = q/norm(q)
beta = [1; q]*(1/sqrt(1+q'*q))
BbarN = Gibbs2C(q)
function x = QuestComp(f,fd,x0,n)
    x = x0;
    for i = 1:n
        x = x - f(x)/fd(x)
    end
end