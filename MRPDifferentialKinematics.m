clear all; clc;
h = 0.01; 
sigma = [0.4; 0.2; -.1 ];
i = [1 0 0; 0 1 0; 0 0 1];
for n= 1:1:4201
    ssm = [0 -sigma(3) sigma(2); sigma(3) 0 -sigma(1); -sigma(2) sigma(1) 0];
    sigma = sigma+0.25*((1-norm(sigma)^2)*i+2*ssm+2*sigma*sigma')*[sin(0.1*n*h); 0.01; cos(0.1*n*h)]*h*deg2rad(20);
    if norm(sigma) > 1
        sigma = (-sigma)/norm(sigma)^2;
    end
end

norm(sigma)