function y = fu_function_rho_theta(x,u)
%%
lambda = x(1);
gamma_M = x(3);
gamma_dot = function_radome_3(u);
T_alpha = 1;
alpha_M = T_alpha*gamma_dot;
theta_M = alpha_M + gamma_M;
lambda_s = lambda - theta_M;
%%
%%
r=-0.5;
zeta=4;
p=1;

if lambda_s<0
    a=exp(zeta*lambda_s);
    b=zeta*sin(2*pi*lambda_s/p);
    c=2*pi/p*cos(2*pi*lambda_s/p);
else
    a=exp(-zeta*lambda_s);
    b=-zeta*sin(2*pi*lambda_s/p);
    c=2*pi/p*cos(2*pi*lambda_s/p);
end

y=r*a*(b+c);
y=y/rad2deg(1);

%%
% y=0.025;
% y=0.05;
% y=-0.025;
end