function Am_dot=fu_function_input(x,u,rho)

N=4;
V_M=500;
tau=0.1;
T_alpha=1;

% lambda=x(1);
% R=x(2);
% gamma=x(3);

lambda_dot = function_radome_1(x);
R_dot = abs(function_radome_2(x));
% R_dot = -function_radome_2(x);

Am_dot = (N*V_M*(1+rho)*R_dot*lambda_dot)/(N*T_alpha*rho*R_dot+tau*V_M) ...
    - (N*rho*R_dot+V_M)/(N*T_alpha*rho*R_dot+tau*V_M)*u;