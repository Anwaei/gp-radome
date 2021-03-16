clear
clc
close all

N=4;
V_M=500;
tau=0.1;
T_alpha=1;

lambda_0 = deg2rad(5.7106);
R_0 = 10050;
gamma_M_0 = deg2rad(0);
A_M_0 = 0;
x_0 = [lambda_0, R_0, gamma_M_0];
rho_0 = fu_function_rho_theta(x_0,A_M_0);
theta_M_0 = gamma_M_0 + T_alpha*A_M_0/V_M;
lambda_r_0 = rho_0*(lambda_0-theta_M_0);
lambda_ME_0 = lambda_0 + lambda_r_0;
