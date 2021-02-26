function [x_curr, u_curr]=fu_function_runge_kutta(x_prev,u_prev,delta_t_disc)

x1_prev=x_prev(1);
x2_prev=x_prev(2);
x3_prev=x_prev(3);

x=[x1_prev;x2_prev;x3_prev]; u=u_prev;
f_11=function_radome_1(x);
f_21=function_radome_2(x);
f_31=function_radome_3(u);
rho=fu_function_rho_theta(x,u);
f_u1=fu_function_input(x,u,rho);

x=[x1_prev;x2_prev;x3_prev]+delta_t_disc/2*[f_11;f_21;f_31];
u=u_prev+delta_t_disc/2*f_u1;
f_12=function_radome_1(x);
f_22=function_radome_2(x);
f_32=function_radome_3(u);
rho=fu_function_rho_theta(x,u);
f_u2=fu_function_input(x,u,rho);

x=[x1_prev;x2_prev;x3_prev]+delta_t_disc/2*[f_12;f_22;f_32];
u=u_prev+delta_t_disc/2*f_u2;
f_13=function_radome_1(x);
f_23=function_radome_2(x);
f_33=function_radome_3(u);
rho=fu_function_rho_theta(x,u);
f_u3=fu_function_input(x,u,rho);

x=[x1_prev;x2_prev;x3_prev]+delta_t_disc*[f_13;f_23;f_33];
u=u_prev+delta_t_disc*f_u2;
f_14=function_radome_1(x);
f_24=function_radome_2(x);
f_34=function_radome_3(u);
rho=fu_function_rho_theta(x,u);
f_u4=fu_function_input(x,u,rho);

x1_curr=x1_prev+delta_t_disc/6*(f_11+2*f_12+2*f_13+f_14);
x2_curr=x2_prev+delta_t_disc/6*(f_21+2*f_22+2*f_23+f_24);
x3_curr=x3_prev+delta_t_disc/6*(f_31+2*f_32+2*f_33+f_34);
u_curr=u_prev+delta_t_disc/6*(f_u1+2*f_u2+2*f_u3+f_u4);

x_curr=[x1_curr;x2_curr;x3_curr];