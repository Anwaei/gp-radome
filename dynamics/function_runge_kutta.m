function x_curr=function_runge_kutta(x_prev,delta_t_disc,delta_t_cont,u_cont,i_disc)

x1_prev=x_prev(1);
x2_prev=x_prev(2);
x3_prev=x_prev(3);

x=[x1_prev;x2_prev;x3_prev];
f_11=function_radome_1(x);
f_21=function_radome_2(x);
%
index=1+(i_disc-2)*delta_t_disc/delta_t_cont;
index=round(index);
%
u_disc=u_cont(index);
f_31=function_radome_3(u_disc);

x=[x1_prev;x2_prev;x3_prev]+delta_t_disc/2*[f_11;f_21;f_31];
f_12=function_radome_1(x);
f_22=function_radome_2(x);
%
index=1+(i_disc-2)*delta_t_disc/delta_t_cont+0.5*delta_t_disc/delta_t_cont;
index=round(index);
%
u_disc=u_cont(index);
f_32=function_radome_3(u_disc);

x=[x1_prev;x2_prev;x3_prev]+delta_t_disc/2*[f_12;f_22;f_32];
f_13=function_radome_1(x);
f_23=function_radome_2(x);
%
index=1+(i_disc-2)*delta_t_disc/delta_t_cont+0.5*delta_t_disc/delta_t_cont;
index=round(index);
%
u_disc=u_cont(index);
f_33=function_radome_3(u_disc);

x=[x1_prev;x2_prev;x3_prev]+delta_t_disc*[f_13;f_23;f_33];
f_14=function_radome_1(x);
f_24=function_radome_2(x);
%
index=1+(i_disc-1)*delta_t_disc/delta_t_cont;
index=round(index);
%
u_disc=u_cont(index);
f_34=function_radome_3(u_disc);

x1_curr=x1_prev+delta_t_disc/6*(f_11+2*f_12+2*f_13+f_14);
x2_curr=x2_prev+delta_t_disc/6*(f_21+2*f_22+2*f_23+f_24);
x3_curr=x3_prev+delta_t_disc/6*(f_31+2*f_32+2*f_33+f_34);

x_curr=[x1_curr;x2_curr;x3_curr];