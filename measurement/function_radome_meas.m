function measurement=function_radome_meas(state,delta_t_disc,delta_t_cont,rho_theta_est,T_alpha,V_M,u_cont,i_disc)

index=1+(i_disc-1)*delta_t_disc/delta_t_cont;
index=round(index);

rho_theta_curr=rho_theta_est(index);
u_curr=u_cont(index);

x1=state(1);
x3=state(3);

alpha_M=T_alpha*u_curr/V_M;

measurement=(1+rho_theta_curr)*x1-rho_theta_curr*x3-rho_theta_curr*alpha_M;