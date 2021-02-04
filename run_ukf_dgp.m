clc
close all

%% read CONTinuous data from simulink
t=out.tout;
R=out.yout{1}.Values.Data;
lambda=out.yout{2}.Values.Data;
dot_gamma_M=out.yout{4}.Values.Data;
A_M=out.yout{6}.Values.Data;
gamma_M=out.yout{7}.Values.Data;
rho_theta=out.yout{8}.Values.Data;
alpha_M=out.yout{9}.Values.Data;
theta_s=out.yout{10}.Values.Data;
theta_M=out.yout{11}.Values.Data;
lambda_r=out.yout{12}.Values.Data;
lambda_ME=out.yout{13}.Values.Data;

%% dynamics
% (1) read CONTinuous data from simulink
delta_t_cont=t(2)-t(1);
x1_cont=lambda;
x2_cont=R;
x3_cont=gamma_M;
u_cont=A_M;

% (2) calculate DISCrete data with 4th order Runge Kutta rule
delta_t_disc=0.001;
% delta_t_disc=0.005;
% delta_t_disc=0.01;
% delta_t_disc=0.02;
% delta_t_disc=0.5;
T=t(end)/delta_t_disc;
x1_disc=zeros(T+1,1);
x2_disc=zeros(T+1,1);
x3_disc=zeros(T+1,1);
for i=1:length(x1_disc)
    if i==1
        x1_disc(i)=x1_cont(i);
        x2_disc(i)=x2_cont(i);
        x3_disc(i)=x3_cont(i);
    else
        x_prev=[x1_disc(i-1);x2_disc(i-1);x3_disc(i-1)];
        x_curr=function_runge_kutta(x_prev,delta_t_disc,delta_t_cont,u_cont,i);
        % archive
        x1_disc(i)=x_curr(1);
        x2_disc(i)=x_curr(2);
        x3_disc(i)=x_curr(3);
    end
end

% (3) testing, plot figure (CONTinuous VS DISCrete)
figure
subplot(2,3,1)
plot(t,rad2deg(x1_cont),0:delta_t_disc:t(end),rad2deg(x1_disc))
xlabel('t (s)')
ylabel('\lambda (deg)')
legend('continuous','discrete')
title('LOS angle')
grid on
axis square
subplot(2,3,2)
plot(t,x2_cont,0:delta_t_disc:t(end),x2_disc)
xlabel('t (s)')
ylabel('R (m)')
legend('continuous','discrete')
title('relative distance between missile & target')
grid on
axis square
subplot(2,3,3)
plot(t,rad2deg(x3_cont),0:delta_t_disc:t(end),rad2deg(x3_disc))
xlabel('t (s)')
ylabel('\gamma_M (deg)')
legend('continuous','discrete')
title('flight path angle')
grid on
axis square
subplot(2,3,4)
error_lambda=zeros(T+1,1);
for i=1:T+1
    index=1+(i-1)*delta_t_disc/delta_t_cont;
    index=round(index);
    error_lambda(i)=x1_cont(index)-x1_disc(i);
end
plot(0:delta_t_disc:t(end),rad2deg(error_lambda))
xlabel('t (s)')
ylabel('\Delta \lambda (deg)')
title('discretization error of LOS angle')
grid on
axis square
subplot(2,3,5)
error_R=zeros(T+1,1);
for i=1:T+1
    index=1+(i-1)*delta_t_disc/delta_t_cont;
    index=round(index);
    error_R(i)=x2_cont(index)-x2_disc(i);
end
plot(0:delta_t_disc:t(end),error_R)
xlabel('t (s)')
ylabel('\Delta R (m)')
title('discretization error of relative distance')
grid on
axis square
subplot(2,3,6)
error_gamma_M=zeros(T+1,1);
for i=1:T+1
    index=1+(i-1)*delta_t_disc/delta_t_cont;
    index=round(index);
    error_gamma_M(i)=x3_cont(index)-x3_disc(i);
end
plot(0:delta_t_disc:t(end),rad2deg(error_gamma_M))
xlabel('t (s)')
ylabel('\Delta \gamma_M (deg)')
title('discretization error of flight path angle')
grid on
axis square

%% measurement
% continuous
z_cont=lambda_ME;

% discrete
z_disc=zeros(T+1,1);
for i=1:T+1
    index=1+(i-1)*delta_t_disc/delta_t_cont;
    index=round(index);
    z_disc(i)=lambda_ME(index);
end

% filter
% (1) set variance of measurement noise
rho_theta_est=zeros(T+1,1);
for i=1:length(rho_theta_est)
    rho_theta_est(i)=rho_theta(i);
%     rho_theta_est(i)=2*rho_theta(i);
end

% R_filter=(deg2rad(5*1e-1))^2;
% R_filter=(deg2rad(1e-1))^2;
% R_filter=(deg2rad(5*1e-2))^2;
% R_filter=(deg2rad(1e-2))^2;
% R_filter=(deg2rad(1e-3))^2;
R_filter=(deg2rad(5*1e-3))^2;
% R_filter=(deg2rad(1e-4))^2;
% R_filter=0;

% (2) calculate measurement from true state

z_est=zeros(T+1,1);
for i=1:length(z_est)
    state=zeros(3,1);
    state(1)=x1_disc(i);
    state(2)=x2_disc(i);
    state(3)=x3_cont(i);
    z_est(i)=function_radome_meas(state,delta_t_disc,delta_t_cont,rho_theta_est,T_alpha,V_M,u_cont,i);
end

% (3) add measurement noise
for i=1:length(z_est)
    v=normrnd(0,sqrt(R_filter));
    z_est(i)=z_est(i)+v;
end

% testing, plot figure LOS angle measurement(CONTinuous VS DISCrete VS ESTimate)
figure
subplot(2,1,1)
plot(t,rad2deg(z_cont))
hold on
plot(0:delta_t_disc:t(end),rad2deg(z_disc))
plot(0:delta_t_disc:t(end),rad2deg(z_est))
xlabel('t (s)')
ylabel('\lambda_{ME} (deg)')
title('LOS angle measurement')
grid on
legend('continuous','discrete','estimate')
axis square
subplot(2,1,2)
plot(0:delta_t_disc:t(end),rad2deg(z_disc-z_est))
xlabel('t (s)')
ylabel('\Delta \lambda_{ME} (deg)')
title('modeling error of LOS angle measurement')
grid on
axis square

%% filter
% (1) set covariance of process noise

% Q_1=(deg2rad(1e-2))^2;
% Q_1=(deg2rad(5*1e-3))^2;
Q_1=(deg2rad(1e-3))^2;
% Q_1=(deg2rad(1e-4))^2;

% Q_2=1e1;
% Q_2=1e2;
% Q_2=2*1e2;
% Q_2=3*1e2;
% Q_2=4*1e2;
Q_2=5*1e2;
% Q_2=1e3;
% Q_2=5*1e3;
% Q_2=1e4;
% Q_2=1e5;
% Q_2=1e6;

% Q_3=(deg2rad(10))^2;
% Q_3=(deg2rad(2))^2;
% Q_3=(deg2rad(1))^2;
% Q_3=(deg2rad(5*1e-1))^2;

% Q_3=(deg2rad(1e-1))^2;
% Q_3=(deg2rad(5*1e-2))^2;

% Q_3=(deg2rad(1e-2))^2;
% Q_3=(deg2rad(5*1e-3))^2;
Q_3=(deg2rad(1e-3))^2;
% Q_3=(deg2rad(5*1e-4))^2;
% Q_3=(deg2rad(1e-4))^2;

Q_filter=diag([Q_1,Q_2,Q_3]);

% (2) initial
x_est=zeros(T+1,3);
P_est=zeros(T+1,3,3);
x_0=zeros(3,1);
x_0(1)=deg2rad(7.7106);
x_0(2)=12050;
x_0(3)=deg2rad(2);
P_0=zeros(3,3);
P_0(1,1)=(deg2rad(1))^2;
P_0(2,2)=(1e3)^2;
P_0(3,3)=(deg2rad(1))^2;
% P_0(1,1)=(deg2rad(2))^2;
% P_0(2,2)=(2*1e3)^2;
% P_0(3,3)=(deg2rad(2))^2;

% (3) main loop
for i=1:length(x_est)
    i
    if i==1
        x_est(i,:)=x_0';
        P_est(i,:,:)=P_0;
    else
        x_prev=(x_est(i-1,:))';
        P_prev=shiftdim(P_est(i-1,:,:));
        
  % MATLAB built in     
%         f=@(state)function_runge_kutta(state,delta_t_disc,delta_t_cont,u_cont,i);
%         h=@(state)function_radome_meas(state,delta_t_disc,delta_t_cont,rho_theta,T_alpha,V_M,u_cont,i);
%         filter=trackingUKF(f,h,'State',x_prev,'StateCovariance',P_prev,'ProcessNoise',Q_filter,'MeasurementNoise',R_filter);  
%         [x_pred,P_pred]=predict(filter);
%         [x_corr,P_corr]=correct(filter,z_disc(i));
   % MATLAB built in
        [x_pred,P_pred,x_corr,P_corr]=function_ukf(x_prev,P_prev,z_disc(i),Q_filter,R_filter,delta_t_disc,delta_t_cont,rho_theta_est,T_alpha,V_M,u_cont,i);       
        x_est(i,:)=x_corr';
        P_est(i,:,:)=P_corr;
    end
end

% (4) testing, plot figure (DISCrete VS ESTimate)
filter_error=zeros(T+1,size(x_est,2));
for i=1:length(filter_error)
    filter_error(i,1)=x1_disc(i)-x_est(i,1);
    filter_error(i,2)=x2_disc(i)-x_est(i,2);
    filter_error(i,3)=x3_disc(i)-x_est(i,3);
end
figure
subplot(2,3,1)
plot(t,rad2deg(x1_cont),0:delta_t_disc:t(end),rad2deg(x_est(:,1)))
xlabel('t (s)')
ylabel('\lambda (deg)')
legend('discrete','estimate')
title('LOS angle')
grid on
axis square
subplot(2,3,2)
plot(0:delta_t_disc:t(end),x2_disc,0:delta_t_disc:t(end),x_est(:,2))
xlabel('t (s)')
ylabel('R (m)')
legend('discrete','estimate')
title('relative distance between missile & target')
grid on
axis square
subplot(2,3,3)
plot(0:delta_t_disc:t(end),rad2deg(x3_disc),0:delta_t_disc:t(end),rad2deg(x_est(:,3)))
xlabel('t (s)')
ylabel('\gamma_M (deg)')
legend('discrete','estimate')
title('flight path angle')
grid on
axis square
subplot(2,3,4)
plot(0:delta_t_disc:t(end),rad2deg(filter_error(:,1)))
xlabel('t (s)')
ylabel('\Delta \lambda (deg)')
title('estimation error: \lambda')
grid on
axis square
subplot(2,3,5)
plot(0:delta_t_disc:t(end),filter_error(:,2))
xlabel('t (s)')
ylabel('\Delta R (m)')
title('estimation error: R')
grid on
axis square
subplot(2,3,6)
plot(0:delta_t_disc:t(end),rad2deg(filter_error(:,3)))
xlabel('t (s)')
ylabel('\Delta \gamma_M (deg)')
title('estimation error: \gamma_M')
grid on
axis square

%% Gaussian process regression
% (1) calculate look angle
% (1.1) calculate alpha_M_disc
alpha_M_disc=zeros(T+1,1);
for i=1:length(alpha_M_disc)
    index=1+(i-1)*delta_t_disc/delta_t_cont;
    index=round(index);
    alpha_M_disc(i)=alpha_M(index);
end
% (1.2) calculate look angle (continuous, discrete, estimate)
% (1.2.1) discrete
theta_s_disc=zeros(T+1,1);
for i=1:length(theta_s_disc)
    theta_s_disc(i)=x1_disc(i)-(alpha_M_disc(i)+x3_disc(i));
end
% (1.2.2) estimate
theta_M_est=zeros(T+1,1);
for i=1:length(theta_M_est)
    theta_M_est(i)=x_est(i,3)+alpha_M_disc(i);
end
theta_s_est=zeros(T+1,1);
for i=1:length(theta_s_est)
    theta_s_est(i)=x_est(i,1)-theta_M_est(i);
end
% (1.3) testing, plot look angle (CONTinuous VS DISCrete VS ESTimate)
figure
subplot(3,3,1)
plot(t,rad2deg(theta_s))
hold on
plot(0:delta_t_disc:t(end),rad2deg(theta_s_disc))
plot(0:delta_t_disc:t(end),rad2deg(theta_s_est))
title('look angle')
xlabel('t (s)')
ylabel('\theta_s (deg)')
legend('continuous','discrete','estimate')
grid on
axis square

% (2) calculate radome error (continuous, discrete, estimate)
% (2.1) discrete
lambda_r_disc=zeros(T+1,1);
for i=1:length(lambda_r_disc)
    index=1+(i-1)*delta_t_disc/delta_t_cont;
    index=round(index);
    lambda_r_disc(i)=lambda_ME(index)-x1_disc(i);
end
% (2.2) estimate
lambda_r_est=zeros(T+1,1);
for i=1:length(lambda_r_est)
    lambda_r_est(i)=z_disc(i)-x_est(i,1);
end
% (2.3) testing, plot radome error (CONTinuous VS DISCrete VS ESTimate)
subplot(3,3,2)
plot(t,rad2deg(lambda_r))
hold on
plot(0:delta_t_disc:t(end),rad2deg(lambda_r_disc))
plot(0:delta_t_disc:t(end),rad2deg(lambda_r_est))
title('radome error')
xlabel('t (s)')
ylabel('\lambda_r (deg)')
legend('continuous','discrete','estimate')
grid on
axis square
% (2.4) testing, plot radome error VS look angle
subplot(3,3,3)
index_m=2/delta_t_disc:100:length(theta_s);
index_s=2/delta_t_disc:75:length(theta_s);
plot(rad2deg(theta_s(index_m)),rad2deg(lambda_r(index_m)))
hold on
index_m_disc=2/delta_t_disc:100:length(theta_s_disc);
index_s_disc=2/delta_t_disc:75:length(theta_s_disc);
plot(rad2deg(theta_s_disc(index_m_disc)),rad2deg(lambda_r_disc(index_m_disc)))
index_m_est=2/delta_t_disc:100:length(theta_s_est);
index_s_est=2/delta_t_disc:75:length(theta_s_est);
plot(rad2deg(theta_s_est(index_m_est)),rad2deg(lambda_r_est(index_m_est)))
title('radome error VS look angle')
xlabel('\theta_s (deg)')
ylabel('\lambda_r (deg)')
legend('continuous','discrete','estimate')
grid on
axis square

% (2.4) calculate radome slope via Gaussian process (continuous, discrete, estimate)
% (2.4.1) continuous

Xm=rad2deg(theta_s(index_m));
Xm=Xm';
fmh=rad2deg(lambda_r(index_m));

% tuning hyperparameters
lf=1;
lx=1;

sfm=0.1;
% sfm=10;

% hyp=function_gp_hyperparam(Xm,fmh);
% lx = sqrt(hyp(1));
% lf = sqrt(hyp(2));
% sfm = sqrt(hyp(3));

Xs=rad2deg(theta_s(index_s));
Xs=Xs';

% GPR
[mPost,SPost,mdPost,SdPost]=function_gp_dgp(Xs,Xm,fmh,lf,lx,sfm);
sPost=sqrt(diag(SPost)); % These are the posterior standard deviations.
sdPost=sqrt(diag(SdPost)); % These are the posterior standard deviations of the derivative.
% GP plot
% [h_1,h_2]=function_gp_plot(Xs,Xm,fmh,mPost,sPost,mdPost,sdPost);
subplot(3,3,4)
hold on;
grid on;
patch([Xs, fliplr(Xs)],[mPost-2*sPost; flipud(mPost+2*sPost)], 1, 'FaceColor', [0.9,0.9,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
patch([Xs, fliplr(Xs)],[mPost-sPost; flipud(mPost+sPost)], 1, 'FaceColor', [0.8,0.8,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(Xs, mPost, 'b-', 'LineWidth', 1); % We plot the mean line.
plot(Xm, fmh, 'ro'); % We plot the measurement points.
title('radome error VS look angle (GPR continuous)')
xlabel('\theta_s (deg)')
ylabel('\lambda_r (deg)')
axis square
subplot(3,3,7)
hold on;
grid on;
patch([Xs, fliplr(Xs)],[mdPost-2*sdPost; flipud(mdPost+2*sdPost)], 1, 'FaceColor', [0.9,0.9,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
patch([Xs, fliplr(Xs)],[mdPost-sdPost; flipud(mdPost+sdPost)], 1, 'FaceColor', [0.8,0.8,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(Xs, mdPost, 'b-', 'LineWidth', 1); % We plot the mean line.
title('radome slope VS look angle (GPR continuous)')
xlabel('\theta_s (deg)')
ylabel('\rho_{\theta} (deg/deg)')
axis square

% (2.4.2) discrete


Xm_disc=rad2deg(theta_s_disc(index_m_disc));
Xm_disc=Xm_disc';
fmh_disc=rad2deg(lambda_r_disc(index_m_disc));

% tuning hyperparameters
lf_disc=1;
lx_disc=1;
sfm_disc=0.1;
% hyp=function_gp_hyperparam(Xm_disc,fmh_disc);
% lx_disc = sqrt(hyp(1));
% lf_disc = sqrt(hyp(2));
% sfm_disc = sqrt(hyp(3));

Xs_disc=rad2deg(theta_s_disc(index_s_disc));
Xs_disc=Xs_disc';

% GPR
[mPost_disc,SPost_disc,mdPost_disc,SdPost_disc]=function_gp_dgp(Xs_disc,Xm_disc,fmh_disc,lf_disc,lx_disc,sfm_disc);
sPost_disc=sqrt(diag(SPost_disc)); % These are the posterior standard deviations.
sdPost_disc=sqrt(diag(SdPost_disc)); % These are the posterior standard deviations of the derivative.
% GP plot
% [h_3,h_4]=function_gp_plot(Xs,Xm,fmh,mPost_disc,sPost_disc,mdPost_disc,sdPost_disc);
subplot(3,3,5)
hold on;
grid on;
patch([Xs, fliplr(Xs)],[mPost_disc-2*sPost_disc; flipud(mPost_disc+2*sPost_disc)], 1, 'FaceColor', [0.9,0.9,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
patch([Xs, fliplr(Xs)],[mPost_disc-sPost_disc; flipud(mPost_disc+sPost_disc)], 1, 'FaceColor', [0.8,0.8,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(Xs, mPost_disc, 'b-', 'LineWidth', 1); % We plot the mean line.
plot(Xm, fmh, 'ro'); % We plot the measurement points.
title('radome error VS look angle (GPR discrete)')
xlabel('\theta_s (deg)')
ylabel('\lambda_r (deg)')
axis square
subplot(3,3,8)
hold on;
grid on;
patch([Xs, fliplr(Xs)],[mdPost_disc-2*sdPost_disc; flipud(mdPost_disc+2*sdPost_disc)], 1, 'FaceColor', [0.9,0.9,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
patch([Xs, fliplr(Xs)],[mdPost_disc-sdPost_disc; flipud(mdPost_disc+sdPost_disc)], 1, 'FaceColor', [0.8,0.8,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(Xs, mdPost_disc, 'b-', 'LineWidth', 1); % We plot the mean line.
title('radome slope VS look angle (GPR discrete)')
xlabel('\theta_s (deg)')
ylabel('\rho_theta (deg/deg)')
axis square

% (2.4.3) estimate

Xm_est=rad2deg(theta_s_est(index_m_est));
Xm_est=Xm_est';
fmh_est=rad2deg(lambda_r_est(index_m_est));

% tuning hyperparameters
lf_est=1;
lx_est=1;
% sfm_est=0.1;
sfm_est=0.5;
% hyp=function_gp_hyperparam(Xm_est,fmh_est);
% lx_est = sqrt(hyp(1));
% lf_est = sqrt(hyp(2));
% sfm_est = sqrt(hyp(3));

Xs_est=rad2deg(theta_s_est(index_s_est));
Xs_est=Xs_est';

% GPR
[mPost_est,SPost_est,mdPost_est,SdPost_est]=function_gp_dgp(Xs_est,Xm_est,fmh_est,lf_est,lx_est,sfm_est);
sPost_est=sqrt(diag(SPost_est)); % These are the posterior standard deviations.
sdPost_est=sqrt(diag(SdPost_est)); % These are the posterior standard deviations of the derivative.
% GP plot
% [h_5,h_6]=function_gp_plot(Xs,Xm,fmh,mPost_est,sPost_est,mdPost_est,sdPost_est);
subplot(3,3,6)
hold on;
grid on;
patch([Xs, fliplr(Xs)],[mPost_est-2*sPost_est; flipud(mPost_est+2*sPost_est)], 1, 'FaceColor', [0.9,0.9,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
patch([Xs, fliplr(Xs)],[mPost_est-sPost_est; flipud(mPost_est+sPost_est)], 1, 'FaceColor', [0.8,0.8,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(Xs, mPost_est, 'b-', 'LineWidth', 1); % We plot the mean line.
plot(Xm, fmh, 'ro'); % We plot the measurement points.
title('radome error VS look angle (GPR estimate)')
xlabel('\theta_s (deg)')
ylabel('\lambda_r (deg)')
axis square
subplot(3,3,9)
hold on;
grid on;
patch([Xs, fliplr(Xs)],[mdPost_est-2*sdPost_est; flipud(mdPost_est+2*sdPost_est)], 1, 'FaceColor', [0.9,0.9,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
patch([Xs, fliplr(Xs)],[mdPost_est-sdPost_est; flipud(mdPost_est+sdPost_est)], 1, 'FaceColor', [0.8,0.8,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(Xs, mdPost_est, 'b-', 'LineWidth', 1); % We plot the mean line.
title('radome slope VS look angle (GPR estimate)')
xlabel('\theta_s (deg)')
ylabel('\rho_theta (deg/deg)')
axis square