% function [y,Y,P,Y1]=ut(f,X,Wm,Wc,n,R)
function [y,Y,P,Y1]=ut_h(X,Wm,Wc,n,R,delta_t_disc,delta_t_cont,rho_theta_est,T_alpha,V_M,u_cont,i_disc)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed smapling points
%        P: transformed covariance
%       Y1: transformed deviations

L=size(X,2);
y=zeros(n,1);
Y=zeros(n,L);
for k=1:L                   
%     Y(:,k)=f(X(:,k));   
    state=X(:,k);
    Y(:,k)=function_radome_meas(state,delta_t_disc,delta_t_cont,rho_theta_est,T_alpha,V_M,u_cont,i_disc);
    y=y+Wm(k)*Y(:,k);       
end
Y1=Y-y(:,ones(1,L));
P=Y1*diag(Wc)*Y1'+R;   