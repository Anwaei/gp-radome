function [x_imm_corr,P_imm_corr,mu_curr]=function_imm(x_imm_prev,P_imm_prev,mu_prev,z,Q,R,delta_t_disc,delta_t_cont,rho_theta,T_alpha,V_M,u_cont,i_disc)


nModes=size(rho_theta,2);

x_imm_corr=cell(nModes,1);
P_imm_corr=cell(nModes,1);

mu_curr=zeros(nModes,1);

mu0=ones(nModes,1)/nModes;

p=zeros(nModes,nModes);
for i=1:nModes
%     p(i,i)=0.9;
    p(i,i)=0.99;
end
for i=1:nModes
    for j=1:nModes
        if i~=j
            p(i,j)=(1-p(i,i))/(nModes-1);
        end
    end
end

%% step 1
if i_disc==1
    mu_km1=mu0;
else
    mu_km1=mu_prev;
end

%% calculate the normalization constants c_bar
c_bar=zeros(nModes,1);
for j=1:nModes
    for k=1:nModes
        if k==1
            c_bar(j)=p(k,j)*mu_km1(k);
        else
            c_bar(j)=c_bar(j)+p(k,j)*mu_km1(k);
        end
    end
end

%% calculate the mixing probabilities
mu_cond_km1=zeros(nModes,nModes);
for j=1:nModes
    for k=1:nModes
        mu_cond_km1(j,k)=p(j,k)*mu_km1(j)/c_bar(k);
    end
end

%% step 2
%% get previous estimate matched to mode j
x_km1=cell(nModes,1);
P_km1=cell(nModes,1);
for j=1:nModes
        x_km1{j}=x_imm_prev{j};     
        P_km1{j}=P_imm_prev{j};
end

%% initial mixed initial condition
x0_km1=cell(nModes,1);
% initial
for j=1:nModes
   x0_km1{j}=zeros(3,1);
end
% update
for j=1:nModes
    for k=1:nModes
        if k==1
            x0_km1{j}=mu_cond_km1(k,j)* x_km1{k};
        else
            x0_km1{j}=x0_km1{j}+mu_cond_km1(k,j)*x_km1{k};
        end
    end
end
P0_km1=cell(nModes,1);
% initial
for j=1:nModes
    P0_km1{j}=zeros(3,3);
end
% update
for j=1:nModes
    for k=1:nModes
        if k==1
            P0_km1{j}=mu_cond_km1(k,j)*(P_km1{k}+(x_km1{k}-x0_km1{j})*(x_km1{k}-x0_km1{j})');
        else
            P0_km1{j}=P0_km1{j}+mu_cond_km1(k,j)*(P_km1{k}+(x_km1{k}-x0_km1{j})*(x_km1{k}-x0_km1{j})');
        end
    end
end

%% step 3
% mode-matched filtering
Lambda=zeros(nModes,1);
for j=1:nModes
    [x_pred,P_pred,x_corr,P_corr,likelihood]=function_ukf(x0_km1{j},P0_km1{j},z,Q,R,delta_t_disc,delta_t_cont,rho_theta(:,j),T_alpha,V_M,u_cont,i_disc); 
    % calculate the likelihood functions matched to mode j
    Lambda(j)=likelihood;
    % log
    x_imm_corr{j}=x_corr;
    P_imm_corr{j}=P_corr;
end

%% step 4
% initial
c=0;
% calculate the normalization constants c
for j=1:nModes
    if j==1
        c=Lambda(j)*c_bar(j);
    else
        c=c+Lambda(j)*c_bar(j);
    end
end

%% mode probability update
for j=1:nModes
    mu_curr(j)=Lambda(j)*c_bar(j)/c;
end