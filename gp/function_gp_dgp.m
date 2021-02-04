function [mPost,SPost,mdPost,SdPost] = function_gp_dgp(Xs,Xm,fmh,lf,lx,sfm)

%% We now set up the (squared exponential) covariance matrix and related terms.
nm = size(Xm,2); % This is the number of measurement points.
ns = size(Xs,2); % This is the number of trial points.
X = [Xm,Xs]; % We merge the measurement and trial points.
n = size(X,2); % This is the number of points.
diff = repmat(X,n,1) - repmat(X',1,n); % This is matrix containing differences between input points.
K = lf^2*exp(-1/2*diff.^2/lx^2); % This is the covariance matrix. It contains the covariances of each combination of points.
Kmm = K(1:nm,1:nm);
Kms = K(1:nm,nm+1:end);
Ksm = Kms';
Kss = K(nm+1:end,nm+1:end);
Sfm = sfm^2*eye(nm); % This is the noise covariance matrix.
mm = zeros(nm,1); % This is the mean vector m(Xm). We assume a zero mean function.
ms = zeros(ns,1); % This is the mean vector m(Xs). We assume a zero mean function.

%% Next, we apply GP regression.
mPost = ms + Ksm/(Kmm + Sfm)*(fmh - mm); % This is the posterior mean vector.
SPost = Kss - Ksm/(Kmm + Sfm)*Kms; % This is the posterior covariance matrix.

%% We now take derivatives of the posterior distribution.
% mdPost
alpha=inv(Kmm+Sfm)*(fmh-mm);
Gamma=lx^2;
mdPost=zeros(length(Xs),1);
for i=1:length(mdPost)
    xs=Xs(i);
    X_tilde_s=zeros(nm,1);
    for j=1:length(X_tilde_s)
        X_tilde_s(j)=xs-Xm(j);
    end
    mdPost(i)=-inv(Gamma)*X_tilde_s'*(Ksm(i,:)'.*alpha);
end
%d2Kss
d2Kss=zeros(length(Xs),length(Xs));
for i=1:length(Xs)
    for j=1:length(Xs)
        xs_i=Xs(i);
        xs_j=Xs(j);
        d2Kss(i,j)=inv(Gamma)*(1-(xs_i-xs_j)*(xs_i-xs_j)*inv(Gamma))*Kss(i,j);
    end
end
% dKms
dKms=zeros(length(Xm),length(Xs));
for i=1:length(Xm)
    for j=1:length(Xs)
        xm_i=Xm(i);
        xs_j=Xs(j);
        dKms(i,j)=inv(Gamma)*(xm_i-xs_j)*Kms(i,j);
    end
end
% dKsm
dKsm=zeros(length(Xs),length(Xm));
for i=1:length(Xs)
    for j=1:length(Xm)
        xs_i=Xs(i);
        xm_j=Xm(j);
        dKsm(i,j)=-inv(Gamma)*(xs_i-xm_j)*Ksm(i,j);
    end
end

assert(norm(dKms-dKsm')<1e-2)

SdPost = d2Kss - dKsm/(Kmm + Sfm)*dKms; % These are the posterior covariance of the derivative.