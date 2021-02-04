%% We set up the workspace, ready for executing scripts.
close all
clear all; % Empty the workspace.
clc; % Empty the command window.
exportFigs = 0; % Do we export figures? 0 for no, 1 (or anything else) for yes.
useColor = 1; % Should we set up plots for colored output (1) or black-and-white output (0)?
addpath('C:\Users\admin\Documents\GitHub\radome\gp\ExportFig'); % We add the functions for exporting figures.

%% We define data.
lf = 1; % This is the output length scale.
lx = 1; % This is the input length scale.
sfm = 0.1; % This is the output noise scale.
Xs = 0:0.01:4; % These are the trial points.
Xm = [1,2.5,3.7,0.1]; % These are the measurement points.
fmh = cos(Xm - 2)' + sfm*randn(size(Xm))'; % These are the measurement values, corrupted by noise.

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
sPost = sqrt(diag(SPost)); % These are the posterior standard deviations.

%% We set up the GP plot.
figure(1);
clf(1);
hold on;
grid on;
xlabel('Input');
ylabel('Output');
if useColor == 0
	patch([Xs, fliplr(Xs)],[mPost-2*sPost; flipud(mPost+2*sPost)], 1, 'FaceColor', [1,1,1]*0.9, 'EdgeColor', 'none'); % This is the grey area in the plot.
	patch([Xs, fliplr(Xs)],[mPost-sPost; flipud(mPost+sPost)], 1, 'FaceColor', [1,1,1]*0.8, 'EdgeColor', 'none'); % This is the grey area in the plot.
	set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
	plot(Xs, mPost, 'k-', 'LineWidth', 1); % We plot the mean line.
	plot(Xm, fmh, 'ko'); % We plot the measurement points.
else
	patch([Xs, fliplr(Xs)],[mPost-2*sPost; flipud(mPost+2*sPost)], 1, 'FaceColor', [0.9,0.9,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
	patch([Xs, fliplr(Xs)],[mPost-sPost; flipud(mPost+sPost)], 1, 'FaceColor', [0.8,0.8,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
	set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
	plot(Xs, mPost, 'b-', 'LineWidth', 1); % We plot the mean line.
	plot(Xm, fmh, 'ro'); % We plot the measurement points.
end
axis([0,4,-1.5,2.5]);
if exportFigs ~= 0
	export_fig('OriginalGP.png','-transparent');
end

%% We now take derivatives of the posterior distribution.
d2Kss = Kss.*(1/lx^2 - (diff(nm+1:end,nm+1:end).^2/lx^4)); % This is the derivative d^2 k(x,x') / dx dx' for the trial points.
dKms = Kms.*diff(1:nm,nm+1:end)/lx^2; % This is the derivative dk(x,x') / dx' where the first input is the set of measurement points and the second input is the set of trial points.
dKsm = dKms'; % This is the derivative dk(x,x') / dx where the first input is the set of trial points and the second input is the set of measurement points.
mdPost = -dKsm/(Kmm + Sfm)*(fmh - mm); % These are the posterior mean of the derivative.
SdPost = d2Kss - dKsm/(Kmm + Sfm)*dKms; % These are the posterior covariance of the derivative.
sdPost = sqrt(diag(SdPost)); % These are the posterior standard deviations of the derivative.

%% We set up the GP plot.
figure(2);
clf(2);
hold on;
grid on;
xlabel('Input');
ylabel('Output');
if useColor == 0
	patch([Xs, fliplr(Xs)],[mdPost-2*sdPost; flipud(mdPost+2*sdPost)], 1, 'FaceColor', [1,1,1]*0.9, 'EdgeColor', 'none'); % This is the grey area in the plot.
	patch([Xs, fliplr(Xs)],[mdPost-sdPost; flipud(mdPost+sdPost)], 1, 'FaceColor', [1,1,1]*0.8, 'EdgeColor', 'none'); % This is the grey area in the plot.
	set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
	plot(Xs, mdPost, 'k-', 'LineWidth', 1); % We plot the mean line.
else
	patch([Xs, fliplr(Xs)],[mdPost-2*sdPost; flipud(mdPost+2*sdPost)], 1, 'FaceColor', [0.9,0.9,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
	patch([Xs, fliplr(Xs)],[mdPost-sdPost; flipud(mdPost+sdPost)], 1, 'FaceColor', [0.8,0.8,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
	set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
	plot(Xs, mdPost, 'b-', 'LineWidth', 1); % We plot the mean line.
end
axis([0,4,-2,2]);
if exportFigs ~= 0
	export_fig('DerivativeGP.png','-transparent');
end