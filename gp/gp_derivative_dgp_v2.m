% rng('default')
%% We set up the workspace, ready for executing scripts.
close all
clear; % Empty the workspace.
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

%% GPR
[mPost,SPost,mdPost,SdPost] = function_gp_dgp(Xs,Xm,fmh,lf,lx,sfm);
sPost = sqrt(diag(SPost)); % These are the posterior standard deviations.
sdPost = sqrt(diag(SdPost)); % These are the posterior standard deviations of the derivative.

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

%% We set up the dGP plot.
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