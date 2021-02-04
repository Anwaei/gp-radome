function [h_1,h_2]=function_gp_plot(Xs,Xm,fmh,mPost,sPost,mdPost,sdPost)

h_1=figure;
clf(h_1);
hold on;
grid on;
patch([Xs, fliplr(Xs)],[mPost-2*sPost; flipud(mPost+2*sPost)], 1, 'FaceColor', [0.9,0.9,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
patch([Xs, fliplr(Xs)],[mPost-sPost; flipud(mPost+sPost)], 1, 'FaceColor', [0.8,0.8,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(Xs, mPost, 'b-', 'LineWidth', 1); % We plot the mean line.
plot(Xm, fmh, 'ro'); % We plot the measurement points.
h_2=figure;
clf(h_2);
hold on;
grid on;
patch([Xs, fliplr(Xs)],[mdPost-2*sdPost; flipud(mdPost+2*sdPost)], 1, 'FaceColor', [0.9,0.9,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
patch([Xs, fliplr(Xs)],[mdPost-sdPost; flipud(mdPost+sdPost)], 1, 'FaceColor', [0.8,0.8,1], 'EdgeColor', 'none'); % This is the grey area in the plot.
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(Xs, mdPost, 'b-', 'LineWidth', 1); % We plot the mean line.