function plot_trajectory_centres(centres)
%PLOT_TRAJECTORY_CENTRES Summary of this function goes here
%   Detailed explanation goes here

centres = squeeze(centres);

figure; hold on;
for i = 1:size(centres,1)
    plot(centres(i,1:2:end),centres(i,2:2:end))
end
