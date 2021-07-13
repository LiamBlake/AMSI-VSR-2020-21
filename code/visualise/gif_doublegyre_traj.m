addpath visualise

conFigure(29, 4/3);

% Load data
load('../data/doublegyre.mat')

X = doublegyre.traj.X;
Y = doublegyre.traj.Y;
T = linspace(doublegyre.res.tlim(1),doublegyre.res.tlim(2),doublegyre.res.Nt);

h = figure('color', [1 1 1]);
axis tight manual
filename = '../figures/doublegyre.gif';
for t = 1:length(T)
    % Draw plot for y = x.^n
    plot(X(:,:,t), Y(:,:,t), 'b.');
    xlim([0 2]); ylim([0 1]);
    xlabel("$x$", "Interpreter", "LaTeX"); ylabel("$y$", "Interpreter", "LaTeX");
    drawnow
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if t == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', .01);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'delaytime', .01);
    end
end
