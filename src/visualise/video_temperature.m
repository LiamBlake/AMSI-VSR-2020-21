% Required directories
addpath visualise


h = figure;
axis tight manual
filename = '../figures/atlantic_temp.gif';
for t = 1:size(temp,3)
    % Draw plot for y = x.^n
    cbar = colourplot(lon(xilim(1):xilim(2)),lat(yilim(1):yilim(2)),temp(:,:,t),"Atlantic Temperature", "Longitude", "Latitude", "(Normalised) Temperature",jet);
    cbar.Limits = [0 1];
    drawnow
    
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if t == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end
