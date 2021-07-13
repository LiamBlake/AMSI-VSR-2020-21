function [u,v,x0,y0,t] = load_checkerboard()
%LOAD_CHECKERBOARD Summary of this function goes here
%   Detailed explanation goes here

% Code adapted from 

% Directory containing .mat files
files = dir('D:\liaml\Documents\Projects\LCS-Clustering\data\checkerboard\0*.mat');
region=[0 0 2300 1700];

[x0,y0]=meshgrid(region(1):5:region(1)+region(3), ...
    region(2):5:region(2)+region(4));

% Constants
FR = 60;     % Framerate

% Load each file
for ii=1:numel(files)
  disp(['Processing ' files(ii).name]);
  [pvx,pvy,px,py]=vp('D:\liaml\Documents\Projects\LCS-Clustering\data\checkerboard\',3+ii,region,0);
  u(:,:,ii) = griddata(px(:),py(:),pvx,x0,y0,'cubic');
  v(:,:,ii) = griddata(px(:),py(:),pvy,x0,y0,'cubic');
  
end

t = 1/FR*1:size(u,3);


end

