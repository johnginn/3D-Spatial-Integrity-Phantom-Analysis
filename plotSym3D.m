% This function plots the volumetric data in 3D with symmetric axes so that 
% the data does not look falsely warped in the plotting 
%
% Input:
% xSphere X locations of the spheres
% ySphere Y locations of the spheres
% zSphere Z locations of the spheres
% xGndTruth X ground-truth locations of the spheres
% yGndTruth Y ground-truth locations of the spheres
% zGndTruth Z ground-truth locations of the spheres
%
% Output:
%
% John Ginn
% Created: 8/11/16
% Modified: 8/11/16

function [] = plotSym3D(xSphere,ySphere,zSphere,xGndTruth,yGndTruth,zGndTruth)

figure;
scatter3(xSphere,ySphere,zSphere)
hold on
scatter3(xGndTruth,yGndTruth,zGndTruth,'r+')
title('Sphere and Ground Truth Comparison','FontSize',20)
legend('Sphere Locations','Ground Truth')
xlabel('x-axis','FontSize',20)
ylabel('y-axis','FontSize',20)
zlabel('z-axis','FontSize',20)
% make sure axes have the same dimensions so the data doesn't look
% falsely warped
xSimPlotRange = [min([min(xSphere),min(xGndTruth)])...
    max([max(xSphere),max(xGndTruth)])];
ySimPlotRange = [min([min(ySphere),min(yGndTruth)])...
    max([max(ySphere),max(yGndTruth)])];
zSimPlotRange = [min([min(zSphere),min(zGndTruth)])...
    max([max(zSphere),max(zGndTruth)])];
maxSimPlotRange = max([(xSimPlotRange(2) - xSimPlotRange(1)),...
    (ySimPlotRange(2) - ySimPlotRange(1)),...
    (zSimPlotRange(2) - zSimPlotRange(1))]);
plotDist = ceil(maxSimPlotRange/2 + 0.1*maxSimPlotRange);
xSimPlotRange = [round(mean(xSimPlotRange)) - plotDist,...
    round(mean(xSimPlotRange)) + plotDist];
ySimPlotRange = [round(mean(ySimPlotRange)) - plotDist,...
    round(mean(ySimPlotRange)) + plotDist];
zSimPlotRange = [round(mean(zSimPlotRange)) - plotDist,...
    round(mean(zSimPlotRange)) + plotDist];
xlim(xSimPlotRange)
ylim(ySimPlotRange)
zlim(zSimPlotRange)
    
end