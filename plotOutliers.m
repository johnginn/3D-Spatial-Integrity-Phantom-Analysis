% This function plots the outliers that the user is going to remove in red
% as well as shows the volumetric dataset so that the locations of theses
% spheres can be confirmed in the phantom volume.
%
% Input:
% xSphere (index) The found x-location of the spheres after corrections are applied
% ySphere (index) The found y-location of the spheres after corrections are applied
% zSphere (index) The found z-location of the spheres after corrections are applied
% diffX (index) The difference btwn the x-loctation of the found and ground truth pos
% diffY (index) The difference btwn the y-loctation of the found and ground truth pos
% diffZ (index) The difference btwn the z-loctation of the found and ground truth pos
% xGndTruthFinal (index) The final ground-truth x-position
% yGndTruthFinal (index) The final ground-truth y-position
% zGndTruthFinal (index) The final ground-truth z-position 
% radiusAll (mm) The distance of the spheres from isocenter
% finalDiffCorr (mm) Deviation of found location to ground truth location in 3D
% radiusThreshold1 (mm) Distance specifying radius 1 region for plotting 
% volData The volumetric dataset
% sphRangeX (index) The radius of the sphere model in the x-dimension 
% sphRangeY (index) The radius of the sphere model in the y-dimension 
% sphRangeZ (index) The radius of the sphere model in the z-dimension 
% sphereVol The sphere template volume
% radiusRemove [min max] (mm) The radii specifying the region of spheres to be removed
% rangeDevRemove [min max] (mm) The deviations specifying the region of spheres to be removed
% plotGUI (y=1,n=0) Whether or not you want the gui to be plotted as well.
%
% Output:
% numRemoved The number of spheres being removed
%
% John Ginn
% Created: 8/1/16
% Modified: 9/22/16

function [numRemoved] = plotOutliers(xSphere,ySphere,zSphere,diffX,diffY,diffZ,...
    xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,radiusAll,finalDiffCorr,...
    radiusThreshold1,volData,sphRangeX,sphRangeY,sphRangeZ,sphereVol,radiusRemove,rangeDevRemove,plotGUI)


countNorm = 0;
countHighlight = 0;
sphereImg = zeros(length(volData(:,1,1)),length(volData(1,:,1)),length(volData(1,1,:)));
blankImg = sphereImg;
sphereFitPtsPass = sphereImg; % store just the location of the spheres
sphereFitPtsFail = sphereImg; % store just the location of the spheres
grndTruthPts = sphereImg; % the location of the "ground truth" points
for step = 1:length(radiusAll)
        % round locations and flip x,y for plotting
    data = round([ySphere(step) xSphere(step) zSphere(step)]);
    xLoc = (data(2)-sphRangeX):1:(data(2)+sphRangeX);
    yLoc = (data(1)-sphRangeY):1:(data(1)+sphRangeY); 
    zLoc = (data(3)-sphRangeZ):1:(data(3)+sphRangeZ); 
    xMarkerLoc = (data(1)-1):1:(data(1)+1); % make + symbol on image
    yMarkerLoc = (data(2)-1):1:(data(2)+1);% make + symbol on image
    if (radiusAll(step)>=radiusRemove(1))&&(radiusAll(step)<=radiusRemove(2))...
            &&(finalDiffCorr(step)>=rangeDevRemove(1))&&(finalDiffCorr(step)<=rangeDevRemove(2)) 
        % markers based on whether or not in selected range
        sphereFitPtsFail(xMarkerLoc,yMarkerLoc,data(3)) = makeShape(length(xMarkerLoc),'x','thick'); % make a x        
        % plot as red to show these locations on all the plots
        countHighlight = countHighlight + 1;
        xDistHighlight(countHighlight) = ((xSphere(step))); % index in the image
        xDiffHighlight(countHighlight) = diffX(step);
        yDistHighlight(countHighlight) = ((ySphere(step)));
        yDiffHighlight(countHighlight) = diffY(step);
        zDistHighlight(countHighlight) = ((zSphere(step)));
        zDiffHighlight(countHighlight) = diffZ(step);    
        radiusHighlight(countHighlight) = radiusAll(step);
        totDiffHighlight(countHighlight) = finalDiffCorr(step);
    else
        % markers based on whether or not in selected range
        sphereFitPtsPass(xMarkerLoc,yMarkerLoc,data(3)) = makeShape(length(xMarkerLoc),'square','thin'); % make a square
        % plot as blue to show these locations on all the plots
        countNorm = countNorm + 1;
        xDistNorm(countNorm) = ((xSphere(step)));
        xDiffNorm(countNorm) = diffX(step);
        yDistNorm(countNorm) = ((ySphere(step)));
        yDiffNorm(countNorm) = diffY(step);
        zDistNorm(countNorm) = ((zSphere(step)));
        zDiffNorm(countNorm) = diffZ(step);
        radiusNorm(countNorm) = radiusAll(step);
        totDiffNorm(countNorm) = finalDiffCorr(step);
    end
    % store location of ground truth for spheres
    data = round([yGndTruthFinal(step) xGndTruthFinal(step) zGndTruthFinal(step)]);
    xMarkerLoc = (data(1)-1):1:(data(1)+1); % make + symbol on image
    
    yMarkerLoc = (data(2)-1):1:(data(2)+1);% make + symbol on image
    grndTruthPts(xMarkerLoc,data(2),data(3)) = 1; % - portion of + marker for ground truth of spheres
    grndTruthPts(data(1),yMarkerLoc,data(3)) = 1; % - portion of + marker for ground truth of spheres

end

figure;
maxDiff = max([max(diffX),max(diffY),max(diffZ),finalDiffCorr]);

subplot(2,2,1)
plot(xDistNorm,xDiffNorm,'.','MarkerSize',8)
hold on
plot(xDistHighlight,xDiffHighlight,'.r','MarkerSize',8)
xlabel('Image position (index)','FontSize',22)
ylabel('Deviation in X (index)','FontSize',22)
title('Deviation in X-Direction','FontSize',22)
% ylim([0 maxDiff*1.1])
% xlim([0 maxPosX*1.1])
subplot(2,2,2)
plot(yDistNorm,yDiffNorm,'.','MarkerSize',8)
hold on
plot(yDistHighlight,yDiffHighlight,'.r','MarkerSize',8)
xlabel('Image position (index)','FontSize',22)
ylabel('Deviation in Y (index)','FontSize',22)
title('Deviation in Y-Direction','FontSize',22)
% ylim([0 maxDiff*1.1])
% xlim([0 maxPosY*1.1])
subplot(2,2,3)
plot(zDistNorm,zDiffNorm,'.','MarkerSize',8)
hold on
plot(zDistHighlight,zDiffHighlight,'.r','MarkerSize',8)
xlabel('Image position (index)','FontSize',22)
ylabel('Deviation in Z (index)','FontSize',22)
title('Deviation in Z-Direction','FontSize',22)
% ylim([0 maxDiff*1.1])
% xlim([0 maxPosZ*1.1])
subplot(2,2,4)
plot(radiusNorm,totDiffNorm,'.','MarkerSize',8)
hold on
plot(radiusHighlight,totDiffHighlight,'.r','MarkerSize',8)
% ylim([0 maxDiff*1.1])
% xlim([0 maxPos*1.1])
xlabel('Distance from Isocenter (index)','FontSize',22)
ylabel('Deviation (index)','FontSize',22)
% title('3D Deviation','FontSize',22)
numRemoved = length(totDiffHighlight);


% markers based on deviation tolerance
% % Plot the final sphere locations with the image
% plotTolerance = radiusThreshold1; % the tolerance that defines when spheres
% % will be plotted as red or green on the CheckCorrGUI
% count = 0;
% sphereImg = zeros(length(volData(:,1,1)),length(volData(1,:,1)),length(volData(1,1,:)));
% blankImg = sphereImg;
% sphereFitPtsPass = sphereImg; % store just the location of the spheres
% sphereFitPtsFail = sphereImg; % store just the location of the spheres
% grndTruthPts = sphereImg; % the location of the "ground truth" points
% for stepGui = 1:length(xSphere)
%     % round locations and flip x,y for plotting
%     data = round([ySphere(stepGui) xSphere(stepGui) zSphere(stepGui)]);
%     xLoc = (data(2)-sphRangeX):1:(data(2)+sphRangeX);
%     yLoc = (data(1)-sphRangeY):1:(data(1)+sphRangeY); 
%     zLoc = (data(3)-sphRangeZ):1:(data(3)+sphRangeZ); 
%     xMarkerLoc = (data(1)-1):1:(data(1)+1); % make + symbol on image
%     yMarkerLoc = (data(2)-1):1:(data(2)+1);% make + symbol on image
%     sphereImg(yLoc,xLoc,zLoc) = sphereVol; % store the sphere volume
%     if finalDiffCorr(stepGui) <= plotTolerance
%         sphereFitPtsPass(xMarkerLoc,yMarkerLoc,data(3)) = makeShape(length(xMarkerLoc),'square','thin'); % make a square
%     else
%         sphereFitPtsFail(xMarkerLoc,yMarkerLoc,data(3)) = makeShape(length(xMarkerLoc),'x','thick'); % make a x
%     end
%     % store location of ground truth for spheres
%     data = round([yGndTruthFinal(stepGui) xGndTruthFinal(stepGui) zGndTruthFinal(stepGui)]);
%     xMarkerLoc = (data(1)-1):1:(data(1)+1); % make + symbol on image
%     
%     yMarkerLoc = (data(2)-1):1:(data(2)+1);% make + symbol on image
%     grndTruthPts(xMarkerLoc,data(2),data(3)) = 1; % - portion of + marker for ground truth of spheres
%     grndTruthPts(data(1),yMarkerLoc,data(3)) = 1; % - portion of + marker for ground truth of spheres
%     % ground truth points that did not meet correlation requirement
% end
% Use a GUI to verify the spheres were appropriately found
if plotGUI == 1
    guiFindSphereData{1} = volData; % phantom image data
    guiFindSphereData{2} = blankImg; % sphere image data
    guiFindSphereData{3} = sphereFitPtsPass; % marker positions for fusing images
    guiFindSphereData{4} = sphereFitPtsFail; % marker positions for fusing images
    guiFindSphereData{5} = grndTruthPts; % marker positions for the ground truth
    CheckCorrGUI2(guiFindSphereData)
end

end