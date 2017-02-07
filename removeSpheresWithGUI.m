% This function removes spheres with the RemoveSpheresGUI.
%
% Input:
% finalSphereDataXCorr The sphere x-locations after corrections
% finalSphereDataYCorr The sphere y-locations after corrections
% finalSphereDataYCorr The sphere z-locations after corrections
% diffX (index) The difference btwn the x-loctation of the found and ground truth pos
% diffY (index) The difference btwn the y-loctation of the found and ground truth pos
% diffZ (index) The difference btwn the z-loctation of the found and ground truth pos
% xGndTruthFinal (index) The final ground-truth x-position
% yGndTruthFinal (index) The final ground-truth y-position
% zGndTruthFinal (index) The final ground-truth z-position 
% radiusAll (mm) The distance of the spheres from isocenter
% finalDiffCorr (mm) Deviation of found location to ground truth location in 3D
% volData The volumetric dataset
% xSphere The sphere x-locations before corrections
% ySphere The sphere y-locations before corrections
% zSphere The sphere z-locations before corrections
%
% Output:
% numRemoved The number of spheres being removed
% xDataCorrRemove Sphere x-locations after correction, after certain spheres have been removed
% yDataCorrRemove Sphere y-locations after correction, after certain spheres have been removed
% zDataCorrRemove Sphere z-locations after correction, after certain spheres have been removed
% xGndTruthRemove Ground truth x-locations after certain spheres have been removed
% yGndTruthRemove Ground truth y-locations after certain spheres have been removed
% zGndTruthRemove Ground truth z-locations after certain spheres have been removed
% radiusRemove Radius of spheres to isocenter after certain spheres have been removed
% diffRemove Deviation of spheres from ground truth after certain spheres have been removed
% diffXRemove Deviation of spheres from ground truth along x after certain spheres have been removed
% diffYRemove Deviation of spheres from ground truth along x after certain spheres have been removed
% diffZRemove Deviation of spheres from ground truth along x after certain spheres have been removed
% xSphereRemove Sphere x-locations before correction, after certain spheres have been removed
% ySphereRemove Sphere y-locations before correction, after certain spheres have been removed
% zSphereRemove Sphere z-locations before correction, after certain spheres have been removed
%
% John Ginn
% Created: 8/1/16
% Modified: 9/22/16

function [numRemoved,xDataCorrRemove,yDataCorrRemove,zDataCorrRemove,...
    xGndTruthRemove,yGndTruthRemove,zGndTruthRemove,radiusRemove,diffRemove,...
    diffXRemove,diffYRemove,diffZRemove,xSphereRemove,ySphereRemove,zSphereRemove] = ...
    removeSpheresWithGUI(finalSphereDataXCorr,finalSphereDataYCorr,finalSphereDataZCorr,diffX,diffY,diffZ,...
    xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,radiusAll,diffCorr,...
    volData,xSphere,ySphere,zSphere)


sphereImg = zeros(length(volData(:,1,1)),length(volData(1,:,1)),length(volData(1,1,:)));
blankImg = sphereImg;
sphereFitPtsPass = sphereImg; % store just the location of the spheres
sphereFitPtsFail = sphereImg; % store just the location of the spheres
grndTruthPts = sphereImg; % the location of the "ground truth" points
for step = 1:length(radiusAll)
    % round locations and flip x,y for plotting
    data = round([finalSphereDataYCorr(step) finalSphereDataXCorr(step) finalSphereDataZCorr(step)]);
    xMarkerLoc = (data(1)-1):1:(data(1)+1); % make + symbol on image
    yMarkerLoc = (data(2)-1):1:(data(2)+1);% make + symbol on image
    % markers based on whether or not in selected range
    sphereFitPtsPass(xMarkerLoc,yMarkerLoc,data(3)) = makeShape(length(xMarkerLoc),'square','thin'); % make a square
    % store location of ground truth for spheres
    data = round([yGndTruthFinal(step) xGndTruthFinal(step) zGndTruthFinal(step)]);
    xMarkerLoc = (data(1)-1):1:(data(1)+1); % make + symbol on image
    yMarkerLoc = (data(2)-1):1:(data(2)+1);% make + symbol on image
    grndTruthPts(xMarkerLoc,data(2),data(3)) = 1; % - portion of + marker for ground truth of spheres
    grndTruthPts(data(1),yMarkerLoc,data(3)) = 1; % - portion of + marker for ground truth of spheres
    
end

figure;
maxDiff = max([max(diffX),max(diffY),max(diffZ),diffCorr]);

subplot(2,2,1)
plot(finalSphereDataXCorr,diffX,'.','MarkerSize',8)
xlabel('Image position (index)','FontSize',22)
ylabel('Deviation in X (index)','FontSize',22)
title('Deviation in X-Direction','FontSize',22)
% ylim([0 maxDiff*1.1])
% xlim([0 maxPosX*1.1])
subplot(2,2,2)
plot(finalSphereDataYCorr,diffY,'.','MarkerSize',8)
xlabel('Image position (index)','FontSize',22)
ylabel('Deviation in Y (index)','FontSize',22)
title('Deviation in Y-Direction','FontSize',22)
% ylim([0 maxDiff*1.1])
% xlim([0 maxPosY*1.1])
subplot(2,2,3)
plot(finalSphereDataZCorr,diffZ,'.','MarkerSize',8)
xlabel('Image position (index)','FontSize',22)
ylabel('Deviation in Z (index)','FontSize',22)
title('Deviation in Z-Direction','FontSize',22)
% ylim([0 maxDiff*1.1])
% xlim([0 maxPosZ*1.1])
subplot(2,2,4)
plot(radiusAll,diffCorr,'.','MarkerSize',8)
% ylim([0 maxDiff*1.1])
% xlim([0 maxPos*1.1])
xlabel('Distance from Isocenter (index)','FontSize',22)
ylabel('Deviation (index)','FontSize',22)
title('3D Deviation','FontSize',22)


% Use a GUI to remove spheres from the analysis
guiFindSphereData{1} = volData; % phantom image data
guiFindSphereData{2} = blankImg; % sphere image data
guiFindSphereData{3} = sphereFitPtsPass; % marker positions for fusing images
guiFindSphereData{4} = sphereFitPtsFail; % marker positions for fusing images
guiFindSphereData{5} = grndTruthPts; % marker positions for the ground truth
disp('Close the GUI when all spheres have been selected')
spheresToRemove = RemoveSpheresGUI(guiFindSphereData);

if isempty(spheresToRemove)
    numRemoved = 0;
else
    numRemoved = length(spheresToRemove(:,1));
end

% select out the spheres to remove
if numRemoved > 0
    numOrigData = length(finalSphereDataXCorr);
    spheresToSkipInd = zeros(numRemoved,1);
    xDataCorrRemove = zeros((numOrigData-numRemoved),1);
    yDataCorrRemove = xDataCorrRemove;
    zDataCorrRemove = xDataCorrRemove;
    xGndTruthRemove = xDataCorrRemove;
    yGndTruthRemove = xDataCorrRemove;
    zGndTruthRemove = xDataCorrRemove;
    xSphereRemove = xDataCorrRemove;
    ySphereRemove = xDataCorrRemove;
    zSphereRemove = xDataCorrRemove;
    diffRemove = xDataCorrRemove'; % transpose to make original dimension
    diffXRemove = xDataCorrRemove;
    diffYRemove = xDataCorrRemove;
    diffZRemove = xDataCorrRemove;
    radiusRemove = xDataCorrRemove;
    for stepSpheres = 1:numRemoved;
        % calculate the distance to all the spheres
        distToSpheres = abs((finalSphereDataXCorr - spheresToRemove(stepSpheres,1)).^2 + ...
            (finalSphereDataYCorr - spheresToRemove(stepSpheres,2)).^2 + ...
            (finalSphereDataZCorr - spheresToRemove(stepSpheres,3)).^2);
        % find the closest sphere
        [minDist, sphereInd] = min(distToSpheres);
        % store the index of this sphere
        spheresToSkipInd(stepSpheres) = sphereInd;
    end
    % remove the spheres
    count = 0;
    for step = 1:length(xSphere)
        removeSphere = 0;
        for stepRemoveList = 1:numRemoved
            if spheresToSkipInd(stepRemoveList) == step
                % remove the sphere
                removeSphere = 1;
                break
            end
        end
        if removeSphere == 0
            count = count + 1;
            xDataCorrRemove(count) = finalSphereDataXCorr(step);
            yDataCorrRemove(count) = finalSphereDataYCorr(step);
            zDataCorrRemove(count) = finalSphereDataZCorr(step);
            xGndTruthRemove(count) = xGndTruthFinal(step);
            yGndTruthRemove(count) = yGndTruthFinal(step);
            zGndTruthRemove(count) = zGndTruthFinal(step);
            xSphereRemove(count) = xSphere(step);
            ySphereRemove(count) = ySphere(step);
            zSphereRemove(count) = zSphere(step);
            diffXRemove(count) = diffX(step);
            diffYRemove(count) = diffY(step);
            diffZRemove(count) = diffZ(step);
            diffRemove(count) = diffCorr(step);
            radiusRemove(count) = radiusAll(step);
        end
    end
else
    % don't remove any spheres, just return the original data
    xDataCorrRemove = finalSphereDataXCorr;
    yDataCorrRemove = finalSphereDataYCorr;
    zDataCorrRemove = finalSphereDataZCorr;
    xGndTruthRemove = xGndTruthFinal;
    yGndTruthRemove = yGndTruthFinal;
    zGndTruthRemove = zGndTruthFinal;
    xSphereRemove = xSphere;
    ySphereRemove = ySphere;
    zSphereRemove = zSphere;
    diffRemove = diffCorr;
    diffXRemove = diffX;
    diffYRemove = diffY;
    diffZRemove = diffZ;
    radiusRemove = radiusAll;
end