% Function to compare this analysis to the analysis code provided by
% ViewRay using the center transverse slice data.
%
% Input:
% sliceSphereX The x-location of the spheres in the slice extracted
% sliceSphereY The y-location of the spheres in the slice extracted
% sliceSphereZ The z-location of the spheres in the slice extacted
% sliceGroundX The x-location of the ground truth in the slices
% sliceGroundY The y-location of the ground truth in the slices
% sliceGroundZ The z-location of the ground truth in the slices
% voxelHeight (mm) The height of the voxels
% voxelWidth (mm) The width of the voxels
% voxelLength (mm) The length of the voxels
% centSliceImg The image data for the center slice
% shiftCorrection [sagittal,coronal,transverse] The correction required by fitting each center plane
% radiusThreshold1 The deviation tolerance first threshold to calculate the percent passing 
% radiusThreshold2 The deviation tolerance second threshold to calculate the percent passing 
% radiusThreshold3 The deviation tolerance third threshold to calculate the percent passing 
% radiusThreshold4 The deviation tolerance fourth threshold to calculate the percent passing 
% centerCol The index of the center column of the phantom
% centerSlice The index of the center slice of the phantom
% centerRow The index of the center row of the phantom
% radiusSearch1 (mm) The radius from isocenter that defines the subvolume for the first threshold
% radiusSearch2 (mm) The radius from isocenter that defines the subvolume for the second threshold
% radiusSearch3 (mm) The radius from isocenter that defines the subvolume for the third threshold
% radiusSearch3 (mm) The radius from isocenter that defines the subvolume for the fourth threshold
%
% Output:
% compareDeviationXY The deviation between each sphere location
% and ground truth in the X,Y plane
% compareDeviationXYZ The deviation between each sphere location
% and ground truth taking into account all three dimensions
% avgDevCompViewRayXY The average deviation between the sphere location
% and ground truth in the X,Y plane
% avgDevCompViewRayXYZ The average deviation between the sphere location
% and ground truth taking into account all three dimensions
% percentPassXY1 The percentage points passing only considering deviation in two dimensions (1st threshold)
% percentPassXYZ1 The percentage points passing  considering deviation in 3D (1st threshold)
% percentPassXY2 The percentage points passing only considering deviation in two dimensions (2nd threshold)
% percentPassXYZ2 The percentage points passing  considering deviation in 3D (2nd threshold)
% percentPassXY3 The percentage points passing only considering deviation in two dimensions (3rd threshold)
% percentPassXYZ3 The percentage points passing  considering deviation in 3D (3rd threshold)
% avgDevRadius1XY Average 2D deviation for spheres within radius 1
% avgDevRadius1XYZ Average 3D deviation for spheres within radius 1
% avgDevRadius2XY Average 2D deviation for spheres within radius 2
% avgDevRadius2XYZ Average 3D deviation for spheres within radius 2
% avgDevRadius3XY Average 2D deviation for spheres within radius 3
% avgDevRadius3XYZ Average 3D deviation for spheres within radius 3
% maxDevRadius1XY Maximum 2D deviation for spheres within radius 1
% maxDevRadius1XYZ Maximum 3D deviation for spheres within radius 1
% maxDevRadius2XY Maximum 2D deviation for spheres within radius 2
% maxDevRadius2XYZ Maximum 3D deviation for spheres within radius 2
% maxDevRadius3XY Maximum 2D deviation for spheres within radius 3
% maxDevRadius3XYZ Maximum 3D deviation for spheres within radius 3
% avgDevRadius4XY Average 2D deviation for spheres within radius 4
% avgDevRadius4XYZ Average 3D deviation for spheres within radius 4
% percentPassXY4 The percentage points passing only considering deviation in two dimensions (4th threshold)
% percentPassXYZ4 The percentage points passing  considering deviation in 3D (4th threshold)
% maxDevRadius4XY Maximum 2D deviation for spheres within radius 4
% maxDevRadius4XYZ Maximum 3D deviation for spheres within radius 4
%
% John Ginn
% Created: 7/7/16
% Modified: 1/24/17
function [compareDeviationXY,compareDeviationXYZ,avgDevCompViewRayXY,avgDevCompViewRayXYZ,...
    percentPassXY1,percentPassXYZ1,percentPassXY2,percentPassXYZ2,percentPassXY3,percentPassXYZ3,avgDevRadius1XY,...
    avgDevRadius1XYZ,avgDevRadius2XY,avgDevRadius2XYZ,avgDevRadius3XY,avgDevRadius3XYZ,...
    maxDevRadius1XY,maxDevRadius1XYZ,maxDevRadius2XY,maxDevRadius2XYZ,maxDevRadius3XY,maxDevRadius3XYZ,...
    avgDevRadius4XY,avgDevRadius4XYZ,percentPassXY4,percentPassXYZ4,maxDevRadius4XY,maxDevRadius4XYZ] = ...
    compViewRay(sliceSphereX,sliceSphereY,sliceSphereZ,sliceGroundX,sliceGroundY,sliceGroundZ,...
    voxelHeight,voxelWidth,voxelLength,centSliceImg,shiftCorrection,radiusThreshold1,radiusThreshold2,...
    radiusThreshold3,radiusThreshold4,centerCol,centerSlice,centerRow,radiusSearch1,radiusSearch2,radiusSearch3,radiusSearch4)
%% scale factor for upscaling the image
plotImage = 0; % y = 1, n = 0 Whether or not to plot the image
scaleFactor = 3;
plotThresh = '2D'; % use 2D or 3D deviation on plot
% calcualte the deviation for the center transverse slice to compare to the
% ViewRay analysis software
avgDevCompViewRayXY = 0;
avgDevCompViewRayXYZ = 0;
compareDeviationXY = zeros(1,length(sliceSphereX));
compareDeviationXYZ = zeros(1,length(sliceSphereX));
radius = zeros(1,length(sliceSphereX));
sphereFitPtsPass = zeros(scaleFactor*length(centSliceImg(:,1)),scaleFactor*length(centSliceImg(1,:)));
sphereFitPtsFail = sphereFitPtsPass;
grndTruthPts = sphereFitPtsPass;
plotDev = zeros(1,length(sliceSphereX));
countPassXY1 = 0; % 1st threshold
countPassXYZ1 = 0; % 1st threshold
countPassXY2 = 0; % 2nd threshold
countPassXYZ2 = 0; % 2nd threshold
countPassXY3 = 0; % 3rd threshold
countPassXYZ3 = 0; % 3rd threshold
countPassXY4 = 0; % 4th threshold
countPassXYZ4 = 0; % 4th threshold
countSpheresInRadius1 = 0; % the number of spheres within the 1st specified radius 
countSpheresInRadius2 = 0; % the number of spheres within the 2nd specified radius 
countSpheresInRadius3 = 0; % the number of spheres within the 3rd specified radius 
countSpheresInRadius4 = 0; % the number of spheres within the 4th specified radius 
% relative to isocenter
for step = 1:length(sliceSphereX)
    % shift the slice location according to what the model requires
    sliceSphereX(step) = shiftCorrection(1) + sliceSphereX(step);
    sliceSphereY(step) = shiftCorrection(2) + sliceSphereY(step);
    sliceSphereZ(step)= shiftCorrection(3) + sliceSphereZ(step);
    
    % current distance from isocenter
    radius(step) = sqrt((voxelWidth*(sliceSphereX(step) - centerCol)).^2 + ...
        (voxelHeight*(sliceSphereY(step) - centerRow)).^2 +...
        (voxelLength*(sliceSphereZ(step) - centerSlice)).^2);
    
    % calculate x,y deviation
    compareDeviationXY(step) =  sqrt((voxelWidth*(sliceSphereX(step) - sliceGroundX(step))).^2 + ...
        (voxelHeight*(sliceSphereY(step) - sliceGroundY(step))).^2);
    avgDevCompViewRayXY = avgDevCompViewRayXY + compareDeviationXY(step);

    % calcualte x,y,z deviation
    compareDeviationXYZ(step) =  sqrt((voxelWidth*(sliceSphereX(step) - sliceGroundX(step))).^2 + ...
        (voxelHeight*(sliceSphereY(step) - sliceGroundY(step))).^2 +...
        (voxelLength*(sliceSphereZ(step) - sliceGroundZ(step))).^2);
    avgDevCompViewRayXYZ = avgDevCompViewRayXYZ + compareDeviationXYZ(step);
    % check to see if point is within specified distance of isocenter
    if radius(step) <= radiusSearch1
        countSpheresInRadius1 = countSpheresInRadius1 + 1;
        radiusSmall1(countSpheresInRadius1) = radius(step);
        devSmallXY1(countSpheresInRadius1) = compareDeviationXY(step);
        devSmallXYZ1(countSpheresInRadius1) = compareDeviationXYZ(step);
        % check XY deviation to see if it passes the threshold
        if compareDeviationXY(step) <= radiusThreshold1
            countPassXY1 = countPassXY1 + 1;
        end
        % check XYZ deviation to see if it passes the threshold
        if compareDeviationXYZ(step) <= radiusThreshold1
            countPassXYZ1 = countPassXYZ1 + 1;
        end
    end
    if radius(step) <= radiusSearch2
        countSpheresInRadius2 = countSpheresInRadius2 + 1;
        radiusSmall2(countSpheresInRadius2) = radius(step);
        devSmallXY2(countSpheresInRadius2) = compareDeviationXY(step);
        devSmallXYZ2(countSpheresInRadius2) = compareDeviationXYZ(step);
        % check XY deviation to see if it passes the threshold
        if compareDeviationXY(step) <= radiusThreshold2
            countPassXY2 = countPassXY2 + 1;
        end
        % check XYZ deviation to see if it passes the threshold
        if compareDeviationXYZ(step) <= radiusThreshold2
            countPassXYZ2 = countPassXYZ2 + 1;
        end
    end
    if radius(step) <= radiusSearch3
        countSpheresInRadius3 = countSpheresInRadius3 + 1;
        radiusSmall3(countSpheresInRadius3) = radius(step);
        devSmallXY3(countSpheresInRadius3) = compareDeviationXY(step);
        devSmallXYZ3(countSpheresInRadius3) = compareDeviationXYZ(step);
        % check XY deviation to see if it passes the threshold
        if compareDeviationXY(step) <= radiusThreshold3
            countPassXY3 = countPassXY3 + 1;
        end
        % check XYZ deviation to see if it passes the threshold
        if compareDeviationXYZ(step) <= radiusThreshold3
            countPassXYZ3 = countPassXYZ3 + 1;
        end
    end
    if radius(step) <= radiusSearch4
        countSpheresInRadius4 = countSpheresInRadius4 + 1;
        radiusSmall4(countSpheresInRadius4) = radius(step);
        devSmallXY4(countSpheresInRadius4) = compareDeviationXY(step);
        devSmallXYZ4(countSpheresInRadius4) = compareDeviationXYZ(step);
        % check XY deviation to see if it passes the threshold
        if compareDeviationXY(step) <= radiusThreshold4
            countPassXY4 = countPassXY4 + 1;
        end
        % check XYZ deviation to see if it passes the threshold
        if compareDeviationXYZ(step) <= radiusThreshold4
            countPassXYZ4 = countPassXYZ4 + 1;
        end
    end
    % round locations for plotting, reverse points because of imshow
    % indexing
    data = round([sliceSphereY(step) sliceSphereX(step) sliceSphereZ(step)].*scaleFactor);
    xMarkerLoc = (data(1)-scaleFactor - 1):1:(data(1)+scaleFactor + 1); % make + symbol on image
    yMarkerLoc = (data(2)-scaleFactor - 1):1:(data(2)+scaleFactor + 1);% make + symbol on image
    
    % choose dataset for determining plot pass/fail
    if (strcmp(plotThresh,'2D') == 1)
        plotDev(step) = compareDeviationXY(step);
    else
        plotDev(step) = compareDeviationXYZ(step);
    end
    % the marker for the spheres that passed
    if plotDev(step) <= radiusThreshold1
%         sphereFitPtsPass(xMarkerLoc,data(2)) = 1; % - portion of + marker for fit of spheres
%         sphereFitPtsPass(data(1),yMarkerLoc) = 1; % | portion of + marker for fit of spheres
        sphereFitPtsPass(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'square','thick'); % make a square
    else % the marker for the spheres that failed
%         sphereFitPtsFail(xMarkerLoc,data(2)) = 1; % - portion of + marker for fit of spheres
%         sphereFitPtsFail(data(1),yMarkerLoc) = 1; % | portion of + marker for fit of spheres
         sphereFitPtsFail(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'x','thick'); % make a x
    end
    % Marker for fusing on the image (add on data in each slice)
    % marker{sliceNumber} = [current data; new data]
    % store location of ground truth for spheres
    data = round([sliceGroundY(step) sliceGroundX(step) sliceGroundZ(step)].*scaleFactor);
    xMarkerLoc = (data(1)-scaleFactor):1:(data(1)+scaleFactor); % make + symbol on image
    yMarkerLoc = (data(2)-scaleFactor):1:(data(2)+scaleFactor);% make + symbol on image
    grndTruthPts(xMarkerLoc,data(2)) = 1; % - portion of + marker for ground truth of spheres
    grndTruthPts(data(1),yMarkerLoc) = 1; % - portion of + marker for ground truth of spheres
    grndTruthPts(xMarkerLoc,data(2)-1) = 1; % - portion of + marker for ground truth of spheres
    grndTruthPts(data(1)-1,yMarkerLoc) = 1; % - portion of + marker for ground truth of spheres
    grndTruthPts(xMarkerLoc,data(2)+1) = 1; % - portion of + marker for ground truth of spheres
    grndTruthPts(data(1)+1,yMarkerLoc) = 1; % - portion of + marker for ground truth of spheres

end
avgDevRadius1XY = sum(devSmallXY1)/length(devSmallXY1);
avgDevRadius1XYZ = sum(devSmallXYZ1)/length(devSmallXYZ1);
avgDevRadius2XY = sum(devSmallXY2)/length(devSmallXY2);
avgDevRadius2XYZ = sum(devSmallXYZ2)/length(devSmallXYZ2);
avgDevRadius3XY = sum(devSmallXY3)/length(devSmallXY3);
avgDevRadius3XYZ = sum(devSmallXYZ3)/length(devSmallXYZ3);
avgDevRadius4XY = sum(devSmallXY4)/length(devSmallXY4);
avgDevRadius4XYZ = sum(devSmallXYZ4)/length(devSmallXYZ4);
avgDevCompViewRayXY = avgDevCompViewRayXY/length(sliceSphereX);
avgDevCompViewRayXYZ = avgDevCompViewRayXYZ/length(sliceSphereX);
percentPassXY1 = 100*countPassXY1/countSpheresInRadius1;
percentPassXYZ1 = 100*countPassXYZ1/countSpheresInRadius1;
percentPassXY2 = 100*countPassXY2/countSpheresInRadius2;
percentPassXYZ2 = 100*countPassXYZ2/countSpheresInRadius2;
percentPassXY3 = 100*countPassXY3/countSpheresInRadius3;
percentPassXYZ3 = 100*countPassXYZ3/countSpheresInRadius3;
percentPassXY4 = 100*countPassXY4/countSpheresInRadius4;
percentPassXYZ4 = 100*countPassXYZ4/countSpheresInRadius4;
maxDevRadius1XY = max(devSmallXY1);
maxDevRadius1XYZ = max(devSmallXYZ1);
maxDevRadius2XY = max(devSmallXY2);
maxDevRadius2XYZ = max(devSmallXYZ2);
maxDevRadius3XY = max(devSmallXY3);
maxDevRadius3XYZ = max(devSmallXYZ3);
maxDevRadius4XY = max(devSmallXY4);
maxDevRadius4XYZ = max(devSmallXYZ4);

% upscale the image for plotting
% original information
origX = 1:length(centSliceImg(1,:));
origY = 1:length(centSliceImg(:,1));
[origMeshX,origMeshY] = meshgrid(origX,origY);
% new data
volUpscaleX = 1/scaleFactor:1/scaleFactor:length(centSliceImg(1,:));
volUpscaleY = 1/scaleFactor:1/scaleFactor:length(centSliceImg(:,1));
[upMeshX, upMeshY] = meshgrid(volUpscaleX,volUpscaleY);
% the upscaled data
SigDataUpscale = interp2(origMeshX,origMeshY,centSliceImg,...
    upMeshX,upMeshY,'cubic');


if plotImage == 1;
    % Just the raw image
    figure;
    % fit marker plotting
    imshow(SigDataUpscale,[])
    % custom axes to show distances (the locations of the axes)
    numOfTicks = 10;
    xAxisSpacing = floor(length(SigDataUpscale(1,:))/numOfTicks);
    yAxisSpacing = floor(length(SigDataUpscale(:,1))/numOfTicks);
    
    % xAxisLocation = linspace(1,length(SigDataUpscale(1,:)),numOfTicks);
    % yAxisLocation = linspace(1,length(SigDataUpscale(:,1)),numOfTicks);
    xAxisLocation = 1:xAxisSpacing:(xAxisSpacing*numOfTicks+1);
    yAxisLocation = 1:yAxisSpacing:(yAxisSpacing*numOfTicks+1);
    % the values on the axes
    xAxisValue = voxelWidth/scaleFactor.*(xAxisLocation - round(median(xAxisLocation)));
    yAxisValue = voxelHeight/scaleFactor.*(yAxisLocation - round(median(yAxisLocation)));
    yAxisLocation(length(yAxisLocation)) = length(SigDataUpscale(:,1)); % special case for plotting for paper
    axis on
    axisHandle = gca;
    set(gca,'XTickMode','manual')
    set(gca,'YTickMode','manual')
    set(gca,'XTick',xAxisLocation)
    set(gca,'YTick',yAxisLocation)
    set(gca,'XTickLabel',xAxisValue)
    set(gca,'YTickLabel',yAxisValue)
    set(gca,'FontSize',16)
    xlabel('x-position (mm)','FontSize',20)
    ylabel('z-position (mm)','FontSize',20)
    
    % Image with passing and failing markers
    figure;
    % fit marker plotting
    imshow(SigDataUpscale,[])
    % cat(3,r,g,b)
    green = cat(3, zeros(size(sphereFitPtsPass)),ones(size(sphereFitPtsPass)), zeros(size(sphereFitPtsPass)));
    % bright = 0.9;
    % white version
    % green = bright.*cat(3, ones(size(sphereFitPtsPass)),ones(size(sphereFitPtsPass)), ones(size(sphereFitPtsPass)));
    hold on
    hGreen = imshow(green);
    hold off
    set(hGreen, 'AlphaData', sphereFitPtsPass) % make color sheet only show markers
    % the spheres that failed
    red = cat(3, ones(size(sphereFitPtsFail)),zeros(size(sphereFitPtsFail)), zeros(size(sphereFitPtsFail)));
    % white version
    % red = bright.*cat(3, ones(size(sphereFitPtsFail)),ones(size(sphereFitPtsFail)), ones(size(sphereFitPtsFail)));
    hold on
    hRed = imshow(red);
    hold off
    set(hRed, 'AlphaData', sphereFitPtsFail) % make color sheet only show markers
    
    % ground truth marker plotting
    % cat(3,r,g,b)
    rgb= [64,224,208]/255;
    white = cat(3, ones(size(grndTruthPts)),ones(size(grndTruthPts)),ones(size(grndTruthPts)));
    hold on
    hBlue = imshow(white);
    hold off
    set(hBlue, 'AlphaData', grndTruthPts) % make color sheet only show markers
    % title('Spatial Integrity Phantom Center Slice','FontSize',20)
    xlabel('x-position (mm)','FontSize',20)
    ylabel('z-position (mm)','FontSize',20)
    
    % custom axes to show distances (the locations of the axes)
    numOfTicks = 10;
    xAxisSpacing = floor(length(SigDataUpscale(1,:))/numOfTicks);
    yAxisSpacing = floor(length(SigDataUpscale(:,1))/numOfTicks);
    
    % xAxisLocation = linspace(1,length(SigDataUpscale(1,:)),numOfTicks);
    % yAxisLocation = linspace(1,length(SigDataUpscale(:,1)),numOfTicks);
    xAxisLocation = 1:xAxisSpacing:(xAxisSpacing*numOfTicks+1);
    yAxisLocation = 1:yAxisSpacing:(yAxisSpacing*numOfTicks+1);
    % the values on the axes
    xAxisValue = voxelWidth/scaleFactor.*(xAxisLocation - round(median(xAxisLocation)));
    yAxisValue = voxelHeight/scaleFactor.*(yAxisLocation - round(median(yAxisLocation)));
    yAxisLocation(length(yAxisLocation)) = length(SigDataUpscale(:,1)); % special case for plotting for paper
    axis on
    axisHandle = gca;
    set(gca,'XTickMode','manual')
    set(gca,'YTickMode','manual')
    set(gca,'XTick',xAxisLocation)
    set(gca,'YTick',yAxisLocation)
    set(gca,'XTickLabel',xAxisValue)
    set(gca,'YTickLabel',yAxisValue)
    set(gca,'FontSize',16)
    
end


end