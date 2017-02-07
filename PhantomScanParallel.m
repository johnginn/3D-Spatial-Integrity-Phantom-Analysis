% This script cycles through the different spheres in the spatial integrity
% phantom in order to calculate their location. This function utilizes
% MATLAB's parallel toolbox to speed up the analysis
%
% Input:
% volSigData The volume of signal data (dimensions: x,y,z)
% sphereData The volume of the sphere (dimensions: x,y,z)
% voxelHeight (mm/voxel) The height of the pixel
% voxelWidth (mm/voxel) The width of the pixel
% voxelLength (mm/voxel) The length of the pixel
% slices The array of slice indices
% centerCol The index of the center column
% centerRow The index of the center row
% centerSlice The index of the center slice
% slicesGnd The array of ground truth slices
% nRows The number of rows
% nSlices The number of slices
% spheresInRow Number of spheres in each row from bottom to top, for SI phantom [29 29 29 29 29 29 29 29 27 27 27 27 25 25 23 23 21 19 17 13 7]
% centerRowNum The number of the center row (not index, ie. 7 rows --> 4)
% searchDist (mm) distance in each direction to search for the sphere
% searchDistWeight (mm) distance in each direction from the location
% sphereBoundary  y = 1, n = 0 Whether or not a boundary of 0's surrounds the sphere model
% upScale y = 1, n = 0 search for the sphere center based off of upscaled data
% MRorCT ('mr' or 'ct') Determines whether the template contrast exhibits a CT or MR scan 
% xSpacing (mm) The x-spacing distance between the spheres
% ySpacing (mm) The y-spacing distance between the spheres
% zSpacing (mm) The z-spacing distance between the spheres
%
% Output:
% coordSphereData [xPos yPos zPos] The position of the sphere located by program 
% coordSpherePlot [yPos xPos zPos] The position of the sphere located by program 
% optCorrelation  The optimal correlation coefficient for each sphere
% sphereGoundTruth [yPos xPos zPos] The starting location of the search area (taken to be
% the ground truth)
% weightCoord [xPos yPos zPos] The position of the sphere by the weighted sum method
% weightCoordPlot [yPos xPos zPos] The position of the sphere by the weighted sum method
%
% coordSphereDataUp [xInd,yInd,zInd] The coordinates of the sphere from upscaled
% analysis as the appear to user in imshow
% 
% coordSpherePlotUp [yInd,xInd,zInd] Nearst voxel in original image using
% the upscaled sphere location
%
% optCorrelationUp The optimal correlation coefficient from upscaled
% analysis
%
% volSearchRegion Volume of 1's and 0's showing the search region used in
% the analysis
%
% volAllGndTruth Volume of 1's and 0's forming + shape for checking ground
% truth locations
%
% volWeightSum Volume of 1's and 0's forming + shape for checking weighted
% sum calculation regions
%
% John Ginn
% Created: 10/20/16
% Modified: 12/13/16

function [coordSphereData, coordSpherePlot, optCorrelation, sphereGoundTruth,...
    weightCoord, weightCoordPlot, coordSphereDataUp,coordSpherePlotUp, optCorrelationUp,...
    volSearchRegion,volAllGndTruth,volWeightSum] = ...
    PhantomScanParallel(volSigData, sphereData, voxelHeight,voxelWidth, voxelLength,slices,...
    centerCol,centerRow,centerSlice,slicesGnd,nRows,nSlices,spheresInRow,centerRowNum,searchDist,...
    searchDistWeight,upScale,sphereBoundary,MRorCT,xSpacing,ySpacing,zSpacing)

% phantom information
distBtwnSpheres = 16; % (mm) The distance between the spheres in each direction
xIndPerSph = xSpacing/voxelWidth; % number of pixels in x-direction between the spheres
yIndPerSph = ySpacing/voxelHeight; % number of pixels in y-direction between the spheres
zIndPerSph = zSpacing/voxelLength; % number of pixels in z-direction between the spheres
searchDistWeightSmall = 4; % (mm) the weighed sum calculation distance for spheres near the wall of the phantom


% calculate the locations of the center sphere and six surrounding spheres
surSphere(1,:) = [centerRow centerCol centerSlice];
surSphere(2,:) = [(centerRow-yIndPerSph) centerCol centerSlice];
surSphere(3,:) = [(centerRow+yIndPerSph) centerCol centerSlice];
surSphere(4,:) = [centerRow (centerCol-xIndPerSph) centerSlice];
surSphere(5,:) = [centerRow (centerCol+xIndPerSph) centerSlice];
surSphere(6,:) = [centerRow centerCol (centerSlice-zIndPerSph)];
surSphere(7,:) = [centerRow centerCol (centerSlice+zIndPerSph)];
countFound = 0; % used for debugging to calculate the number of spheres found
nTotData = nSlices.*sum(spheresInRow);

% determine the progress of finding the spheres
percentFound = 10:10:100;
numSpheresPercentFound = floor(nTotData./10.*(1:10)); % array of number of spheres corresponding to percent found
countPercentFound = 1;

% initialize data for efficiency
coordSphereData = cell(1,nTotData);
coordSpherePlot = cell(1,nTotData);
optCorrelation = cell(1,nTotData);
sphereGoundTruth = cell(1,nTotData);
coordSphereDataUp = cell(1,nTotData);
coordSpherePlotUp = cell(1,nTotData);
optCorrelationUp = cell(1,nTotData);
weightCoord  = cell(1,nTotData);
weightCoordPlot  = cell(1,nTotData);
searchIndX = cell(1,nTotData);
searchIndY = cell(1,nTotData);
searchIndZ = cell(1,nTotData);
% array of the weighted sum distances used to calculate the locations
currWeightSumRange = searchDistWeight.*ones(1,nTotData);
% cycle through slices
yDimVol = length(volSigData(:,1,1)); % dimension of the y-volume
saveCount = 1; % used to save the data
for sliceStep = 1:nSlices
    % the current slice
    currentSlice = slices(sliceStep);
    gndTruthSlice = slicesGnd(sliceStep);
    % step through the rows
    for rowStep = 1:nRows
            % move upward, inverted because of indexing in plot
            % the count should start at 1 and move up (thus subtract 11)
            rowPos = (rowStep - centerRowNum);
            currentRow = centerRow - round(yIndPerSph*rowPos);
            gndTruthRow = centerRow - yIndPerSph*rowPos;
        % spheres remaining to step through in the row
        spheresRemaining = spheresInRow(rowStep); 
        stepRight = 0; % reset the step right counter
        % step through the columns
        currentCol = centerCol; % start at the center
        for colStep = 1:spheresInRow(rowStep)
            % find the center sphere
            if spheresRemaining == (spheresInRow(rowStep))
                % do nothing since aready at the center location
                 gndTruthCol = centerCol;
                if (rowStep == 1)
                    % at the lowest row that borders the phantom base
                    currWeightSumRange(saveCount) = searchDistWeightSmall;
                end
            % move through spheres to the left 
            elseif (spheresRemaining > (spheresInRow(rowStep) - 1)/2 &&...
                spheresRemaining ~= (spheresInRow(rowStep)))
                stepAmount = colStep - 1; % remove one since the first step 
                % was used to calculate the center location
                % the current column location
                currentCol = centerCol - round(xIndPerSph*stepAmount);
                gndTruthCol = centerCol - xIndPerSph*stepAmount;
                % 
                % change range of weighted sum for spheres near the wall
                if (spheresInRow(rowStep) == 29)&&(stepAmount==((spheresInRow(rowStep)-1)/2))
                    % at the boundary in the lowest 8 rows
                    currWeightSumRange(saveCount) = searchDistWeightSmall;
                end
                if (rowStep == 1)
                    % at the lowest row that borders the phantom base
                    currWeightSumRange(saveCount) = searchDistWeightSmall;
                end
                if ((nRows - rowStep)<= 5)&&(stepAmount==((spheresInRow(rowStep)-1)/2))
                    % at the top six rows, sphere furthest from isocenter
                    currWeightSumRange(saveCount) = searchDistWeightSmall;
                end
            else
                % the current column location
                stepRight = stepRight + 1;
                currentCol = centerCol + round(xIndPerSph*stepRight);  
                gndTruthCol = centerCol + xIndPerSph*stepRight;
                % 
                % change range of weighted sum for spheres near the wall
                if (spheresInRow(rowStep) == 29)&&(stepRight==14)
                    % at the boundary in the lowest rows
                    currWeightSumRange(saveCount) = searchDistWeightSmall;
                end
                if (rowStep == 1)
                    % at the lowest row that borders the phantom base
                    currWeightSumRange(saveCount) = searchDistWeightSmall;
                end
                if ((nRows - rowStep)<= 5)&&(stepRight==((spheresInRow(rowStep)-1)/2))
                    % at the top six rows, sphere furthest from isocenter
                    currWeightSumRange(saveCount) = searchDistWeightSmall;
                end
            end
            % check to see if the current sphere is either the center
            % sphere or one of the 6 surrounding spheres
            % calculate the locations of the center sphere and six surrounding spheres
            %
            % set the weighted sum range to 4 mm if so
            for stepSurSphere = 1:length(surSphere)
                if (surSphere(stepSurSphere,1) == gndTruthRow)&&...
                        (surSphere(stepSurSphere,2) == gndTruthCol)&&...
                        (surSphere(stepSurSphere,3) == gndTruthSlice)
                        currWeightSumRange(saveCount) = searchDistWeightSmall; % 4 mm
                        countFound = countFound + 1;
                end
            end
            startLocation(saveCount,:) = [currentCol currentRow currentSlice];
            % ground truth locations (note rotation for plotting) [y x z]
            sphereGoundTruth{saveCount} = [gndTruthRow gndTruthCol gndTruthSlice];
            saveCount = saveCount + 1; % used to save the data
            % one more sphere has been found
            spheresRemaining = spheresRemaining - 1;
        end
    end
end

% countSteps = [];
parfor stepParallel = 1:nTotData
    %     countSteps = [countSteps,1];
    %     if countSteps == ceil(nTotData/4)
    %         disp('25% of sphere locations determined')
    %     end
    % Find the sphere
    [coordSphereData{stepParallel},coordSpherePlot{stepParallel}, optCorrelation{stepParallel},...
        weightCoord{stepParallel}, weightCoordPlot{stepParallel},...
        coordSphereDataUp{stepParallel}, coordSpherePlotUp{stepParallel}, optCorrelationUp{stepParallel},...
        searchIndX{stepParallel},searchIndY{stepParallel},searchIndZ{stepParallel}] = ...
        findSphereCenter(volSigData, sphereData,startLocation(stepParallel,:),voxelHeight,voxelWidth, voxelLength,...
        searchDist,currWeightSumRange(stepParallel),upScale,sphereBoundary,MRorCT);
end


% used for troubleshooting
% for x = 1:length(coordSpherePlot)
% testSphere(x,:) = weightCoordPlot{x};
% testSphereX(x) = testSphere(x,1);
% testSphereY(x) = testSphere(x,2);
% testSphereZ(x) = testSphere(x,3);
% testGnd(x,:) = sphereGoundTruth{x};
% testGndX(x) = testGnd(x,1);
% testGndY(x) = testGnd(x,2);
% testGndZ(x) = testGnd(x,3);
% end
% plotSym3D(testSphereX,testSphereY,testSphereZ,...
%     testGndX,testGndY,testGndZ);

% initialize volume of data for search region using correlation coeff
volDimX = length(volSigData(:,1,1));
volDimY = length(volSigData(1,:,1));
volDimZ = length(volSigData(1,1,:));
volSearchRegion = zeros(volDimX,volDimY,volDimZ);
volAllGndTruth = volSearchRegion;
volWeightSum = volSearchRegion;
% step through each of the spheres and construct a border to show search
% region
for step = 1:(nTotData - 1)
    % make volume of + shapes for plotting
    gndTruthCurrent = round(sphereGoundTruth{step});
    gndTruthX = (gndTruthCurrent(1)-1):1:(gndTruthCurrent(1)+1);
    gndTruthY = (gndTruthCurrent(2)-1):1:(gndTruthCurrent(2)+1);
    % no need to rotate x,y because that was already done
    volAllGndTruth(gndTruthX,gndTruthCurrent(2),gndTruthCurrent(3)) = 1;
    volAllGndTruth(gndTruthCurrent(1),gndTruthY,gndTruthCurrent(3)) = 1;
    % search region by correlation coefficent
    xRegion = searchIndX{step};
    xMin = min(xRegion);
    xMax = max(xRegion);
    yRegion = searchIndY{step};
    yMin = min(yRegion);
    yMax = max(yRegion);
    zRegion = searchIndZ{step};
    zMin = min(zRegion);
    zMax = max(zRegion);
    for stepX = xRegion
        for stepY = yRegion
            for stepZ = zRegion
                % if at one of the boundaries of the search region store a
                % 1 for plotting
                if (((stepX == xMin) || (stepX == xMax))||...
                        ((stepY == yMin) || (stepY == yMax)))
                    % ensure integer values
                    stepX = round(stepX);
                    stepY = round(stepY);
                    stepZ = round(stepZ);
                    % reverse for plotting
                    volSearchRegion(stepY,stepX,stepZ) = 1;
                end
            end
        end
    end
    
    % the location found by the correlation coefficent for the current
    % sphere
    sphereCurrent = round(coordSpherePlot{step});
    % weighted sum region calculation
    xWeightDist = round(currWeightSumRange(step)/voxelWidth);
    yWeightDist = round(currWeightSumRange(step)/voxelHeight);
    zWeightDist = round(currWeightSumRange(step)/voxelLength);
    % x-bounds
    if (sphereCurrent(1)-xWeightDist) < 1;
        xMinWeight = 1;
    else
        xMinWeight = (sphereCurrent(1)-xWeightDist);
    end
    if (sphereCurrent(1)+xWeightDist) > volDimX;
        xMaxWeight = xMax;
    else
        xMaxWeight = (sphereCurrent(1)+xWeightDist); % already ensured this does not go beyond boundary
    end
    % y-bounds
    if (sphereCurrent(2)-yWeightDist) < 1;
        yMinWeight = 1;
    else
        yMinWeight = (sphereCurrent(2)-yWeightDist);
    end
    if (sphereCurrent(2)+yWeightDist) > volDimY;
        yMaxWeight = yMax;
    else
        yMaxWeight = (sphereCurrent(2)+yWeightDist); % already ensured this does not go beyond boundary
    end
    
    if (sphereCurrent(3)-zWeightDist) < 1;
        zMinWeight = 1;
    else
        zMinWeight = (sphereCurrent(3)-zWeightDist);
    end
    if (sphereCurrent(3)+zWeightDist) > volDimZ;
        zMaxWeight = 1;
    else
        zMaxWeight = (sphereCurrent(3)+zWeightDist); % already ensured this does not go beyond boundary
    end
    % make sure min is not > max
    if xMinWeight > xMaxWeight
        xMinWeight = xMaxWeight;
    end
    if yMinWeight > yMaxWeight
        yMinWeight = yMaxWeight;
    end
    if zMinWeight > zMaxWeight
        zMinWeight = zMaxWeight;
    end
    
    xWeight = xMinWeight:1:xMaxWeight;
    yWeight = yMinWeight:1:yMaxWeight;
    zWeight = zMinWeight:1:zMaxWeight;
    for stepX = xWeight
        for stepY = yWeight
            for stepZ = zWeight
                % if at one of the boundaries of the search region store a
                % 1 for plotting
                if (((stepX == xMinWeight) || (stepX == xMaxWeight))||...
                        ((stepY == yMinWeight) || (stepY == yMaxWeight)))
                    % ensure integer values
                    stepX = round(stepX);
                    stepY = round(stepY);
                    stepZ = round(stepZ);
                    volWeightSum(stepX,stepY,stepZ) = 1;
                end
            end
        end
    end
end


end