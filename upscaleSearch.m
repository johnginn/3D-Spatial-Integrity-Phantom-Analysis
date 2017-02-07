% This function takes the locations of the spheres that were found using
% the native resolution images, upscales the volume data and extracts a
% more precise location of the spheres.
%
% Input:
% volSigData The volume of signal data (dimensions: y,x,z)
% centerSearch [xInd,yInd,zInd] The center of the search volume
% voxelHeight (mm/voxel) The height of the pixel
% voxelWidth (mm/voxel) The width of the pixel
% voxelLength (mm/voxel) The length of the pixel
% sphereBoundary  y = 1, n = 0 Whether or not a boundary of 0's surrounds the sphere model
% MRorCT ('mr' or 'ct') Determines whether the template contrast exhibits a CT or MR scan 
%
% Output:
% coordSphereUp [xInd,yInd,zInd] The coordinates of the sphere from upscaled
% analysis as the appear to user in imshow
%
% optCorrelationUp The optimal correlation coefficient from upscaled
% analysis
%
% NOTE: The input index of the x and y corrdinates for centerSearch are
% reversed for this function to findSphereCenter.m This is due to the way
% that interp3 works and plotting in imshow is indexed.
%
% John Ginn
% Created: 6/29/16
% Modified: 11/3/16

function [coordSphereDataUp, optCorrelationUp] =...
    upscaleSearch(volSigData,centerSearch,voxelHeight,voxelWidth, voxelLength,sphereBoundary,MRorCT)
troubleShoot = 0; % troubleshooting? y = 1, n = 0

sphereDiameter = 8; % (mm)
% (mm) distance in each direction to search for the sphere. We changed this
% value to search only within one voxel of the first location determined by
% template matching at the native resolution of the images
searchDistX = voxelWidth; 
searchDistY = voxelHeight; 
searchDistZ = voxelLength; 

upscaleFactor = 3;
upscaleVoxX = voxelWidth/upscaleFactor; % (mm) the desired upscale voxel x-dimension
upscaleVoxY = voxelHeight/upscaleFactor; % (mm) the desired upscale voxel y-dimension
upscaleVoxZ = voxelLength/upscaleFactor; % (mm) the desired upscale voxel z-dimension
scaleFactorX = upscaleFactor;
scaleFactorY = upscaleFactor;
scaleFactorZ = upscaleFactor;


% calculate dimensions of upscaled data (calculating entire volume would
% take too long and require more memory)
volXDimOrig = length(volSigData(1,:,1)); % before scaling
volYDimOrig = length(volSigData(:,1,1)); % before scaling
volZDimOrig = length(volSigData(1,1,:)); % before scaling
volXDim = scaleFactorX.*volXDimOrig; % after scaling
volYDim = scaleFactorY.*volYDimOrig; % after scaling
volZDim = scaleFactorZ.*volZDimOrig; % after scaling

% find necessary range to sample
xRangeOrig = round(searchDistX/voxelWidth); % before scaling
yRangeOrig = round(searchDistY/voxelHeight); % before scaling
zRangeOrig = round(searchDistZ/voxelLength); % before scaling
xRange = round(searchDistX/upscaleVoxX); % after scaling
yRange = round(searchDistY/upscaleVoxY); % after scaling
zRange = round(searchDistZ/upscaleVoxZ); % after scaling

% the upscaled sphere volume makeSphere(height,width,length,diameter)
[sphereVolUpscale] = ...
    makeSphere(upscaleVoxY, upscaleVoxX, upscaleVoxZ, sphereDiameter,sphereBoundary,MRorCT);
[sphereVolOrig] = ...
    makeSphere(voxelHeight, voxelWidth, voxelLength, sphereDiameter,sphereBoundary,MRorCT);

% range of the sphere (will subtract from the range needed to be tested because
% the sphere itself occupies some volume within the search area)
% upscale
sphereXWidth = ceil(length(sphereVolUpscale(1,:,1))/2);
sphereYWidth = ceil(length(sphereVolUpscale(:,1,1))/2);
sphereZWidth = ceil(length(sphereVolUpscale(1,1,:))/2);
sphereXDim = length(sphereVolUpscale(1,:,1));
sphereYDim = length(sphereVolUpscale(:,1,1));
sphereZDim = length(sphereVolUpscale(1,1,:));
% original
sphereXWidthOrig = ceil(length(sphereVolOrig(1,:,1))/2);
sphereYWidthOrig = ceil(length(sphereVolOrig(:,1,1))/2);
sphereZWidthOrig = ceil(length(sphereVolOrig(1,1,:))/2);
sphereXDimOrig = length(sphereVolOrig(1,:,1));
sphereYDimOrig = length(sphereVolOrig(:,1,1));
sphereZDimOrig = length(sphereVolOrig(1,1,:));

% the location of the upscaled center search region
upscCenterSearch(1) = centerSearch(1).*scaleFactorX;
upscCenterSearch(2) = centerSearch(2).*scaleFactorY;
upscCenterSearch(3) = centerSearch(3).*scaleFactorZ;

% normalize sphere to signal in the entire region that will be sampled
normX = (centerSearch(1) - xRangeOrig - sphereXWidthOrig):1:...
    (centerSearch(1) + xRangeOrig + sphereXWidthOrig);
normY = (centerSearch(2) - yRangeOrig - sphereYWidthOrig):1:...
    (centerSearch(2) + yRangeOrig + sphereYWidthOrig);
normZ = (centerSearch(3) - zRangeOrig - sphereZWidthOrig):1:...
(centerSearch(3) + zRangeOrig + sphereZWidthOrig);
% avoid going beyond bounds of image
if max(normX) > volXDimOrig
    normX = (centerSearch(1) - xRangeOrig):1:volXDimOrig;
end
if max(normY) > volYDimOrig
    normY = (centerSearch(2) - yRangeOrig):1:volYDimOrig;
end
if max(normZ) > volZDimOrig
    normZ = (centerSearch(3) - zRangeOrig):1:volZDimOrig;
end
if min(normX) < 1
    normX = 1:1:(centerSearch(1) + xRangeOrig);
end
if min(normY) < 1
    normY = 1:1:(centerSearch(2) + yRangeOrig);
end
if min(normZ) < 1
    normZ = 1:1:(centerSearch(3) + zRangeOrig);
end
% normalize sphere to signal in the volume to be searched
normFactor = max(max(max(volSigData(normY,normX,normZ))));
normSphere = sphereVolUpscale.*normFactor;

% the bounds will be determined by the center location, the volume that
% will be searched, and the volume of the sphere
xMin = (upscCenterSearch(1) - xRange - sphereXWidth);
xMax = (upscCenterSearch(1) + xRange - sphereXWidth);
yMin = (upscCenterSearch(2) - yRange - sphereYWidth);
yMax = (upscCenterSearch(2) + yRange - sphereYWidth);
zMin = (upscCenterSearch(3) - zRange - sphereZWidth);
zMax = (upscCenterSearch(3) + zRange - sphereZWidth);

% make sure to stay within bounds of the image and avoid errors
if xMax > (volXDim - sphereXDim + 1)
    xMax = volXDim - sphereXDim + 1; % avoid going out of image volume
end
if xMin < 1
    xMin = 1; % avoid going out of image volume
elseif xMin > xMax
    xMin = xMax;
end
if yMax > (volYDim - sphereYDim + 1)
    yMax = volYDim - sphereYDim + 1; % avoid going out of image volume
end
if yMin < 1
    yMin = 1; % avoid going out of image volume
elseif yMin > yMax
    yMin = yMax;
end
if zMax > (volZDim - sphereZDim + 1)
    zMax = volZDim - sphereZDim + 1; % avoid going out of image volume
end
if zMin < 1
    zMin = 1; % avoid going out of image volume
elseif zMin > zMax
    zMin = zMax;
end


% calculate the upscaled volume data that will be searched
%
% these bounds are the bounds of the whole subvolume, need to obtain
% addtional region to account for width of sphere
xMinOrig = (centerSearch(1) - xRangeOrig - sphereXWidthOrig); 
xMaxOrig = (centerSearch(1) + xRangeOrig + sphereXWidthOrig);
yMinOrig = (centerSearch(2) - yRangeOrig - sphereYWidthOrig);
yMaxOrig = (centerSearch(2) + yRangeOrig + sphereYWidthOrig);
zMinOrig = (centerSearch(3) - zRangeOrig - sphereZWidthOrig);
zMaxOrig = (centerSearch(3) + zRangeOrig + sphereZWidthOrig); 
% avoid errors going outside the bounds
if xMaxOrig > (volXDimOrig)
    xMaxOrig = volXDimOrig; % avoid going out of image volume
end
if xMinOrig < 1
    xMinOrig = 1; % avoid going out of image volume
elseif xMinOrig > xMaxOrig
    xMinOrig = xMaxOrig;
end
if yMaxOrig > (volYDimOrig)
    yMaxOrig = volYDimOrig; % avoid going out of image volume
end
if yMinOrig < 1
    yMinOrig = 1; % avoid going out of image volume
elseif yMinOrig > yMaxOrig
    yMinOrig = yMaxOrig;
end
if zMaxOrig > (volZDimOrig)
    zMaxOrig = volZDimOrig; % avoid going out of image volume
end
if zMinOrig < 1
    zMinOrig = 1; % avoid going out of image volume
elseif zMinOrig > zMaxOrig
    zMinOrig = zMaxOrig;
end


%% Upscale the data
origX = xMinOrig:1:xMaxOrig;
origY = yMinOrig:1:yMaxOrig;
origZ = zMinOrig:1:zMaxOrig;
[origMeshX,origMeshY,origMeshZ] = meshgrid(origX,origY,origZ);

% % select out data to be upscaled (store y,x,z like in image data)
volSigOrig = zeros(length(origY),length(origX),length(origZ));
countX = 1;
countY = 1;
countZ = 1;
for stepX = origX
    countY = 1; % reset the count
    for stepY = origY       
        countZ = 1; % reset the count
        for stepZ = origZ
            % the selected data
            volSigOrig(countY,countX,countZ) = volSigData(stepY,stepX,stepZ);
            countZ = countZ + 1;
        end
        countY = countY + 1;
    end
    countX = countX + 1;
end

% note: do not start at 1/scaleFactor because original data starts at 1 
% starting less than 1 would sample outside of bounds of the volume selected
%
% these are the locations of the interpolation, NOT indexes thus the range
% will be larger than expected by a factor of the width of the sphere
volUpscaleX = xMinOrig:1/scaleFactorX:xMaxOrig;
volUpscaleY = yMinOrig:1/scaleFactorY:yMaxOrig;
volUpscaleZ = zMinOrig:1/scaleFactorZ:zMaxOrig;
[upMeshX, upMeshY, upMeshZ] = meshgrid(volUpscaleX,volUpscaleY,volUpscaleZ);
% the upscaled data
volSigDataUpscale = interp3(origMeshX,origMeshY,origMeshZ,volSigOrig,...
    upMeshX,upMeshY,upMeshZ,'cubic');


% cycle through the volume and calculate correlation coefficients -->
% start beginning correlation as max correlation --> cycle through and
% replace correlation if corrcoef yeilds a higher value

% initialize loop
optCorrelationUp = 0; % optimal correlation
coordSphereDataUp = centerSearch; % optimal coordinates
currentVol = zeros(sphereXDim,sphereYDim,sphereZDim);
% will return initial guess w/ optCorrelation = 0 if corrcoeff = nan
xCoord = centerSearch(1); 
yCoord = centerSearch(2);
zCoord = centerSearch(3);
% plus one on bounds because otherwise last volume location will be skipped
for X = 1:1:(length(volSigDataUpscale(1,:,1)) - sphereXDim + 1);
    % x-indices to search
    xSearch = (X:1:(X+sphereXDim - 1));
    for Y = 1:1:(length(volSigDataUpscale(:,1,1)) - sphereYDim + 1);
        % y-indices to search
        ySearch = (Y:1:(Y+sphereYDim - 1));
        for Z = 1:1:(length(volSigDataUpscale(1,1,:)) - sphereZDim + 1);
            % z-indices to search      
            zSearch = (Z:1:(Z+sphereZDim - 1));
            % calculate the correlation between the sphere and the current
            % volume selected in the loop (minus 1 otherwise dimensions
            % will not agree)
            currentVol = volSigDataUpscale(ySearch,xSearch,zSearch);
%             currentCorrelation = corrcoef(normSphere,currentVol);
%             if (currentCorrelation(2,1) > optCorrelationUp)
%                 % a new optimal correlation, update sphere location
%                 optCorrelationUp = currentCorrelation(2,1);
            currentCorrelation = corrCoeff3D(normSphere,currentVol);
            if (currentCorrelation > optCorrelationUp)
                % a new optimal correlation, update sphere location
                optCorrelationUp = currentCorrelation;
                % the corrdinates are at the median of the region sampled
                xCoordIndex = median(xSearch);
                yCoordIndex = median(ySearch); 
                zCoordIndex = median(zSearch);
                % because the upscaled sub-volume indexing is different from
                % original volume's index, grab optimal location from meshgrid
                % Note that coordinates are reversed because of meshgrid
                xCoord = upMeshX(1,xCoordIndex,1);
                yCoord = upMeshY(yCoordIndex,1,1);
                zCoord = upMeshZ(1,1,zCoordIndex); 
            end
        end
    end
end
% reversed compared to findSphereCenter
coordSphereDataUp = [yCoord,xCoord,zCoord]; % The coordinates of the sphere as they appear to user in image



%% troubleshooting plots and GUI
if troubleShoot == 1;
    % check the upscale algorithm
    guiData{1} = volSigOrig;
    guiData{2} = volSigDataUpscale;
    guiData{3} = scaleFactorZ;
    UpscaleGUI(guiData)
    % plot the results
end
end