% This function searches for the center of the sphere in the phantom given
% the volume data of the phantom an array containing the volume of the 
% sphere. 
%
% Input:
% volSigData The volume of signal data (dimensions: x,y,z)
% sphereData The volume of the sphere (dimensions: x,y,z)
% centerSearch [xInd,yInd,zInd] The center of the search volume
% voxelHeight (mm/voxel) The height of the pixel
% voxelWidth (mm/voxel) The width of the pixel
% voxelLength (mm/voxel) The length of the pixel
% searchDist (mm) distance in each direction to search for the sphere
% searchDistWeight (mm) distance in each direction from the location
% sphereBoundary  y = 1, n = 0 Whether or not a boundary of 0's surrounds the sphere model
% upScale y = 1, n = 0 search for the sphere center based off of upscaled data
% MRorCT ('mr' or 'ct') Determines whether the template contrast exhibits a CT or MR scan 
%
%
% Output:
% coordSphere [xInd,yInd,zInd] The coordinates of the sphere as the appear to user in imshow
% 
% coordSpherePlot [yInd,xInd,zInd] The coordinates of the sphere 
% for plotting because the coordinates are reversed in imshow
%
% optCorrelation The optimal correlation coefficient
%
% weightCoord [xInd,yInd,zInd] The coordinates determined by the weighted sum method
% 
% weightCoordPlot [yInd,xInd,zInd] The coordinates determined by the weighted sum method
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
% normX The x-indicies of the volume searched
% normY The y-indicies of the volume searched
% normZ The z-indicies of the volume searched
%
% John Ginn
% Created: 6/24/16
% Modified: 1/24/17
%
% Edited to allow interpolation to non-voxel search distances

function [coordSphereData,coordSpherePlot, optCorrelation,...
    weightCoord, weightCoordPlot,coordSphereDataUp, coordSpherePlotUp,...
    optCorrelationUp, normX, normY, normZ] = ...
    findSphereCenter(volSigData, sphereData,centerSearch,voxelHeight,voxelWidth, voxelLength,...
    searchDist,searchDistWeight,upScale,sphereBoundary,MRorCT)
troubleShoot = 0; % troubleshooting? y = 1, n = 0

% search volume (note: sphere separation is 16 mm in the phantom, sphere 
% sphere diameter is 8 mm)

% determined by the non-upscaled correlation coefficient search to calculate 
% the weigthed sum

% find necessary range to sample
xRange = round(searchDist/voxelWidth);
yRange = round(searchDist/voxelHeight);
zRange = round(searchDist/voxelLength);
xRangeWeight = round(searchDistWeight/voxelWidth);
yRangeWeight = round(searchDistWeight/voxelHeight);
zRangeWeight = round(searchDistWeight/voxelLength);


volXDim = length(volSigData(1,:,1));
volYDim = length(volSigData(:,1,1));
volZDim = length(volSigData(1,1,:));

% range of the sphere (will subtract from the range needed to be tested because
% the sphere itself occupies some volume within the search area)
sphereXWidth = ceil(length(sphereData(1,:,1))/2);
sphereYWidth = ceil(length(sphereData(:,1,1))/2);
sphereZWidth = ceil(length(sphereData(1,1,:))/2);
sphereXDim = length(sphereData(1,:,1));
sphereYDim = length(sphereData(:,1,1));
sphereZDim = length(sphereData(1,1,:));

% % normalize sphere to signal in the entire region that will be sampled
normX = (centerSearch(1) - xRange - sphereXWidth):1:(centerSearch(1) + xRange + sphereXWidth);
normY = (centerSearch(2) - yRange - sphereYWidth):1:(centerSearch(2) + yRange + sphereYWidth);
normZ = (centerSearch(3) - zRange - sphereZWidth):1:(centerSearch(3) + zRange + sphereZWidth);
% % avoid going beyond bounds of image
if max(normX) > volXDim
    normX = (centerSearch(1) - xRange):1:volXDim;
end
if max(normY) > volYDim
    normY = (centerSearch(2) - yRange):1:volYDim;
end
if max(normZ) > volZDim
    normZ = (centerSearch(3) - zRange):1:volZDim;
end
if min(normX) < 1
    normX = 1:1:(centerSearch(1) + xRange);
end
if min(normY) < 1
    normY = 1:1:(centerSearch(2) + yRange);
end
if min(normZ) < 1
    normZ = 1:1:(centerSearch(3) + zRange);
end

% normalize sphere to signal in the volume to be searched
% normFactor = max(max(max(volSigData(normY,normX,normZ))));
normFactor = 1; % no normalization to signal in the region
normSphere = sphereData.*normFactor;

% the bounds will be determined by the center location, the volume that
% will be searched, and the volume of the sphere
xMin = (centerSearch(1) - xRange - sphereXWidth);
xMax = (centerSearch(1) + xRange - sphereXWidth);
yMin = (centerSearch(2) - yRange - sphereYWidth);
yMax = (centerSearch(2) + yRange - sphereYWidth);
zMin = (centerSearch(3) - zRange - sphereZWidth);
zMax = (centerSearch(3) + zRange - sphereZWidth);

% make sure to stay within bounds of the image and avoid errors
if xMax > (volXDim - sphereXDim + 1)
    xMax = volXDim - sphereXDim + 1; % avoid going out of image volume
end
if xMax < 1
    xMax = 1;
end
if xMin < 1
    xMin = 1; % avoid going out of image volume
end
if yMax > (volYDim - sphereYDim + 1)
    yMax = volYDim - sphereYDim + 1; % avoid going out of image volume
end
if yMax < 1
    yMax = 1;
end
if yMin < 1
    yMin = 1; % avoid going out of image volume
end
if zMax > (volZDim - sphereZDim + 1)
    zMax = volZDim - sphereZDim + 1; % avoid going out of image volume
end
if zMax < 1
    zMax = 1;
end
if zMin < 1
    zMin = 1; % avoid going out of image volume
end
% make sure xMin not > xMax
if xMin > xMax
    xMin = xMax;
end
if yMin > yMax
    yMin = yMax;
end
if zMin > zMax
    zMin = zMax;
end

% make sure the bounds are integer values
xMin = round(xMin);
xMax = round(xMax);
yMin = round(yMin);
yMax = round(yMax);
zMin = round(zMin);
zMax = round(zMax);
% initialize loop
optCorrelation = 0; % optimal correlation
% will return initial guess w/ optCorrelation = 0 if corrcoeff = nan
xCoord = centerSearch(1); 
yCoord = centerSearch(2);
zCoord = centerSearch(3);
for X = xMin:1:xMax;
    % x-indices to search
    xSearch = (X:1:(X+sphereXDim - 1));
    for Y = yMin:1:yMax;
        % y-indices to search
        ySearch = (Y:1:(Y+sphereYDim - 1));
        for Z = zMin:1:zMax;
            % z-indices to search      
            zSearch = (Z:1:(Z+sphereZDim - 1));
            % calculate the correlation between the sphere and the current
            % volume selected in the loop (minus 1 otherwise dimensions
            % will not agree)
            % rotated because of way data is stored
            currentVol = volSigData(ySearch,xSearch,zSearch);
%             currentCorrelation = corrcoef(normSphere,currentVol);
%             if (currentCorrelation(2,1) > optCorrelation)
%                 % a new optimal correlation, update sphere location
%                 optCorrelation = currentCorrelation(2,1);
            currentCorrelation = corrCoeff3D(normSphere,currentVol);
            if (currentCorrelation > optCorrelation)
                % a new optimal correlation, update sphere location
                optCorrelation = currentCorrelation;
                % the corrdinates are at the median of the region sampled
                xCoord = median(xSearch); 
                yCoord = median(ySearch); 
                zCoord = median(zSearch);
                plotPhantomVol = currentVol; % store the data for plotting
                optXInd = xSearch; % optimal x-indices
                optYInd = ySearch; % optimal y-indices
                optZInd = zSearch; % optimal z-indices
            end
        end
    end
end
coordSpherePlot = [yCoord,xCoord,zCoord]; % The coordinates of the sphere for plotting 
coordSphereData = [xCoord,yCoord,zCoord]; % The coordinates of the sphere as they appear to user in image


% calculate location of center based on upscaled volume data%
%
% NOTE: The input index of the x and y corrdinates for centerSearch are
% reversed for this function to findSphereCenter.m This is due to the way
% that interp3 works and plotting in imshow is indexed.
if upScale == 1;
    [coordSpherePlotUp, optCorrelationUp] = ...
        upscaleSearch(volSigData,coordSphereData,voxelHeight,voxelWidth, voxelLength,sphereBoundary,MRorCT);
    % nearst voxel in original image using the upscaled sphere location
    coordSphereDataUp = ([coordSpherePlotUp(2),coordSpherePlotUp(1),coordSpherePlotUp(3)]);
else
    coordSphereDataUp = nan;
    optCorrelationUp = nan;
    coordSpherePlotUp = nan;
end

if upScale == 1;
    startWeight = coordSphereDataUp;
    upscaleFactor = 3;
    % use the ceil function here to ensure interpolation does not go "out
    % of bounds" for these voxels
    xRangeWeight = round(searchDistWeight/voxelWidth);
    yRangeWeight = round(searchDistWeight/voxelHeight);
    zRangeWeight = round(searchDistWeight/voxelLength);
else
    startWeight = coordSphereData;
end

% the weighted sum Method
xMinWeight = floor(startWeight(1) - xRangeWeight);
if xMinWeight < 1
   xMinWeight = 1; 
end

xMaxWeight = ceil(startWeight(1) + xRangeWeight);
if xMaxWeight > volXDim
    xMaxWeight = volXDim;
end
if xMaxWeight < 1
    xMaxWeight = 1;
end
yMinWeight = floor(startWeight(2) - yRangeWeight);
if yMinWeight < 1
   yMinWeight = 1; 
end
yMaxWeight = ceil(startWeight(2) + yRangeWeight);
if yMaxWeight > volYDim
    yMaxWeight = volYDim;
end
if yMaxWeight < 1
    yMaxWeight = 1;
end
zMinWeight = floor(startWeight(3) - zRangeWeight);
if zMinWeight < 1
   zMinWeight = 1; 
end
zMaxWeight = ceil(startWeight(3) + zRangeWeight);
if zMaxWeight > volZDim
    zMaxWeight = volZDim;
end
if zMaxWeight < 1
    zMaxWeight = 1;
end
% make sure min is not greater than max
if xMinWeight > xMaxWeight
    xMinWeight = xMaxWeight;
end
if yMinWeight > yMaxWeight
    yMinWeight = yMaxWeight;
end
if zMinWeight > zMaxWeight
    zMinWeight = zMaxWeight;
end
xSearchWeight = xMinWeight:1:xMaxWeight;
ySearchWeight = yMinWeight:1:yMaxWeight;
zSearchWeight = zMinWeight:1:zMaxWeight;
sigVolumeWeight = volSigData(ySearchWeight,xSearchWeight,zSearchWeight); % reversed data because of way it is stored
weightCoord = weightSumLoc(xSearchWeight,ySearchWeight,zSearchWeight,sigVolumeWeight,...
    startWeight,searchDistWeight,voxelWidth,voxelHeight,voxelLength,...
    volXDim,volYDim,volZDim,upScale,upscaleFactor);
weightCoordPlot = [weightCoord(2),weightCoord(1),weightCoord(3)];
%% troubleshooting plots and GUI
if troubleShoot == 1;
    % plot the results
    sphereCenterSlice = median(1:1:sphereZDim); % find the index of the center sphere slice
    figure;
    windHigh = normFactor; % normalize windowing
    windLow = 0; % normalize windowing
    subplot(3,1,1)
    imshow(plotPhantomVol(:,:,sphereCenterSlice),[windLow windHigh])
    title('plotPhantomVol Image')
    subplot(3,1,2)
    volData = volSigData(optYInd,optXInd,optZInd);
    imshow(volData(:,:,sphereCenterSlice),[windLow windHigh])
    title('Volume Data Image')
    subplot(3,1,3)
    imshow(normSphere(:,:,sphereCenterSlice),[windLow windHigh])
    title('Sphere Image')
    figure;
    imshow(volSigData(:,:,zCoord),[windLow windHigh])
    title(strcat(['Entire Volume, Slice:',num2str(zCoord)]),'FontSize',24)
end

end