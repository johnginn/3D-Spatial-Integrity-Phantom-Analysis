% Function to calculated the location of a sphere based on the weighted sum
% of signals within a given volume.
%
% Input:
% xSearch Array of x-locations used to calculate the lcoation of the sphere (native resolution of volume)
% ySearch Array of y-locations used to calculate the location of the sphere (native resolution of volume)
% zSearch Array of z-locations used to calculate the location of the sphere (native resolution of volume)
% sigVolume The volume containing the signals (must match x-,y-,z-Search)
% startWeight [x,y,z] The starting location for the search (native resolution)
% searchDistWeight The search distance in each direction for the weighted sum method
% voxelWidth The width of the voxel (x-dim)
% voxelHeight The height of the voxel (y-dim)
% voxelLength The length of the voxel (z-dim)
% volXDim The native volume x-dimension
% volYDim The native volume y-dimension
% volZDim The native volume z-dimension
% upscale Whether or not to upscale the weighted sum localization
% upscaleFactor The factor to upscale the image
%
% Output:
% weightCoord = [weightXInd, weightYInd, weightZInd];
%
% John Ginn
% Created: 7/6/16
% Modified: 1/24/17
% Edited to allow interpolation to non-voxel search distances


function [weightCoord] = ...
    weightSumLoc(xSearch,ySearch,zSearch,sigVolume,startWeight,searchDistWeight,...
    voxelWidth,voxelHeight,voxelLength,volXDim,volYDim,volZDim,upscale,upscaleFactor)

if upscale == 1
    [origMeshX,origMeshY,origMeshZ] = meshgrid(xSearch,ySearch,zSearch);
    % range in terms of upscale voxels
    xRangeWeightUp = (searchDistWeight/(voxelWidth/upscaleFactor));
    yRangeWeightUp = (searchDistWeight/(voxelHeight/upscaleFactor));
    zRangeWeightUp = (searchDistWeight/(voxelLength/upscaleFactor));
    % starting and ending location in terms of upscale voxels
    xSearchUpMin = startWeight(1)*upscaleFactor - xRangeWeightUp;
    xSearchUpMax = startWeight(1)*upscaleFactor + xRangeWeightUp;
    ySearchUpMin = startWeight(2)*upscaleFactor - yRangeWeightUp;
    ySearchUpMax = startWeight(2)*upscaleFactor + yRangeWeightUp;
    zSearchUpMin = startWeight(3)*upscaleFactor - zRangeWeightUp;
    zSearchUpMax = startWeight(3)*upscaleFactor + zRangeWeightUp;
    % convert to native resolution positions
    xSearchUpMin = xSearchUpMin./upscaleFactor;
    xSearchUpMax = xSearchUpMax./upscaleFactor;
    ySearchUpMin = ySearchUpMin./upscaleFactor;
    ySearchUpMax = ySearchUpMax./upscaleFactor;
    zSearchUpMin = zSearchUpMin./upscaleFactor;
    zSearchUpMax = zSearchUpMax./upscaleFactor;
    % check to make sure within the appropriate bounds
    if xSearchUpMin < 1
        xSearchUpMin = 1;
    end
    if xSearchUpMax > volXDim
        xSearchUpMax = volXDim;
    end
    if xSearchUpMax < 1
        xSearchUpMax = 1;
    end
    if ySearchUpMin < 1
        ySearchUpMin = 1;
    end
    if ySearchUpMax > volYDim
        ySearchUpMax = volYDim;
    end
    if ySearchUpMax < 1
        ySearchUpMax = 1;
    end
    if zSearchUpMin < 1
        zSearchUpMin = 1;
    end
    if zSearchUpMax > volZDim
        zSearchUpMax = volZDim;
    end
    if zSearchUpMax < 1
        zSearchUpMax = 1;
    end
    if xSearchUpMin > xSearchUpMax
        xSearchUpMin = xSearchUpMax;
    end
    if ySearchUpMin > ySearchUpMax
        ySearchUpMin = ySearchUpMax;
    end
    if zSearchUpMin > zSearchUpMax
        zSearchUpMin = zSearchUpMax;
    end
    % make sure upscale range is within the non-upscale range
    if xSearchUpMin < min(xSearch)
        xSearchUpMin = min(xSearch);
    end
    if xSearchUpMax > max(xSearch)
        xSearchUpMax = max(xSearch);
    end
    if ySearchUpMin < min(ySearch)
        ySearchUpMin = min(ySearch);
    end
    if ySearchUpMax > max(ySearch)
        ySearchUpMax = max(ySearch);
    end
    if zSearchUpMin < min(zSearch)
        zSearchUpMin = min(zSearch);
    end
    if zSearchUpMax > max(zSearch)
        zSearchUpMax = max(zSearch);
    end
    % array of upscale voxels converted back to native resolution positions
    % to select out the appropriate interpolation locations
    volUpscaleX = (xSearchUpMin:1/upscaleFactor:xSearchUpMax);
    volUpscaleY = (ySearchUpMin:1/upscaleFactor:ySearchUpMax);
    volUpscaleZ = (zSearchUpMin:1/upscaleFactor:zSearchUpMax);
    [upMeshX, upMeshY, upMeshZ] = meshgrid(volUpscaleX,volUpscaleY,volUpscaleZ);
    sigVolumeUp = interp3(origMeshX,origMeshY,origMeshZ,sigVolume,...
        upMeshX,upMeshY,upMeshZ,'cubic');
    % calculate the total signal
    totSig = 0;
    for xStep = 1:length(volUpscaleX)
        for yStep = 1:length(volUpscaleY)
            for zStep = 1:length(volUpscaleZ)
                totSig = totSig + sigVolumeUp(yStep,xStep,zStep);
            end
        end
    end
    
    % normalize the signal
    normSig = sigVolumeUp./totSig;
    weightXInd = 0;
    weightYInd = 0;
    weightZInd = 0;
    for xStep = 1:length(volUpscaleX)
        for yStep = 1:length(volUpscaleY)
            for zStep = 1:length(volUpscaleZ)
                currentSig = normSig(yStep,xStep,zStep);
                % weight the locations by the signal
                weightXInd = weightXInd + currentSig*volUpscaleX(xStep);
                weightYInd = weightYInd + currentSig*volUpscaleY(yStep);
                weightZInd = weightZInd + currentSig*volUpscaleZ(zStep);
            end
        end
    end
    
    weightCoord = [weightXInd, weightYInd, weightZInd];
    
else
    
    % calculate the total signal
    totSig = 0;
    for xStep = 1:length(xSearch)
        for yStep = 1:length(ySearch)
            for zStep = 1:length(zSearch)
                totSig = totSig + sigVolume(yStep,xStep,zStep);
            end
        end
    end
    
    % normalize the signal
    normSig = sigVolume./totSig;
    weightXInd = 0;
    weightYInd = 0;
    weightZInd = 0;
    for xStep = 1:length(xSearch)
        for yStep = 1:length(ySearch)
            for zStep = 1:length(zSearch)
                currentSig = normSig(yStep,xStep,zStep);
                % weight the locations by the signal
                weightXInd = weightXInd + currentSig*xSearch(xStep);
                weightYInd = weightYInd + currentSig*ySearch(yStep);
                weightZInd = weightZInd + currentSig*zSearch(zStep);
            end
        end
    end
    
    weightCoord = [weightXInd, weightYInd, weightZInd];
end

end