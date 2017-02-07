% Function to calculate the SNR in the centermost sphere. Note, if the
% standard deviation of the background is 0, the code will shift the ROI
% upwards until the standard deviation of the background signal is no
% longer zero.
%
% Input:
% volData The volume data
% voxelHeight (mm) The height of the voxels
% voxelWidth (mm) The width of the voxels
% voxelLength (mm) The length of the voxels
% centerSlice The index for the center slice of spheres
% centerCol The location of center column of spheres
% centerRow The index of the center row
% sphereDiameter (mm) The diameter of the sphere
% xData The sphere x-location after rotation and translation corrections
% yData The sphere y-location after rotation and translation corrections
% zData The sphere z-location after rotation and translation corrections
% xGndTruth The ground-truth x-location
% yGndTruth The ground-truth y-location
% zGndTruth The ground-truth z-location
%
% Output:
% SNR The signal to noise ratio in the centermost sphere in the phantom
% avgSig The average signal in the center sphere volume
% bgStandardDev The standard deviation of the signal in the background volume
% 
% John Ginn
% Created: 8/16/16
% Modified: 8/16/16
function[SNR, avgSig, bgStandardDev] = calcSNR(volData,voxelHeight,voxelWidth,voxelLength,centerSlice,...
centerCol,centerRow,sphereDiameter,xData,yData,zData,xGndTruth,yGndTruth,zGndTruth)
plotSNRImg = 0;

% parameters for the background signal calculation
cubeLength = 10; % length of one of the sides of the cube use to calculate the
% noise (standard deviation of background signal)
xBgRange = round(cubeLength/voxelWidth);
yBgRange = round(cubeLength/voxelHeight);
zBgRange = round(cubeLength/voxelLength);

sphereRadius = sphereDiameter/2; % radius of each sphere
% range to calculate signal in center sphere
xRange = round(sphereRadius/voxelWidth);
yRange = round(sphereRadius/voxelHeight);
zRange = round(sphereRadius/voxelLength);

% dimensions of the volume
xDim = length(volData(1,:,1));
yDim = length(volData(:,1,1));
zDim = length(volData(1,1,:));


% find the location of the center sphere after the rotation and translation
% corrections have been applied
countCenter = 0;
for step = 1:length(xData)
    if (xGndTruth(step) == centerCol)&&...
            (yGndTruth(step) == centerRow)&&...
            (zGndTruth(step) == centerSlice)
        xCenterLoc = round(xData(step));
        yCenterLoc = round(yData(step));
        zCenterLoc = round(zData(step));
        countCenter = countCenter + 1;
    end
end

% calculate volume to determine average signal in the sphere
xMin = xCenterLoc - xRange;
xMax = xCenterLoc + xRange;
if xMin < 1;
    xMin = 1;
end
if xMax > xDim
    xMax = xDim;
end
if xMax < xMin
    xMax = xMin;
end
yMin = yCenterLoc - yRange;
yMax = yCenterLoc + yRange;
if yMin < 1;
    yMin = 1;
end
if yMax > yDim
    yMax = yDim;
end
if yMax < yMin
    yMax = yMin;
end
zMin = zCenterLoc - zRange;
zMax = zCenterLoc + zRange;
if zMin < 1;
    zMin = 1;
end
if zMax > zDim
    zMax = zDim;
end
if zMax < zMin
    zMax = zMin;
end



% image of the center slice
sliceImg = volData(:,:,centerSlice);
signalImg =zeros(length(sliceImg(:,1)),length(sliceImg(1,:)));

% calculate average signal in the center-most sphere
count = 0;
totSig = 0;
for x = xMin:1:xMax
    for y = yMin:1:yMax
        for z = zMin:1:zMax
            count = count + 1;
            totSig = totSig + volData(y,x,z);
            if (x == xMin)||(x == xMax)||(y == yMin)||(y == yMax)
                signalImg(y,x) = 1; % store the region used to calculate the signal
            end
        end
    end
end

% average signal in the sphere
avgSig = totSig/count;

% calculate the standard deviation of the signal in the background
bgStandardDev = 0;
% start at the center, and bottom of the image. Move background region up
% until the standard deviation is no longer equal to zero
shiftUp = 0; % number of voxels from the bottom of the image
while bgStandardDev == 0;
    shiftUp = shiftUp + 1;
    xBgMin = xCenterLoc - xRange;
    xBgMax = xCenterLoc + xRange;
    if xBgMin < 1;
        xBgMin = 1;
    end
    if xBgMax > xDim
        xBgMax = xDim;
    end
    if xBgMax < xBgMin
        xBgMax = xBgMin;
    end
    yBgMin = yDim - 2*yRange - shiftUp;
    yBgMax = yDim - shiftUp;
    if yBgMin < 1;
        yBgMin = 1;
    end
    if yBgMax > yDim
        yBgMax = yDim;
    end
    if yBgMax < yBgMin
        yBgMax = yBgMin;
    end
    zBgMin = zCenterLoc - zRange;
    zBgMax = zCenterLoc + zRange;
    if zBgMin < 1;
        zBgMin = 1;
    end
    if zBgMax > zDim
        zBgMax = zDim;
    end
    if zBgMax < zBgMin
        zBgMax = zBgMin;
    end
    
    backgroundImg = zeros(length(sliceImg(:,1)),length(sliceImg(1,:)));
    countBg = 0;
    for xBg = xBgMin:1:xBgMax
        for yBg = yBgMin:1:yBgMax
            for zBg = zBgMin:1:zBgMax
                countBg = countBg + 1;
                signalBg(countBg) = volData(yBg,xBg,zBg);
                if (xBg == xBgMin)||(xBg == xBgMax)||(yBg == yBgMin)||(yBg == yBgMax)
                    backgroundImg(yBg,xBg) = 1; % store the region used to calculate the signal
                end
            end
        end
    end
    
    % standard deviation of the background signal
    bgStandardDev = std(signalBg);
end

% calculate SNR
SNR = avgSig/bgStandardDev;

if plotSNRImg == 1;
    % plot image of center slice
    figure
    imshow(sliceImg,[])
    title('SNR Signal and BG Regions','FontSize',14)
    
    % create markers for bounds of the regions
    % marker plotting
    % cat(3,r,g,b)
    orange = cat(3, 255/255.*ones(size(backgroundImg)),165/255.*ones(size(backgroundImg)), zeros(size(backgroundImg)));
    hold on
    hOrange = imshow(orange);
    hold off
    set(hOrange, 'AlphaData', backgroundImg) % make color sheet only show markers
    
    % marker plotting
    % cat(3,r,g,b)
    white = cat(3, 102/255.*ones(size(signalImg)),255/255.*ones(size(signalImg)), 0/255.*ones(size(signalImg)));
    hold on
    hWhite = imshow(white);
    hold off
    set(hWhite, 'AlphaData', signalImg) % make color sheet only show markers
    
end

end