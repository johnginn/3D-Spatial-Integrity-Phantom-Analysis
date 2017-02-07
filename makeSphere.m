% This function makes the sphere template used to step through the images to
% calculate a correlation coefficient
%
% Input:
% voxelHeight (mm/voxel) The height of the pixel
% voxelWidth (mm/voxel) The width of the pixel
% voxelLength (mm/voxel) The length of the pixel
% sphereDiameter (mm) The diameter of the sphere
% sphereBoundary Whether or not you want the sphere model to be surrounded
% by a cubic shell of 0's y = 1, n = 0
% MRorCT ('mr' or 'ct') Determines whether the template contrast exhibits a CT or MR scan 
%
% Output:
% sphereVol (pixels) The volume of the sphere (1's and 0's) where the sphere exists
%
% John Ginn
% Created: 6/24/16
% Modified: 10/20/16

function [sphereVol] = makeSphere(voxelHeight, voxelWidth, voxelLength, sphereDiameter,sphereBoundary,MRorCT)
sphereRadius = sphereDiameter/2;

% find necessary range to sample
xRange = round(sphereRadius/voxelWidth);
yRange = round(sphereRadius/voxelHeight);
zRange = round(sphereRadius/voxelLength);

% used for storing the sphere data
countX = 1;
countY = 1;
countZ = 1;

% temporary sphere model until it is determine whether or not the sphere
% needs to be trimmed
sphereVolTemp = zeros(length(-yRange:1:yRange),length(-xRange:1:xRange),length(-zRange:1:zRange));
for X = -xRange:1:xRange;
    % reset the count in the matrix
    countY = 1;
    for Y = -yRange:1:yRange;
        % reset the count in the matrix
        countZ = 1;
        for Z = -zRange:1:zRange;
            % sphere equation --> x^2 + y^2 + z^2 = r^2
            if ((X*voxelWidth)^2 + (Y*voxelHeight)^2 + (Z*voxelLength)^2 <= sphereRadius^2)
                % if inside the sphere, or on boundary assign a 1
                sphereVolTemp(countY,countX,countZ) = 1;
            end
            countZ = countZ + 1;
        end
        countY = countY + 1;
    end
    countX = countX + 1;
end

% count the number of ones on the x boundary
countOnesX = 0;
countOnesY = 0;
countOnesZ = 0;
% dimensions of the sphere
yDim = length(sphereVolTemp(:,1,1));
xDim = length(sphereVolTemp(1,:,1));
zDim = length(sphereVolTemp(1,1,:));
for X = 1:xDim
    for Y = 1:yDim
        for Z = 1:zDim
            if (X == 1)||(X == xDim)
                if (sphereVolTemp(Y,X,Z) == 1)
                    % there is a 1 on one of the x-edges of the volume
                    countOnesX = countOnesX + 1;
                end
            end
            if (Y == 1)||(Y == yDim)
                if (sphereVolTemp(Y,X,Z) == 1)
                    % there is a 1 on one of the y-edges of the volume
                    countOnesY = countOnesY + 1;
                end
            end
            if (Z == 1)||(Z == zDim)
                if (sphereVolTemp(Y,X,Z) == 1)
                    % there is a 1 on one of the z-edges of the volume
                    countOnesZ = countOnesZ + 1;
                end
            end
        end
    end
end

% determine if there will be a boundary of zeros or not
if sphereBoundary == 1 % yes, boundary of zeros
    % need to add zeros to each dimension separately, because the voxel
    % size might not always be symmetric
    if (countOnesX > 0)
        % add a boundary
        xSize = xDim + 2;
        % the location in the new dataset where the old data will be stored
        xArray = (1:xDim) + 1;
    else
        % no boundary needed
        xSize = xDim;
        xArray = (1:xDim);
    end
    if (countOnesY > 0)
        % add a boundary
        ySize = yDim + 2;
        % the location in the new dataset where the old data will be stored
        yArray = (1:yDim) + 1;
    else
        % no boundary needed
        ySize = yDim;
        yArray = (1:yDim);
    end
    if (countOnesZ > 0)
        % add a boundary
        zSize = zDim + 2;
        % the location in the new dataset where the old data will be stored
        zArray = (1:zDim) + 1;
    else
        % no boundary needed
        zSize = zDim;
        zArray = (1:zDim);
    end
    % initialize the new data
    sphereVol = zeros(ySize,xSize,zSize);
    sphereVol(yArray,xArray,zArray) = sphereVolTemp;

else % no, no boundary of zeros
    % trim the edges of the data if sphereRadius/voxelWidth is not an integer
    % value. (avoids using 0's surrounding the sphere model in the calculation
    % of the correlation coefficient)
    if(countOnesX > 0)
        % don't trim x-dimension
        xData = 1:length(sphereVolTemp(1,:,1));
    else
        % trim the x-dimension
        xData = 2:(length(sphereVolTemp(1,:,1)) - 1);
    end
    if(countOnesY > 0)
        % don't trim y-dimension
        yData = 1:length(sphereVolTemp(:,1,1));
    else
        % trim the y-dimension
        yData = 2:(length(sphereVolTemp(:,1,1)) - 1);
    end
    if(countOnesZ > 0)
        % don't trim z-dimension
        zData = 1:length(sphereVolTemp(1,1,:));
    else
        % trim the z-dimension
        zData = 2:(length(sphereVolTemp(1,1,:)) - 1);
    end
    % store the final data
    xDataFinal = 1:length(xData); % x-dimension may be a different length, start at 1
    yDataFinal = 1:length(yData); % y-dimension may be a different length, start at 1
    zDataFinal = 1:length(zData); % z-dimension may be a different length, start at 1
    sphereVol(yDataFinal,xDataFinal,zDataFinal) = sphereVolTemp(yData,xData,zData);
end

% if you are using a CT scan, invert the contrast
if (strcmp(MRorCT,'ct') == 1)
   for xStep = 1:length(sphereVol(:,1,1)) 
       for yStep = 1:length(sphereVol(1,:,1))
           for zStep = 1:length(sphereVol(1,1,:))
               if sphereVol(xStep,yStep,zStep) > 0
                   sphereVol(xStep,yStep,zStep) = 0;
               else
                   sphereVol(xStep,yStep,zStep) = 1;
               end
           end
       end
   end
end

end