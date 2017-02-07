% function to simulate a phantom of spheres to test to see if the program is
% working correctly
%
% Input:
% voxelHeight (mm/voxel) The height of the pixel
% voxelWidth (mm/voxel) The width of the pixel
% voxelLength (mm/voxel) The length of the pixel
% sphereDiameter (mm) The diameter of the sphere
% nSpheres Number of spheres in each direction i.e. row, col, slices
% sphereBoundary Whether or not you want the sphere model to be surrounded
% by a cubic shell of 0's y = 1, n = 0
%
% Output: 
% phantomVol The phantom volume
% slices The slices in the phantom (may be rounded)
% sliceGnd The ground truth slices (not rounded)
% gndTruth The ground-truth locations of the spheres in the phantom
%
% John Ginn
% Created: 7/1/16
% Modified: 8/16/16

function [phantomVol, slices, sliceGnd,gndTruth] = ...
    makePhantom(voxelHeight, voxelWidth, voxelLength, sphereDiameter,nSpheres,sphereBoundary)

% an individual sphere
[sphereVol] = makeSphere(voxelHeight, voxelWidth, voxelLength, sphereDiameter,sphereBoundary,'mr');
xDimStep = length(sphereVol(1,:,1));
yDimStep = length(sphereVol(:,1,1));
zDimStep = length(sphereVol(1,1,:));
if mod(xDimStep,2) == 0; % x-dimension of sphere is even
    sphRangeX = xDimStep/2;
else % x-dimension of sphere is odd
    sphRangeX = (xDimStep-1)/2;
end
if mod(yDimStep,2) == 0; % y-dimension of sphere is even
    sphRangeY = yDimStep/2;
else % y-dimension of sphere is odd
    sphRangeY = (yDimStep-1)/2;
end
if mod(zDimStep,2) == 0; % z-dimension of sphere is even
    sphRangeZ = zDimStep/2;
else % z-dimension of sphere is odd
    sphRangeZ = (zDimStep-1)/2;
end
% spacing between the spheres (16 mm between spheres - 8 mm sphere
% diameter)
xSpacing = round(16/voxelWidth) - xDimStep; 
ySpacing = round(16/voxelHeight) - yDimStep;
zSpacing = round(16/voxelLength) - zDimStep;
xIndPerSph = 16/voxelWidth; % number of pixels in x-direction between the spheres
yIndPerSph = 16/voxelHeight; % number of pixels in y-direction between the spheres
zIndPerSph = 16/voxelLength; % number of pixels in z-direction between the spheres
xArray = 1:xDimStep;
yArray = 1:yDimStep;
zArray = 1:zDimStep;
emptyBounds = 3; % number of spacing to avoid going out of bounds

xDimension = 2*emptyBounds*(length(xArray)) + nSpheres*(xSpacing + xDimStep); 
yDimension = 2*emptyBounds*(length(yArray)) + nSpheres*(ySpacing + yDimStep); 
zDimension = 2*emptyBounds*(length(zArray)) + nSpheres*(zSpacing + zDimStep); 
phantomVol = zeros(yDimension,xDimension,zDimension);

% number of the center row;
centerCol = round(median(1:xDimension)); % location of center column
centerRow = round(median(1:yDimension)); % location of the center row
centerSlice = round(median(1:zDimension)); % location of the center slice
centerRowI = round(median(1:nSpheres)); % the center row number
centerSliceI = round(median(1:nSpheres)); % the center row number

count = 1;
sliceCount = 1;

% calculate the slice positions. The actual voume will be constructed in
% another set of loops to ensure that the order of the sphere locations
% calculated in this code correspond directly with those calculated in the
% PhantomScan.m script
for sliceStep = 1:nSpheres
    if sliceStep < centerSliceI;
        % move downward inverted because of matlab indexing in plot
        currentSlice = centerSlice + round(zIndPerSph*sliceStep);
        currentSliceGnd = centerSlice + zIndPerSph*sliceStep;
    elseif sliceStep == centerSliceI; % the center row
        currentSlice = centerSlice;
        currentSliceGnd = centerSlice;
    else
        % move upward, inverted because of indexing in plot
        % the count should start at 1 and move up (thus subtract 11)
        slicePos = (sliceStep - centerSliceI);
        currentSlice = centerSlice - round(zIndPerSph*slicePos);
        currentSliceGnd = centerSlice - zIndPerSph*slicePos;
    end
    % step through rows
    for rowStep = 1:nSpheres
        if rowStep < centerRowI;
            % move downward inverted because of matlab indexing in plot
            currentRow = centerRow + round(yIndPerSph*rowStep);
        elseif rowStep == centerRowI; % the center row
            currentRow = centerRow;
        else
            % move upward, inverted because of indexing in plot
            % the count should start at 1 and move up (thus subtract 11)
            rowPos = (rowStep - centerRowI);
            currentRow = centerRow - round(yIndPerSph*rowPos);
        end
        % spheres remaining to step through in the row
        spheresRemaining = nSpheres;
        stepRight = 1; % reset the step right counter
        % step through the columns
        currentCol = centerCol; % start at the center
        for colStep = 1:nSpheres
            % find the center sphere
            if spheresRemaining == (nSpheres)
                % do nothing since aready at the center location
                
                % move through spheres to the left
            elseif (spheresRemaining > (nSpheres - 1)/2 &&...
                    spheresRemaining ~= (nSpheres))
                stepAmount = colStep - 1; % remove one since the first step
                % was used to calculate the center location
                % the current column location
                currentCol = centerCol - round(xIndPerSph*stepAmount);
            else
                % the current column location
                currentCol = centerCol + round(xIndPerSph*stepRight);
                stepRight = stepRight + 1;
            end
            % one more sphere has been added
            spheresRemaining = spheresRemaining - 1;
            % store the slices
            if (rowStep == 1)&&(colStep == 1)
               slices(sliceCount) = currentSlice;
               sliceGnd(sliceCount) = currentSliceGnd;
               sliceCount = sliceCount + 1;
            end
        end
    end
end


% Now make the phantom data. Calculate the sphere locations in the same
% order as they will be calculated in the PhantomScan.m script
% step through slices
%
% sort the slices so they will be in the correct order
saveCount = 0;
slices = sort(slices);
for sliceStep = 1:nSpheres
    % the current slice
    currentSlice = slices(sliceStep);
    zStore = (currentSlice - sphRangeZ):1:(currentSlice + sphRangeZ);
    % step through the rows
    for rowStep = 1:nSpheres
        % move upward, inverted because of indexing in plot
        % the count should start at 1 and move up (thus subtract 11)
        rowPos = (rowStep - centerSliceI);
        currentRow = centerRow - round(yIndPerSph*rowPos);
        yStore = (currentRow - sphRangeY):1:(currentRow + sphRangeY);
        % spheres remaining to step through in the row
        spheresRemaining = nSpheres;
        stepRight = 0; % reset the step right counter
        % step through the columns
        currentCol = centerCol; % start at the center
        for colStep = 1:nSpheres
            % find the center sphere
            if spheresRemaining == (nSpheres)
                % do nothing since aready at the center location
                % move through spheres to the left
            elseif (spheresRemaining > (nSpheres - 1)/2 &&...
                    spheresRemaining ~= (nSpheres))
                stepAmount = colStep - 1; % remove one since the first step
                % was used to calculate the center location
                % the current column location
                currentCol = centerCol - round(xIndPerSph*stepAmount);
                
            else
                % the current column location
                stepRight = stepRight + 1;
                currentCol = centerCol + round(xIndPerSph*stepRight);
            end
            % ground truth locations (note rotation for plotting) [y x z]
            saveCount = saveCount + 1; % used to save the data
            gndTruth{saveCount} = [currentRow currentRow currentCol];
            xStore = (currentCol - sphRangeX):1:(currentCol + sphRangeX);
            phantomVol(yStore,xStore,zStore) = sphereVol;
            % one more sphere has been found
            spheresRemaining = spheresRemaining - 1;
        end
    end
end

end
