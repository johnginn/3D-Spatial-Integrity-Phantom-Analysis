% Function for calculating the new ground truth locations based off custom
% spacing of the spheres.
%
% Input:
% voxelHeight (mm/voxel) The height of the pixel
% voxelWidth (mm/voxel) The width of the pixel
% voxelLength (mm/voxel) The length of the pixel
% centerCol The index of the center column
% centerRow The index of the center row
% centerSlice The index of the center slice
% nRows The number of rows
% nSlices The number of slices
% spheresInRow Number of spheres in each row from bottom to top, for SI phantom [29 29 29 29 29 29 29 29 27 27 27 27 25 25 23 23 21 19 17 13 7]
% centerRowNum The number of the center row (not index, ie. 7 rows --> 4)
% xSpacing (mm) The x-spacing distance between the spheres
% ySpacing (mm) The y-spacing distance between the spheres
% zSpacing (mm) The z-spacing distance between the spheres
% xGroundTruthRemoved The ground truth x-locations after removing the spheres in a previous analysis
% yGroundTruthRemoved The ground truth y-locations after removing the spheres in a previous analysis
% zGroundTruthRemoved The ground truth z-locations after removing the spheres in a previous analysis
%
% Output:
% xSphereGroundTruthFinal The new ground truth x-locations
% ySphereGroundTruthFinal The new ground truth y-locations
% zSphereGroundTruthFinal The new ground truth z-locations
%
% John Ginn
% Created: 12/7/16
% Modified: 12/7/16

function [xSphereGroundTruthFinal,ySphereGroundTruthFinal,zSphereGroundTruthFinal] = ...
    calcNewGroundTruth(voxelHeight,voxelWidth, voxelLength,...
    centerCol,centerRow,centerSlice,nRows,nSlices,spheresInRow,centerRowNum,...
    xSpacing,ySpacing,zSpacing,xGroundTruthRemoved,yGroundTruthRemoved,zGroundTruthRemoved)

% phantom information
distBtwnSpheres = 16; % (mm) The distance between the spheres in each direction
xIndPerSphNew = xSpacing/voxelWidth; % number of pixels in x-direction between the spheres
yIndPerSphNew = ySpacing/voxelHeight; % number of pixels in y-direction between the spheres
zIndPerSphNew = zSpacing/voxelLength; % number of pixels in z-direction between the spheres
xIndPerSphOld = distBtwnSpheres/voxelWidth; % number of pixels in x-direction between the spheres
yIndPerSphOld = distBtwnSpheres/voxelHeight; % number of pixels in y-direction between the spheres
zIndPerSphOld = distBtwnSpheres/voxelLength; % number of pixels in z-direction between the spheres
searchDistWeightSmall = 3; % (mm) the weighed sum calculation distance for spheres near the wall of the phantom
for step = 1:nSlices
        if step < 5;
            % move downward inverted because of matlab indexing in plot
            sliceCalcNew(step) = centerSlice - round(zIndPerSphNew*step);
            sliceGndNew(step) = centerSlice - zIndPerSphNew*step;
            sliceCalcOld(step) = centerSlice - round(zIndPerSphOld*step);
            sliceGndOld(step) = centerSlice - zIndPerSphOld*step;
        elseif step == 5; % the center slice of spheres
            sliceCalcNew(step) = centerSlice;
            sliceGndNew(step) = centerSlice;
            sliceCalcOld(step) = centerSlice;
            sliceGndOld(step) = centerSlice;
        else
            % move upward, inverted because of indexing in plot
            % the count should start at 1 and move up (thus subtract 5)
            slicePos = (step - 5);
            sliceCalcNew(step) = centerSlice + round(zIndPerSphNew*slicePos);
            sliceGndNew(step) = centerSlice + zIndPerSphNew*slicePos;
            sliceCalcOld(step) = centerSlice + round(zIndPerSphOld*step);
            sliceGndOld(step) = centerSlice + zIndPerSphOld*slicePos;
        end
end
slicesGndNew = sort(sliceGndNew);
slicesGndOld = sort(sliceGndOld);

% calculate the locations of the center sphere and six surrounding spheres
nTotData = nSlices.*sum(spheresInRow);

% initialize data for efficiency
sphereGoundTruthNew = cell(1,nTotData);
sphereGoundTruthOld = cell(1,nTotData);
% array of the weighted sum distances used to calculate the locations
% cycle through slices
saveCount = 0; % used to save the data
for sliceStep = 1:nSlices
    % the current slice
    gndTruthSliceNew = slicesGndNew(sliceStep);
    gndTruthSliceOld = slicesGndOld(sliceStep);
    % step through the rows
    for rowStep = 1:nRows
            % move upward, inverted because of indexing in plot
            % the count should start at 1 and move up (thus subtract 11)
            rowPos = (rowStep - centerRowNum);
            gndTruthRowNew = centerRow - yIndPerSphNew*rowPos;
            gndTruthRowOld = centerRow - yIndPerSphOld*rowPos;
        % spheres remaining to step through in the row
        spheresRemaining = spheresInRow(rowStep); 
        stepRight = 0; % reset the step right counter
        % step through the columns
        for colStep = 1:spheresInRow(rowStep)
            % find the center sphere
            if spheresRemaining == (spheresInRow(rowStep))
                % do nothing since aready at the center location
                 gndTruthColNew = centerCol;
                 gndTruthColOld = centerCol;
            % move through spheres to the left 
            elseif (spheresRemaining > (spheresInRow(rowStep) - 1)/2 &&...
                spheresRemaining ~= (spheresInRow(rowStep)))
                stepAmount = colStep - 1; % remove one since the first step 
                % was used to calculate the center location
                % the current column location
                gndTruthColNew = centerCol - xIndPerSphNew*stepAmount;
                gndTruthColOld = centerCol - xIndPerSphOld*stepAmount;
            else
                % the current column location
                stepRight = stepRight + 1;
                gndTruthColNew = centerCol + xIndPerSphNew*stepRight;
                gndTruthColOld = centerCol + xIndPerSphOld*stepRight;
            end
            
            % ground truth locations (note rotation for plotting) [y x z]
            saveCount = saveCount + 1; % used to save the data
            sphereGoundTruthNew{saveCount} = [gndTruthRowNew gndTruthColNew gndTruthSliceNew];
            sphereGoundTruthOld{saveCount} = [gndTruthRowOld gndTruthColOld gndTruthSliceOld];
            % one more sphere has been found
            spheresRemaining = spheresRemaining - 1;
        end
    end
end

% Remove any locations that were removed in the analysis already
countSameSpheres = 0;
for step = 1:saveCount
    newLoc = sphereGoundTruthNew{step};
    CheckThisLoc = sphereGoundTruthOld{step};
    xCheckThisLoc(step) = CheckThisLoc(2);
    yCheckThisLoc(step) = CheckThisLoc(1);
    zCheckThisLoc(step) = CheckThisLoc(3);
    locFound = 0;
    countSteps = 1;
    % step through the removed spheres to determine new locations
    while (locFound == 0)&&(countSteps <= length(xGroundTruthRemoved))
       if (xCheckThisLoc(step) == xGroundTruthRemoved(countSteps))&&...
               (yCheckThisLoc(step) == yGroundTruthRemoved(countSteps))&&...
               (zCheckThisLoc(step) == zGroundTruthRemoved(countSteps))
           locFound = 1;
           countSameSpheres = countSameSpheres + 1;
           xSphereGroundTruthFinal(countSameSpheres) = newLoc(2);
           ySphereGroundTruthFinal(countSameSpheres) = newLoc(1);
           zSphereGroundTruthFinal(countSameSpheres) = newLoc(3);
       end
        countSteps = countSteps + 1; 
    end
end
