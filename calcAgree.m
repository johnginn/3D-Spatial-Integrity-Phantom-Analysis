% Function to rotate the data, shift the data and calculate the agreement
% between the ground truth and sphere locations.
% Note to user: the length of both the
% ground truth data and the sphere locations must be the same and they must
% correspond to one another (each datapoint in coordSphere must correspond
% to each datapoint in sphereGroundTruth)
% 
% Input:
% xSphere The x-locations of the spheres found by the program
% ySphere The y-locations of the spheres found by the program
% zSphere The z-locations of the spheres found by the program
% xGndTruthFinal The ground truth x-location of the spheres calculated by the program
% yGndTruthFinal The ground truth y-location of the spheres calculated by the program
% zGndTruthFinal The ground truth z-location of the spheres calculated by the program
% shiftForRot [column, row, slice] The shift to center the rotation (from plane fit typically)
% shiftCorrection [sagittal,coronal,transverse] The correction required by fitting each center plane
% RotationAngle The angle to rotate the data (from first fit)
% radiusSearch (mm) The radius from isocenter that defines the subvolume
% voxelHeight (mm) The height of the voxels
% voxelWidth (mm) The width of the voxels
% voxelLength (mm) The length of the voxels
% centerCol The index of the center column of the phantom
% centerRow The index of the center row of the phantom
% centerSlice The index of the center slice of the phantom
% countPass The number of spheres that met the tolerance threshold
% rotMatrixClock The clockwise rotation matrix
% rotMatrixCounter The counter-clockwise rotation matrix
%
%
% Output:
% xSphereRotClock All sphere x-locations rotated clockwise
% ySphereRotClock All sphere y-locations rotated clockwise
% zSphereRotClock All sphere z-locations rotated clockwise
% xSmallRotClock Subvolume within radiusSearch x-locations of the sphere rotated clockwise
% ySmallRotClock Subvolume within radiusSearch y-locations of the sphere rotated clockwise
% zSmallRotClock Subvolume within radiusSearch z-locations of the sphere rotated clockwise
% radiusSmallClock The radius from isocenter for clockwise rotated points within radiusSearch
% xSmallRotClockGnd Ground truth x-locations within radiusSearch corresponding to clockwise rotation
% ySmallRotClockGnd Ground truth y-locations within radiusSearch corresponding to clockwise rotation
% zSmallRotClockGnd Ground truth z-locations within radiusSearch corresponding to clockwise rotation
% avgDiffWithRotClock Avg. variance (sphere-gnd)^2 for all datapoints after clockwise rotation
% avgDiffSmallWithRotClock Avg. variance (sphere-gnd)^2 for points within radiusSearch after clockwise rotation
% currDevRotClock Variance (sphere-gnd)^2 for each point after clockwise rotation 
% currDevSmallRotClock Variance (sphere-gnd)^2 for points within radiusSearch after clockwise rotation
% xSphereRotCount All sphere x-locations rotated counter-clockwise
% ySphereRotCount All sphere y-locations rotated counter-clockwise
% zSphereRotCount All sphere z-locations rotated counter-clockwise
% radiusSmallCount The radius from isocenter for counter-clockwise rotated points within radiusSearch
% xSmallRotCount Subvolume within radiusSearch x-locations of the sphere rotated counter-clockwise
% ySmallRotCount Subvolume within radiusSearch y-locations of the sphere rotated counter-clockwise
% zSmallRotCount Subvolume within radiusSearch z-locations of the sphere rotated counter-clockwise
% xSmallRotCountGnd Ground truth x-locations within radiusSearch corresponding to counter-clockwise rotation
% ySmallRotCountGnd Ground truth y-locations within radiusSearch corresponding to counter-clockwise rotation
% zSmallRotCountGnd Ground truth z-locations within radiusSearch corresponding to counter-clockwise rotation
% avgDiffWithRotCount Avg. variance (sphere-gnd)^2 for all datapoints after clockwise rotation
% avgDiffSmallWithRotCount Avg. variance (sphere-gnd)^2 for points within radiusSearch after clockwise rotation 
% currDevRotCount Variance (sphere-gnd)^2 for each point after clockwise rotation
% currDevRotSmallCount Variance (sphere-gnd)^2 for points within radiusSearch after clockwise rotation
% radiusSmall The distance from isocenter for all points
% devSmallNoRot The deviation for all points with no rotation correction
% avgDiffNoRot Avg. variance (sphere-gnd)^2 for all datapoints without rotation
% avgDiffSmallNoRot Avg. variance (sphere-gnd)^2 for points within radiusSearch with no rotation
% 
% John Ginn
% Created: 7/6/16
% Modified: 8/16/16
function [xSphereRotClock,ySphereRotClock,zSphereRotClock,radiusSmallClock,...
    xSmallRotClock,ySmallRotClock,zSmallRotClock,...
    xSmallRotClockGnd,ySmallRotClockGnd,zSmallRotClockGnd,...
    avgDiffWithRotClock,avgDiffSmallWithRotClock,...
    currDevRotClock,currDevSmallRotClock,...
    xSphereRotCount,ySphereRotCount,zSphereRotCount,radiusSmallCount,...
    xSmallRotCount,ySmallRotCount,zSmallRotCount,...
    xSmallRotCountGnd,ySmallRotCountGnd,zSmallRotCountGnd,...
    avgDiffWithRotCount,avgDiffSmallWithRotCount,...
    currDevRotCount,currDevRotSmallCount,...
    radiusSmall,devSmallNoRot,avgDiffNoRot,avgDiffSmallNoRot] = ...
    calcAgree(xSphere,ySphere,zSphere,xGndTruthFinal,yGndTruthFinal,...
    zGndTruthFinal,shiftForRot,shiftCorrection,RotationAngle,radiusSearch,...
    voxelHeight,voxelWidth,voxelLength,centerCol,centerRow,centerSlice,countPass,...
    rotMatrixClock,rotMatrixCounter)




% Rotatate the spheres and determine how the deviation between the two sets 
% rotation is about y-axis in the image data 
% clockwise
rotationClock = rotMatrixClock(RotationAngle);
rotationCount = rotMatrixCounter(RotationAngle);
% shift the vector prior to rotation and after rotation to determine axis
% of rotation
shiftVectorX = @(shift,vector) [vector(1) + shift,vector(2),vector(3)];
shiftVectorY = @(shift,vector) [vector(1),vector(2) + shift,vector(3)];
shiftVectorZ = @(shift,vector) [vector(1),vector(2),vector(3) + shift];
% clockwise
SphereRotatedClock = zeros(length(xGndTruthFinal),3);
xSphereRotClock = zeros(length(xGndTruthFinal),1);
ySphereRotClock = zeros(length(xGndTruthFinal),1);
zSphereRotClock = zeros(length(xGndTruthFinal),1);
% counter-clockwise
SphereRotatedCount = zeros(length(xGndTruthFinal),3);
xSphereRotCount = zeros(length(xGndTruthFinal),1);
ySphereRotCount = zeros(length(xGndTruthFinal),1);
zSphereRotCount = zeros(length(xGndTruthFinal),1);


% rotate the data and calculate deviation
avgDiffNoRot = 0; % the average difference between the data that was not rotated
avgDiffWithRotClock = 0; % the average difference between the data that was  rotated clockwise
avgDiffWithRotCount = 0; % the average difference between the data that was  rotated counterclockwise
subSearchStep = 0;
subSearchStepClock = 0;
subSearchStepCount = 0;
avgDiffSmallNoRot = 0;
avgDiffSmallWithRotClock = 0;
avgDiffSmallWithRotCount = 0;
for step = 1:length(xGndTruthFinal)
    SphereRotatedClock(step,:) = [xSphere(step), ySphere(step), zSphere(step)];
    SphereRotatedCount(step,:) = [xSphere(step), ySphere(step), zSphere(step)];    
    % shift the vector so the ceter column is at the origin
    SphereRotatedClock(step,:)= shiftVectorX(-shiftForRot(1),SphereRotatedClock(step,:));
    SphereRotatedCount(step,:)= shiftVectorX(-shiftForRot(1),SphereRotatedCount(step,:));
    SphereRotatedClock(step,:)= shiftVectorY(-shiftForRot(2),SphereRotatedClock(step,:));
    SphereRotatedCount(step,:)= shiftVectorY(-shiftForRot(2),SphereRotatedCount(step,:));
    SphereRotatedClock(step,:)= shiftVectorZ(-shiftForRot(3),SphereRotatedClock(step,:));
    SphereRotatedCount(step,:)= shiftVectorZ(-shiftForRot(3),SphereRotatedCount(step,:));
    % for debugging
%     centeredData(step,:) = SphereRotatedCount(step,:);
%     if step == length(xGndTruthFinal)
%         figure
%         scatter3(centeredData(:,1),centeredData(:,2),centeredData(:,3))
%     end
    % now at origin, apply rotation
    SphereRotatedClock(step,:) = rotationClock*SphereRotatedClock(step,:)';
    SphereRotatedCount(step,:) = rotationCount*SphereRotatedCount(step,:)';
    % shift back
    SphereRotatedClock(step,:)= shiftVectorX(shiftForRot(1),SphereRotatedClock(step,:));
    SphereRotatedCount(step,:)= shiftVectorX(shiftForRot(1),SphereRotatedCount(step,:));
    SphereRotatedClock(step,:)= shiftVectorY(shiftForRot(2),SphereRotatedClock(step,:));
    SphereRotatedCount(step,:)= shiftVectorY(shiftForRot(2),SphereRotatedCount(step,:));
    SphereRotatedClock(step,:)= shiftVectorZ(shiftForRot(3),SphereRotatedClock(step,:));
    SphereRotatedCount(step,:)= shiftVectorZ(shiftForRot(3),SphereRotatedCount(step,:));
    % apply the extra shift determined by plane fitting
    SphereRotatedClock(step,:) = shiftVectorX(shiftCorrection(1),SphereRotatedClock(step,:));
    SphereRotatedCount(step,:) = shiftVectorX(shiftCorrection(1),SphereRotatedCount(step,:));
    SphereRotatedClock(step,:) = shiftVectorY(shiftCorrection(2),SphereRotatedClock(step,:));
    SphereRotatedCount(step,:) = shiftVectorY(shiftCorrection(2),SphereRotatedCount(step,:));
    SphereRotatedClock(step,:)= shiftVectorZ(shiftCorrection(3),SphereRotatedClock(step,:));
    SphereRotatedCount(step,:)= shiftVectorZ(shiftCorrection(3),SphereRotatedCount(step,:));
    % location of spheres 
    xSphereRotClock(step) = SphereRotatedClock(step,1); % sphere location x-component
    ySphereRotClock(step) = SphereRotatedClock(step,2); % sphere location y-component
    zSphereRotClock(step) = SphereRotatedClock(step,3); % sphere location z-component
    
    xSphereRotCount(step) = SphereRotatedCount(step,1); % sphere location x-component
    ySphereRotCount(step) = SphereRotatedCount(step,2); % sphere location y-component
    zSphereRotCount(step) = SphereRotatedCount(step,3); % sphere location z-component

    % Calcualte deviation between spheres and the ground truth (in mm)
    currDevNoRot(step) =  sqrt((voxelWidth*(xSphere(step) - xGndTruthFinal(step))).^2 + ...
        (voxelHeight*(ySphere(step) - yGndTruthFinal(step))).^2 + (voxelLength*(zSphere(step) - zGndTruthFinal(step))).^2);
    avgDiffNoRot = avgDiffNoRot + currDevNoRot(step); % avg deviation calc
    currDevRotClock(step) = sqrt((voxelWidth*(xSphereRotClock(step) - xGndTruthFinal(step))).^2 + ...
        (voxelHeight*(ySphereRotClock(step) - yGndTruthFinal(step))).^2 + (voxelLength*(zSphereRotClock(step) - zGndTruthFinal(step))).^2);
    avgDiffWithRotClock = avgDiffWithRotClock + currDevRotClock(step); % avg deviation calc
    currDevRotCount(step) = sqrt((voxelWidth*(xSphereRotCount(step) - xGndTruthFinal(step))).^2 + ...
        (voxelHeight*(ySphereRotCount(step) - yGndTruthFinal(step))).^2 + (voxelLength*(zSphereRotCount(step) - zGndTruthFinal(step))).^2);
    avgDiffWithRotCount = avgDiffWithRotCount + currDevRotCount(step); % avg deviation calc
    
    % calculate deviation between spheres within some radius of the center
    % of the phantom (no rotation)
    if (sqrt((voxelWidth*(xSphere(step) - centerCol)).^2 + (voxelHeight*(ySphere(step) - centerRow)).^2 +...
            (voxelLength*(zSphere(step) - centerSlice)).^2) < radiusSearch);
        subSearchStep = subSearchStep + 1;
        avgDiffSmallNoRot = avgDiffSmallNoRot + currDevNoRot(step);
        devSmallNoRot(subSearchStep) = currDevNoRot(step);
        
        % locations
        xSmall(subSearchStep) = xSphere(step);
        ySmall(subSearchStep) = ySphere(step);
        zSmall(subSearchStep) = zSphere(step);
        xSmallGnd(subSearchStep) = xGndTruthFinal(step);
        ySmallGnd(subSearchStep) = yGndTruthFinal(step);
        zSmallGnd(subSearchStep) = zGndTruthFinal(step);
        % radius from center
        radiusSmall(subSearchStep) = sqrt((voxelWidth*(xSphere(step) - centerCol)).^2 +...
            (voxelHeight*(ySphere(step) - centerRow)).^2 +...
            (voxelLength*(zSphere(step) - centerSlice)).^2);
    end
    % clockwise rotation
    if (sqrt((voxelWidth*(xSphereRotClock(step) - centerCol)).^2 + ...
            (voxelHeight*(ySphereRotClock(step) - centerRow)).^2 +...
            (voxelLength*(zSphereRotClock(step) - centerSlice)).^2) < radiusSearch);
        subSearchStepClock = subSearchStepClock + 1;
        % deviation
        avgDiffSmallWithRotClock = avgDiffSmallWithRotClock + currDevRotClock(step);
        currDevSmallRotClock(subSearchStepClock) = currDevRotClock(step);
        % locations
        xSmallRotClock(subSearchStepClock) = xSphereRotClock(step);
        ySmallRotClock(subSearchStepClock) = ySphereRotClock(step);
        zSmallRotClock(subSearchStepClock) = zSphereRotClock(step);
        xSmallRotClockGnd(subSearchStepClock) = xGndTruthFinal(step);
        ySmallRotClockGnd(subSearchStepClock) = yGndTruthFinal(step);
        zSmallRotClockGnd(subSearchStepClock) = zGndTruthFinal(step);
        % radius from center
        radiusSmallClock(subSearchStepClock) = sqrt((voxelWidth*(xSphereRotClock(step) - centerCol)).^2 + ...
            + (voxelHeight*(ySphereRotClock(step) - centerRow)).^2 +...
            (voxelLength*(zSphereRotClock(step) - centerSlice)).^2);
    end
    if (sqrt((voxelWidth*(xSphereRotCount(step) - centerCol)).^2 + ...
            (voxelHeight*(ySphereRotCount(step) - centerRow)).^2 +...
            (voxelLength*(zSphereRotCount(step) - centerSlice)).^2) < radiusSearch);
        % deviation
        subSearchStepCount = subSearchStepCount + 1;
        avgDiffSmallWithRotCount = avgDiffSmallWithRotCount + currDevRotCount(step);
        currDevRotSmallCount(subSearchStepCount) = currDevRotCount(step);
        % locations
        xSmallRotCount(subSearchStepCount) = xSphereRotCount(step);
        ySmallRotCount(subSearchStepCount) = ySphereRotCount(step);
        zSmallRotCount(subSearchStepCount) = zSphereRotCount(step);
        xSmallRotCountGnd(subSearchStepCount) = xGndTruthFinal(step);
        ySmallRotCountGnd(subSearchStepCount) = yGndTruthFinal(step);
        zSmallRotCountGnd(subSearchStepCount) = zGndTruthFinal(step);
        % radius from center
        radiusSmallCount(subSearchStepCount) = sqrt((voxelWidth*(xSphereRotCount(step) - centerCol)).^2 + ...
            (voxelHeight*(ySphereRotCount(step) - centerRow)).^2 +...
            (voxelLength*(zSphereRotCount(step) - centerSlice)).^2);
    end
end

% calculate the average deviation for all the datapoints
avgDiffNoRot = avgDiffNoRot/countPass;
avgDiffWithRotClock = avgDiffWithRotClock/countPass;
avgDiffWithRotCount = avgDiffWithRotCount/countPass;
% calculate the average deviation datapoints within specified distance from
% isocenter
avgDiffSmallNoRot = avgDiffSmallNoRot/subSearchStep;
avgDiffSmallWithRotClock = avgDiffSmallWithRotClock/subSearchStepClock;
avgDiffSmallWithRotCount = avgDiffSmallWithRotCount/subSearchStepCount;