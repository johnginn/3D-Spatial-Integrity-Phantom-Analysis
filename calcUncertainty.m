% This function calculates the average deviation from expected distance between the sphere aligned
% to isocenter, and the 6 closest surrounding spheres. This metric provides
% information about the combined uncertainty of finding the spheres and
% manufacturing. (assuming homogeneous B0 field) Additionally, the average
% deviation from the ground truth locations is calculated for the 7
% centermost spheres (applying rotation and translation corrections
% separately for this analysis with procrustes method)
%
% Input:
% centerRow The index of the center row of spheres
% centerCol The index of the center column of spheres
% centerSlice The index of the center slice of spheres
% voxelHeight (mm) The height of each voxel
% voxelWidth (mm) The width of each voxel
% voxelLength (mm) The length of each voxel
% xSphereData All sphere x-locations of the sphere after corrections
% ySphereData All sphere y-locations of the sphere after corrections
% zSphereData All sphere z-locations of the sphere after corrections
% xGndData All ground-truth x-locations
% yGndData All ground-truth y-locations
% zGndData All ground-truth z-locations
% phantomSim y=1, n=0 Whether or not this data comes from a simulation
% xSpacing (mm) The x-spacing distance between the spheres
% ySpacing (mm) The y-spacing distance between the spheres
% zSpacing (mm) The z-spacing distance between the spheres
%
% Output:
% devFromExpected (mm) Array of the deviation of the sphere location from
% what is expected given the phantom specifications
%
% avgDevFromExpected (mm) The average deviation of the sphere location from
% what is expected given the phantom specifications
%
% numSpheresCalDevExp The number of spheres found and used to calculate the average deviation
% from what is expected given the phantom specifications
%
% calcSurLoc The calculated location of the surrounding spheres
%
% foundSurSphere The location of the surrounding spheres found by the
% software
%
% rmsDevFromExp The root-mean-square of the deviation of the measured distance 
% to the expected distance between the sphere aligned to isocenter, and the 
% surrounding spheres
%
% rmsDevFromGnd The root-mean-square of the deviation between the sphere
% location and ground truth location for the sphere aligned to isocenter
% and the six surrounding spheres
%
% John Ginn
% Created: 8/1/16
% Modified: 12/13/16

function [devFromExpected,avgDevFromExpected,numSpheresCalDevExp,...
    calcSurLoc,foundSurSphere,rmsDevFromExp,rmsDevFromGnd] = calcUncertainty(centerRow,centerCol,...
    centerSlice,voxelHeight,voxelWidth,voxelLength,xSphereData,ySphereData,zSphereData,...
    xGndData,yGndData,zGndData,phantomSim,xSpacing,ySpacing,zSpacing)

xIndPerSph = xSpacing/voxelWidth; % number of pixels in x-direction between the spheres
yIndPerSph = ySpacing/voxelHeight; % number of pixels in y-direction between the spheres
zIndPerSph = zSpacing/voxelLength; % number of pixels in z-direction between the spheres

% calculate the locations of the center sphere and six surrounding spheres
surSphere(1,:) = [centerRow centerCol centerSlice];
surSphere(2,:) = [(centerRow-yIndPerSph) centerCol centerSlice];
surSphere(3,:) = [(centerRow+yIndPerSph) centerCol centerSlice];
surSphere(4,:) = [centerRow (centerCol-xIndPerSph) centerSlice];
surSphere(5,:) = [centerRow (centerCol+xIndPerSph) centerSlice];
surSphere(6,:) = [centerRow centerCol (centerSlice-zIndPerSph)];
surSphere(7,:) = [centerRow centerCol (centerSlice+zIndPerSph)];

% the distance between the spheres in the simulated data may differ from
% 16 mm in the real phantom because the location of the spheres must be
% assigned an integer value when constructing the simulated phantom
if phantomSim == 1
    % calculate the locations of the center sphere and six surrounding spheres
    surSphere(1,:) = round([centerRow centerCol centerSlice]);
    surSphere(2,:) = round([(centerRow-yIndPerSph) centerCol centerSlice]);
    surSphere(3,:) = round([(centerRow+yIndPerSph) centerCol centerSlice]);
    surSphere(4,:) = round([centerRow (centerCol-xIndPerSph) centerSlice]);
    surSphere(5,:) = round([centerRow (centerCol+xIndPerSph) centerSlice]);
    surSphere(6,:) = round([centerRow centerCol (centerSlice-zIndPerSph)]);
    surSphere(7,:) = round([centerRow centerCol (centerSlice+zIndPerSph)]);
    xGndData = round(xGndData);
    yGndData = round(yGndData);
    zGndData = round(zGndData);
end

countFound = 0; % verify all the spheres were found
foundSurSphere = zeros(7,3);
foundSurGround = zeros(7,3);
% step through the spheres
for stepSphere = 1:length(xSphereData)
    % acquire the data
    gndTruthData = [yGndData(stepSphere), xGndData(stepSphere), zGndData(stepSphere)];
    sphereData = [ySphereData(stepSphere), xSphereData(stepSphere), zSphereData(stepSphere)];
    % step through the locations of the sphere you want to find
    for stepSurSpheres = 1:size(surSphere,1)
        % check if the current sphere is one of the spheres we are
        % looking for. Round to avoid precision error where == does not
        % detect the spheres
        if (round(gndTruthData(3)) == round(surSphere(stepSurSpheres,3)))&&...
                (round(gndTruthData(2)) == round(surSphere(stepSurSpheres,2)))&&...
                (round(gndTruthData(1)) == round(surSphere(stepSurSpheres,1)))
            % note: the index of the sphere locations in surSphere
            % corresponds to the spheres in this array
            countFound = countFound + 1;
            foundSurSphere(stepSurSpheres,:) = sphereData;
            foundSurGround(stepSurSpheres,:) = gndTruthData;
        end
    end
end

if (foundSurGround(1,1) == 0)&&(foundSurGround(1,2) == 0)&&(foundSurGround(1,3) == 0)
    error('center sphere not found!')
else
    centerSphereData = foundSurSphere(1,:);
end

% calculate the deviation of the spheres with respect to the center sphere
% location (minus one because the first sphere is the center location)
devFromExpected = zeros(size(surSphere,1) - 1,1);
deviationSphere = devFromExpected;
deviationGround = devFromExpected;
totDist = 0;
centerGndData = foundSurGround(1,:);
for stepFound = 1:(countFound - 1)
        % plus one location so that you skip the first location which would
        % just give zero
        % surrounding sphere deviation
        deviationSphere(stepFound) =...
        sqrt((voxelHeight*(foundSurSphere(stepFound+1,1) - centerSphereData(1)))^2 + ...
        (voxelWidth*(foundSurSphere(stepFound+1,2) - centerSphereData(2)))^2 +...
        (voxelLength*(foundSurSphere(stepFound+1,3) - centerSphereData(3)))^2);
        % ground-truth deviation
        deviationGround(stepFound) =...
    sqrt((voxelHeight*(foundSurGround(stepFound+1,1) - centerGndData(1)))^2 + ...
        (voxelWidth*(foundSurGround(stepFound+1,2) - centerGndData(2)))^2 +...
        (voxelLength*(foundSurGround(stepFound+1,3) - centerGndData(3)))^2);
    devFromExpected(stepFound) = abs(deviationSphere(stepFound) - deviationGround(stepFound));
    totDist = totDist + devFromExpected(stepFound);
end
% subtract one from the averaging because the center sphere does not
% deviate from itself...
numSpheresCalDevExp = (countFound - 1);
avgDevFromExpected = totDist/numSpheresCalDevExp;
calcSurLoc = surSphere; % the calculated positions of the surrounding spheres

% calculate the root mean square of the data
sumSquare = 0;
for step = 1:length(devFromExpected)
    sumSquare = sumSquare + (devFromExpected(step))^2;
end
rmsDevFromExp = sqrt(1/length(devFromExpected)*sumSquare);

% use procrustes method
proSphereBefore = foundSurSphere;
proGround = surSphere;
[d, proSphereAfter,transform] = procrustes(proGround,proSphereBefore,'scaling',false);

correctedSphere = proSphereAfter;


totDev = 0;
for stepFound = 1:(countFound - 1)
    % plus one location so that you skip the first location which would
    % just give zero
    % surrounding sphere deviation
    deviationSphereCorr(stepFound) =...
        sqrt((voxelHeight*(correctedSphere(stepFound+1,1) - centerSphereData(1)))^2 + ...
        (voxelWidth*(correctedSphere(stepFound+1,2) - centerSphereData(2)))^2 +...
        (voxelLength*(correctedSphere(stepFound+1,3) - centerSphereData(3)))^2);
    % ground-truth deviation
    devGroundToCorrSphere(stepFound) =...
        sqrt((voxelHeight*(foundSurGround(stepFound+1,1) - correctedSphere(stepFound+1,1)))^2 + ...
        (voxelWidth*(foundSurGround(stepFound+1,2) - correctedSphere(stepFound+1,2)))^2 +...
        (voxelLength*(foundSurGround(stepFound+1,3) - correctedSphere(stepFound+1,3)))^2);
    devFromExpectedCorr(stepFound) = abs(deviationSphere(stepFound) - deviationGround(stepFound));
    totDev = totDev + devGroundToCorrSphere(stepFound);
end
avgDev = totDev/length(devGroundToCorrSphere);

% calculate the root mean square of the data
sumSquareCorr = 0;
for step = 1:length(devGroundToCorrSphere)
    sumSquareCorr = sumSquareCorr + (devGroundToCorrSphere(step))^2;
end
rmsDevFromGnd = sqrt(1/length(devGroundToCorrSphere)*sumSquareCorr);
end