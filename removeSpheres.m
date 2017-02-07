% Function to remove spheres that are outliers found in regions that are
% not spheres. This function uses the deviation and distance from isocenter
% to reject certain spheres from the analysis region.
%
% Input:
% finalSphereDataXCorr The sphere x-locations after corrections
% finalSphereDataYCorr The sphere y-locations after corrections
% finalSphereDataYCorr The sphere z-locations after corrections
% xGndTruth The ground-truth x-locations
% yGndTruth The ground-truth y-locations
% zGndTruth The ground-truth z-locations
% radius The distance of the spheres from isocenter
% diffCorr The deviation from the ground truth for the spheres
% radiusBounds [min max] Bounds of sphere distances from isocenter you want to remove
% rangeDiffBounds [min max] Bounds of deviation for the spheres you want to remove
% xSphere The sphere x-locations before corrections
% ySphere The sphere y-locations before corrections
% zSphere The sphere z-locations before corrections
%
% Output:
% xDataCorrRemove Sphere x-locations after correction, after certain spheres have been removed
% yDataCorrRemove Sphere y-locations after correction, after certain spheres have been removed
% zDataCorrRemove Sphere z-locations after correction, after certain spheres have been removed
% xGndTruthRemove Ground truth x-locations after certain spheres have been removed
% yGndTruthRemove Ground truth y-locations after certain spheres have been removed
% zGndTruthRemove Ground truth z-locations after certain spheres have been removed
% radiusRemove Radius of spheres to isocenter after certain spheres have been removed
% diffRemove Deviation of spheres from ground truth after certain spheres have been removed
% xSphereRemove Sphere x-locations before correction, after certain spheres have been removed
% ySphereRemove Sphere y-locations before correction, after certain spheres have been removed
% zSphereRemove Sphere z-locations before correction, after certain spheres have been removed
%
% John Ginn
% Created: 8/16/16
% Modified: 11/8/16
function [xDataCorrRemove,yDataCorrRemove,zDataCorrRemove,...
    xGndTruthRemove,yGndTruthRemove,zGndTruthRemove,radiusRemove,diffRemove,...
    xSphereRemove,ySphereRemove,zSphereRemove] = ...
    removeSpheres(finalSphereDataXCorr,finalSphereDataYCorr,finalSphereDataZCorr,xGndTruth,yGndTruth,zGndTruth,...
    radius,diffCorr,radiusBounds,rangeDiffBounds,xSphere,ySphere,zSphere)

count = 0;
countSkip = 0;
for step = 1:length(finalSphereDataXCorr)
    if (radius(step)>=radiusBounds(1))&&(radius(step)<=radiusBounds(2))&&...
            (diffCorr(step)>=rangeDiffBounds(1))&&(diffCorr(step)<=rangeDiffBounds(2))
        % within the bounds that you want to skip
        countSkip = countSkip + 1;
    else
        % not in the region you want to skp
        count = count + 1;
        xDataCorrRemove(count) = finalSphereDataXCorr(step);
        yDataCorrRemove(count) = finalSphereDataYCorr(step);
        zDataCorrRemove(count) = finalSphereDataZCorr(step);
        xGndTruthRemove(count) = xGndTruth(step);
        yGndTruthRemove(count) = yGndTruth(step);
        zGndTruthRemove(count) = zGndTruth(step);
        xSphereRemove(count) = xSphere(step);
        ySphereRemove(count) = ySphere(step);
        zSphereRemove(count) = zSphere(step);
        diffRemove(count) = diffCorr(step);
        radiusRemove(count) = radius(step);
    end
end


end