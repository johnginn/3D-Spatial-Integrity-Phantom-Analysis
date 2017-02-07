% Function to extract out the cardinal planes for conducting Procrustes
% method
%
% Input:
% centerRow The index of the center row in the phantom
% centerCol The index of the center col in the phantom
% centerSlice The index of the center slice in the phantom
% xSphere The x-location of the spheres found by the software
% ySphere The y-location of the spheres found by the software
% zSphere The z-location of the spheres found by the software
% xGndTruth The calculated ground truth x-location
% yGndTruth The calculated ground truth y-location
% zGndTruth The calculated ground truth z-location
%
% Output:
% xSpherePro The found x-positions of the spheres on the cardinal planes
% ySpherePro The found y-positions of the spheres on the cardinal planes
% zSpherePro The found z-positions of the spheres on the cardinal planes
% xGndPro The ground truth x-positions of the spheres on the cardinal planes
% yGndPro The ground truth y-positions of the spheres on the cardinal planes
% zGndPro The ground truth z-positions of the spheres on the cardinal planes
%
% John Ginn
% Created: 1/24/17
% Modified: 1/24/17

function [xSpherePro,ySpherePro,zSpherePro,xGndPro,yGndPro,zGndPro] =...
    selectCardinalPlanes(centerRow,centerCol,centerSlice,...
    xSphere,ySphere,zSphere,xGndTruth,yGndTruth,zGndTruth)

countCardPlanes = 0;
for step = 1:length(xSphere)
    % if this sphere is located on one of the cardinal planes store this
    % location.
    if (xGndTruth(step) == centerCol)||...
            (yGndTruth(step) == centerRow)||...
            (zGndTruth(step) == centerSlice)
        countCardPlanes = countCardPlanes + 1;
        xSpherePro(countCardPlanes) = xSphere(step);
        ySpherePro(countCardPlanes) = ySphere(step);
        zSpherePro(countCardPlanes) = zSphere(step);
        xGndPro(countCardPlanes) = xGndTruth(step);
        yGndPro(countCardPlanes) = yGndTruth(step);
        zGndPro(countCardPlanes) = zGndTruth(step);
    end
end

end