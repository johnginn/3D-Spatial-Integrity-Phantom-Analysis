% Function to extract out different slices in the data for determining the
% x, y, and z shifts of the phantom at isocenter. This code used to
% calculate rotation corrections based off of a series of fit planes,
% however this method has been replaced with procrustes method.
% Note to user: the length of both the ground truth data and the sphere 
% locations must be the same and they must correspond to one another 
% (each datapoint in coordSphere must correspond to each datapoint in sphereGroundTruth)
% 
% Input:
% xSphere The x-locations of the spheres found by the program
% ySphere The y-locations of the spheres found by the program
% zSphere The z-locations of the spheres found by the program
% xGndTruthFinal The ground truth x-location of the spheres calculated by the program
% yGndTruthFinal The ground truth y-location of the spheres calculated by the program
% zGndTruthFinal The ground truth z-location of the spheres calculated by the program
% radiusSearch1 (mm) The radius of the first tolerance region
% voxelHeight (mm) The height of each voxel
% voxelWidth (mm) The width of each voxel
% voxelLength (mm) The length of each voxel
% centerCol The index of the center column of the phantom
% centerSlice The index of the center slice of the phantom
% centerRow The index of the center row of the phantom
% plotPlaneFit Whether or not you want to plot the fit of the plane
% sliceOpt The slice you want to extract, all slices pass through isocenter 'transverse' 'coronal' 'sagittal'
% countPass The number of spheres that met the tolerance threshold
% calcUncertainty Is this rotation correction being used for calculating
% uncertainty? y = 1, n = 0; If so, use the center sphere location for shift
%
% Output:
% RotationAngle The angle for the rotation correction
% finalSphereDataX The x-location of the spheres after rotation correction (all spheres)
% finalSphereDataY The y-location of the spheres after rotation correction (all spheres)
% finalSphereDataZ The z-location of the spheres after rotation correction (all spheres)
% finalSmallSphereDataX The x-location of the spheres after rotation correction (spheres in small region defined by radiusSearch1)
% finalSmallSphereDataY The y-location of the spheres after rotation correction (spheres in small region defined by radiusSearch1)
% finalSmallSphereDataZ The z-location of the spheres after rotation correction (spheres in small region defined by radiusSearch1)
% finalSmallGndX The x-ground-truth locations in the small volume defined by radiusSearch1
% finalSmallGndY The y-ground-truth locations in the small volume defined by radiusSearch1
% finalSmallGndZ The z-ground-truth locations in the small volume defined by radiusSearch1
% rotationDir The direction of rotation 'clockwise' or 'counter-clockwise'
% finalAvgDiff Avg. variance (sphere-gnd)^2 for all datapoints after rotation
% rotMatrixClock The clockwise rotation matrix
% rotMatrixCounter The counter-clockwise rotation matrix
% shiftForRot [column, row, slice] The shift to center the rotation (from plane fit)
%
% John Ginn
% Created: 8/5/16
% Modified: 9/2/16
function [RotationAngle, finalSphereDataX,finalSphereDataY,finalSphereDataZ,...
    finalSmallSphereDataX,finalSmallSphereDataY,finalSmallSphereDataZ,...
    finalSmallGndX,finalSmallGndY,finalSmallGndZ,rotationDir,finalAvgDiff,...
    rotMatrixClock,rotMatrixCounter,shiftForRot] = ...
    calcRotation(xSphere,ySphere,zSphere,xGndTruthFinal,yGndTruthFinal,...
    zGndTruthFinal,radiusSearch1,voxelHeight,voxelWidth,voxelLength,...
    centerCol,centerSlice,centerRow,plotPlaneFit,sliceOpt,countPass,calcUncertainty)

% use different center location for sphere calculation
if calcUncertainty == 1
    for step = 1:length(xGndTruthFinal)
        if (xGndTruthFinal(step) == centerCol)&&(yGndTruthFinal(step) == centerRow)&&...
                (zGndTruthFinal(step) == centerSlice)
            centerRowSphere = ySphere(step);
            centerColSphere = xSphere(step);
            centerSliceSphere = zSphere(step);
        end
    end
    
else
    centerRowSphere = centerRow;
    centerColSphere = centerCol;
    centerSliceSphere = centerSlice;
end


sliceCount = 0;

if strcmp(sliceOpt,'transverse') == 1
    for step = 1:length(zSphere)
        % extract out center transverse slice data
        if zGndTruthFinal(step) == centerSlice
            sliceCount = sliceCount + 1;
            sliceSphereX(sliceCount,1) = xSphere(step);
            sliceSphereY(sliceCount,1) = ySphere(step);
            sliceSphereZ(sliceCount,1) = zSphere(step);
            sliceGroundX(sliceCount,1) = xGndTruthFinal(step);
            sliceGroundY(sliceCount,1) = yGndTruthFinal(step);
            sliceGroundZ(sliceCount,1) = zGndTruthFinal(step);
        end
    end
    % determine type of fit
    % rotational
    planeModel = fittype( @(a,b,c,x,y) a*(x - b) + y.*0 + c, ...
        'independent', {'x', 'y'}, ...
        'dependent', 'z' );
    surfaceFitGnd = fit([sliceGroundX, sliceGroundY],sliceGroundZ,planeModel,...
        'StartPoint',[0,centerCol,centerSlice]);
    surfaceFitSphere = fit([sliceSphereX, sliceSphereY],sliceSphereZ,planeModel,...
        'StartPoint',[0,centerColSphere,centerSliceSphere]);
    % matrices for rotating the data
    rotMatrixClock = @(theta)[cosd(theta),  0,    -sind(theta);...
                              0,         1,          0;
                           sind(theta),  0,     cosd(theta)];
    % counterclockwise
    rotMatrixCounter =  @(theta)[cosd(theta),  0,    sind(theta);...
                                 0,         1,          0;
                              -sind(theta),  0,     cosd(theta)];
    coeffSph = coeffvalues(surfaceFitSphere);
    % the shift that must be applied to center the data before rotating the data
    shiftColFit = coeffSph(2); 
    shiftRowFit = centerRowSphere;
    shiftSliceFit = coeffSph(3); 
    shiftForRot = [shiftColFit, shiftRowFit, shiftSliceFit]; 
end
if strcmp(sliceOpt,'coronal') == 1
    if calcUncertainty == 1
        % just 7 spheres nearest to isocenter, use center coronal slice
        for step = 1:length(ySphere)
            % extract out center transverse slice data
            if yGndTruthFinal(step) == centerRow
                sliceCount = sliceCount + 1;
                sliceSphereX(sliceCount,1) = xSphere(step);
                sliceSphereY(sliceCount,1) = ySphere(step);
                sliceSphereZ(sliceCount,1) = zSphere(step);
                sliceGroundX(sliceCount,1) = xGndTruthFinal(step);
                sliceGroundY(sliceCount,1) = yGndTruthFinal(step);
                sliceGroundZ(sliceCount,1) = zGndTruthFinal(step);
            end
        end
        % rotational
        planeModel = fittype( @(a,b,c,z,x) a*(z - b) + x.*0 + c, ...
            'independent', {'z', 'x'}, ...
            'dependent', 'y' );
        surfaceFitGnd = fit([sliceGroundZ, sliceGroundX],sliceGroundY,planeModel,...
            'StartPoint',[0,centerSlice,centerRow]);
        surfaceFitSphere = fit([sliceSphereZ, sliceSphereX],sliceSphereY,planeModel,...
            'StartPoint',[0,centerSliceSphere,centerRowSphere]);
        % matrices for rotating the data
        rotMatrixClock = @(theta)[1,            0,              0;...
            0,        cosd(theta),    -sind(theta);
            0,        sind(theta),     cosd(theta)];
        % counterclockwise
        rotMatrixCounter =  @(theta)[1,            0,              0;...
            0,        cosd(theta),    sind(theta);
            0,        -sind(theta),     cosd(theta)];
        coeffSph = coeffvalues(surfaceFitSphere);
        % the shift that must be applied to center the data before rotating the data
        shiftColFit = centerColSphere;
        shiftRowFit = coeffSph(3);
        shiftSliceFit = coeffSph(2);
        shiftForRot = [shiftColFit, shiftRowFit, shiftSliceFit];
    else
        % entire phantom data, use center transverse slice since it has
        % more spheres than the center coronal slice
        for step = 1:length(zSphere)
            % extract out center transverse slice data
            if zGndTruthFinal(step) == centerSlice
                sliceCount = sliceCount + 1;
                sliceSphereX(sliceCount,1) = xSphere(step);
                sliceSphereY(sliceCount,1) = ySphere(step);
                sliceSphereZ(sliceCount,1) = zSphere(step);
                sliceGroundX(sliceCount,1) = xGndTruthFinal(step);
                sliceGroundY(sliceCount,1) = yGndTruthFinal(step);
                sliceGroundZ(sliceCount,1) = zGndTruthFinal(step);
            end
        end
    % rotational
    planeModel = fittype( @(a,b,c,y,x) a*(y - b) + x.*0 + c, ...
        'independent', {'y', 'x'}, ...
        'dependent', 'z' );
    surfaceFitGnd = fit([sliceGroundY, sliceGroundX],sliceGroundZ,planeModel,...
        'StartPoint',[0,centerRow,centerSlice]);
    surfaceFitSphere = fit([sliceSphereY, sliceSphereX],sliceSphereZ,planeModel,...
        'StartPoint',[0,centerRowSphere,centerSliceSphere]);
    % matrices for rotating the data
    rotMatrixClock = @(theta)[1,            0,              0;...
                              0,        cosd(theta),    -sind(theta);
                              0,        sind(theta),     cosd(theta)];
    % counterclockwise
    rotMatrixCounter =  @(theta)[1,            0,              0;...
                              0,        cosd(theta),    sind(theta);
                              0,        -sind(theta),     cosd(theta)];
    coeffSph = coeffvalues(surfaceFitSphere);
    % the shift that must be applied to center the data before rotating the data
    shiftColFit = centerColSphere; 
    shiftRowFit = coeffSph(2);
    shiftSliceFit = coeffSph(3); 
    shiftForRot = [shiftColFit, shiftRowFit, shiftSliceFit];      
    end
end
if strcmp(sliceOpt,'sagittal') == 1
    for step = 1:length(xSphere)
        % extract out center sagittal plane data
        if xGndTruthFinal(step) == centerCol
            sliceCount = sliceCount + 1;
            sliceSphereX(sliceCount,1) = xSphere(step);
            sliceSphereY(sliceCount,1) = ySphere(step);
            sliceSphereZ(sliceCount,1) = zSphere(step);
            sliceGroundX(sliceCount,1) = xGndTruthFinal(step);
            sliceGroundY(sliceCount,1) = yGndTruthFinal(step);
            sliceGroundZ(sliceCount,1) = zGndTruthFinal(step);
        end
    end
    % rotational
    planeModel = fittype( @(a,b,c,y,z) a*(y - b) + z.*0 + c, ...
        'independent', {'y', 'z'}, ...
        'dependent', 'x' );
    surfaceFitGnd = fit([sliceGroundY, sliceGroundZ],sliceGroundX,planeModel,...
        'StartPoint',[0,centerRow,centerCol]);
    surfaceFitSphere = fit([sliceSphereY, sliceSphereZ],sliceSphereX,planeModel,...
        'StartPoint',[0,centerRowSphere,centerColSphere]);
    % matrices for rotating the data
    rotMatrixClock = @(theta)[cosd(theta),  -sind(theta),    0;...
                              sind(theta),   cosd(theta),    0;
                              0,                  0,          1];
    % counterclockwise
    rotMatrixCounter = @(theta)[cosd(theta),    sind(theta),    0;...
                                -sind(theta),   cosd(theta),    0;
                                 0,                  0,          1];
    coeffSph = coeffvalues(surfaceFitSphere);
    % the shift that must be applied to center the data before rotating the data
    shiftColFit = coeffSph(3); 
    shiftRowFit = coeffSph(2);
    shiftSliceFit = centerSliceSphere; 
    shiftForRot = [shiftColFit, shiftRowFit, shiftSliceFit];       
end



% calculate the shift or rotation required
coeffGnd = coeffvalues(surfaceFitGnd);
aSph = coeffSph(1);
bSph = coeffSph(2);
cSph = coeffSph(3);
aGnd = coeffGnd(1);
bGnd = coeffGnd(2);
cGnd = coeffGnd(3);
% calculate rotation. Calculate normal vectors because MATLAB does not have a
% built in function to calculate the normal vectors from this type of fit,
% reconstruct datapoints from surface to find normal vectors using surfnorm
fitModel = @(a,b,c,x,y) a*(x - b) + y.*0 + c;
count = 1;
nPoints = 30;
xData = zeros(nPoints,nPoints);
yData = zeros(nPoints,nPoints);
gndFitData = zeros(nPoints,nPoints);
sphereFitData = zeros(nPoints,nPoints);
% calculate planes predicted by the fits
for x = 1:1:nPoints;
    for y = 1:1:nPoints;
        count = count + 1;
        % sphere fit data from fit
        xData(x,y) = x;
        yData(x,y) = y;
        sphereFitData(x,y) = fitModel(aSph,bSph,cSph,x,y);
        % ground truth fit data from fit
        gndFitData(x,y) = fitModel(aGnd,bGnd,cGnd,x,y);
    end
end
% the normal vectors
[NormSphereX,NormSphereY,NormSphereZ] = surfnorm(xData,yData,sphereFitData);
[NormGndX,NormGndY,NormGndZ] = surfnorm(xData,yData,gndFitData);
% all the normal vectors in the plane are the same, extract one normal
% vector
normVecSphere = [NormSphereX(1),NormSphereY(1),NormSphereZ(1)];
normVecGnd = [NormGndX(1),NormGndY(1),NormGndZ(1)];
sphereCrossGnd = cross(normVecSphere,normVecGnd);
% magnitude of the vectors
magCalc = @(vector) sqrt(vector(1)^2 + vector(2)^2 + vector(3)^2);
magSphereCrossGnd = magCalc(sphereCrossGnd);
magNormVecSphere = magCalc(normVecSphere);
magNormVecGnd = magCalc(normVecGnd);
% calculate the angle of rotation
RotationAngle = asind(magSphereCrossGnd/(magNormVecSphere*magNormVecGnd));


% from fitting
shiftColFit = bSph; % from fitting
shiftSliceFit = cSph; % from fitting
% select out a subvolume within the phantom
extraShift = [0 0 0]; % do not correct for isocenter shift until after rotation
% this will be determined by additional fitting later
%
% rotate data and calculate deviation
% NOTE: x,y,z must be specified appropriately in order to apply the
% rotation correction about the right axis
% if strcmp(sliceOpt,'transverse') == 1


    [xSphereRotClock,ySphereRotClock,zSphereRotClock,radiusSmallClock,...
        xSmallRotClock,ySmallRotClock,zSmallRotClock,...
        xSmallRotClockGnd,ySmallRotClockGnd,zSmallRotClockGnd,...
        avgDiffWithRotClock,avgDiffSmallWithRotClock,...
        currDevRotClock,currDevSmallRotClock,...
        xSphereRotCount,ySphereRotCount,zSphereRotCount,radiusSmallCount,...
        xSmallRotCount,ySmallRotCount,zSmallRotCount,...
        xSmallRotCountGnd,ySmallRotCountGnd,zSmallRotCountGnd,...
        avgDiffWithRotCount] = ...
        calcAgree(xSphere,ySphere,zSphere,xGndTruthFinal,yGndTruthFinal,...
        zGndTruthFinal,shiftForRot,extraShift,RotationAngle,radiusSearch1,...
        voxelHeight,voxelWidth,voxelLength,centerCol,centerRow,centerSlice,countPass,...
        rotMatrixClock,rotMatrixCounter);
% end

% clockwise rotation yeilded best results
if avgDiffWithRotClock < avgDiffWithRotCount
    finalSphereDataX = xSphereRotClock;
    finalSphereDataY = ySphereRotClock;
    finalSphereDataZ = zSphereRotClock;
    finalSmallSphereDataX = xSmallRotClock;
    finalSmallSphereDataY = ySmallRotClock;
    finalSmallSphereDataZ = zSmallRotClock;
    finalSmallGndX = xSmallRotClockGnd;
    finalSmallGndY = ySmallRotClockGnd;
    finalSmallGndZ = zSmallRotClockGnd;
    rotationDir = 'clockwise';
    finalAvgDiff = avgDiffWithRotClock;
else % counter-clockwise yeilds best results
    finalSphereDataX = xSphereRotCount;
    finalSphereDataY = ySphereRotCount;
    finalSphereDataZ = zSphereRotCount;
    finalSmallSphereDataX = xSmallRotCount;
    finalSmallSphereDataY = ySmallRotCount;
    finalSmallSphereDataZ = zSmallRotCount;
    finalSmallGndX = xSmallRotCountGnd;
    finalSmallGndY = ySmallRotCountGnd;
    finalSmallGndZ = zSmallRotCountGnd;
    rotationDir = 'counter-clockwise';
    finalAvgDiff = avgDiffWithRotCount;
end



% extract out the single slice that was rotated from the corrected data to
% check to make sure the correction is actually being implemented
% appropriately
sliceCount = 0;
if strcmp(sliceOpt,'transverse') == 1
    for step = 1:length(zGndTruthFinal)
        % extract out center transverse slice data
        if zGndTruthFinal(step) == centerSlice
            sliceCount = sliceCount + 1;
            xCorrPlot(sliceCount) = finalSphereDataX(step);
            yCorrPlot(sliceCount) = finalSphereDataY(step);
            zCorrPlot(sliceCount) = finalSphereDataZ(step);
        end
    end
end
if strcmp(sliceOpt,'coronal') == 1
    for step = 1:length(zGndTruthFinal)
        % extract out center transverse slice data
        if zGndTruthFinal(step) == centerSlice
            sliceCount = sliceCount + 1;
            xCorrPlot(sliceCount) = finalSphereDataX(step);
            yCorrPlot(sliceCount) = finalSphereDataY(step);
            zCorrPlot(sliceCount) = finalSphereDataZ(step);

        end
    end
end
if strcmp(sliceOpt,'sagittal') == 1
    for step = 1:length(zGndTruthFinal)
        % extract out center saggital slice data
        if xGndTruthFinal(step) == centerCol
            sliceCount = sliceCount + 1;
            xCorrPlot(sliceCount) = finalSphereDataX(step);
            yCorrPlot(sliceCount) = finalSphereDataY(step);
            zCorrPlot(sliceCount) = finalSphereDataZ(step);
        end
    end
end



if plotPlaneFit == 1
    % plot the slice that is extracted
    figure;
    % rotate the plotting for each case so that the fit plot matches the
    % correct axes as in the 
    if strcmp(sliceOpt,'transverse') == 1
        scatter3(sliceSphereX,sliceSphereY,sliceSphereZ)
        hold on
        scatter3(sliceGroundX,sliceGroundY,sliceGroundZ,'r+')
        hold on
        scatter3(xCorrPlot,yCorrPlot,zCorrPlot,'*k')
        legend('Sphere Fit','Ground Truth','Corrected Data')
        xlabel('x-axis (index)','FontSize',20)
        % switch to agree with clinical defintion
        ylabel('z-axis (index)','FontSize',20)
        zlabel('y-axis (index)','FontSize',20)
    elseif strcmp(sliceOpt,'coronal') == 1
        if calcUncertainty == 1;
            scatter3(sliceSphereZ,sliceSphereX,sliceSphereY)
            hold on
            scatter3(sliceGroundZ,sliceGroundX,sliceGroundY,'r+')
            hold on
            scatter3(zCorrPlot,xCorrPlot,yCorrPlot,'*k')
            legend('Sphere Fit','Ground Truth','Corrected Data')
            xlabel('z-axis (index)','FontSize',20)
            % switch to agree with clinical defintion
            ylabel('x-axis (index)','FontSize',20)
            zlabel('y-axis (index)','FontSize',20)
        else
            scatter3(sliceSphereY,sliceSphereX,sliceSphereZ)
            hold on
            scatter3(sliceGroundY,sliceGroundX,sliceGroundZ,'r+')
            hold on
            scatter3(yCorrPlot,xCorrPlot,zCorrPlot,'*k')
            legend('Sphere Fit','Ground Truth','Corrected Data')
            xlabel('z-axis (index)','FontSize',20)
            % switch to agree with clinical defintion
            ylabel('x-axis (index)','FontSize',20)
            zlabel('y-axis (index)','FontSize',20)
        end
    elseif strcmp(sliceOpt,'sagittal') == 1
        scatter3(sliceSphereY,sliceSphereZ,sliceSphereX)
        hold on
        scatter3(sliceGroundY,sliceGroundZ,sliceGroundX,'r+')
        hold on
        scatter3(yCorrPlot,zCorrPlot,xCorrPlot,'*k')
        legend('Sphere Fit','Ground Truth','Corrected Data')
        xlabel('z-axis (index)','FontSize',20)
        % switch to agree with clinical defintion
        ylabel('y-axis (index)','FontSize',20)
        zlabel('x-axis (index)','FontSize',20)
    end
    hold on
    plot(surfaceFitSphere) % the surface fit of the spheres that were found
    hold on
    plot(surfaceFitGnd) % the surface fit of the ground truth
    plotTitle = strcat([sliceOpt,' rotation: ',num2str(RotationAngle),' (degrees)']);
    title(plotTitle,'FontSize',20)
end


end