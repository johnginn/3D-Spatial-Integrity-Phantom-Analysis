% Function to extract out different slices in the data.
% Note to user: the length of both the ground truth data and the sphere 
% locations must be the same and they must correspond to one another 
% (each datapoint in coordSphere must correspond to each datapoint in sphereGroundTruth)
% Fitting code is from an older version of the analysis code
% 
% Input:
% sphereLocX The x-locations of the spheres found by the program
% sphereLocY The y-locations of the spheres found by the program
% sphereLocZ The z-locations of the spheres found by the program
% gndTruthLocX The ground truth x-location of the spheres calculated by the program
% gndTruthLocY The ground truth y-location of the spheres calculated by the program
% gndTruthLocZ The ground truth z-location of the spheres calculated by the program
% centerCol The index of the center column of the phantom
% centerSlice The index of the center slice of the phantom
% centerRow The index of the center row of the phantom
% plotPlaneFit Whether or not you want to plot the fit of the plane
% sliceOpt The slice you want to extract, all slices pass through isocenter
% 'transverse' 'coronal' 'sagittal'
%
% Output:
% sliceShift The shift from fitting between the ground truth and sphere location planes
% sliceSphereX The x-location of the spheres in the slice extracted
% sliceSphereY The y-location of the spheres in the slice extracted
% sliceSphereZ The z-location of the spheres in the slice extracted
% sliceGroundX The x-location of the ground truth in the slices
% sliceGroundY The y-location of the ground truth in the slices
% sliceGroundZ The z-location of the ground truth in the slices
%
% John Ginn
% Created: 7/6/16
% Modified: 10/27/16
function [sliceShift, sliceSphereX,sliceSphereY,sliceSphereZ,...
    sliceGroundX,sliceGroundY,sliceGroundZ] = ...
    extractSlice(sphereLocX,sphereLocY,sphereLocZ,gndTruthLocX,gndTruthLocY,...
    gndTruthLocZ,centerCol,centerSlice,centerRow,plotPlaneFit,sliceOpt)

sliceCount = 0;
if strcmp(sliceOpt,'transverse') == 1
    for step = 1:length(gndTruthLocZ)
        % extract out center transverse slice data
        if gndTruthLocZ(step) == centerSlice
            sliceCount = sliceCount + 1;
            sliceSphereX(sliceCount,1) = sphereLocX(step);
            sliceSphereY(sliceCount,1) = sphereLocY(step);
            sliceSphereZ(sliceCount,1) = sphereLocZ(step);
            sliceGroundX(sliceCount,1) = gndTruthLocX(step);
            sliceGroundY(sliceCount,1) = gndTruthLocY(step);
            sliceGroundZ(sliceCount,1) = gndTruthLocZ(step);
        end
    end
    planeModel = fittype( @(c,x,y) 0*(x) + y.*0 + c, ...
        'independent', {'x', 'y'}, ...
        'dependent', 'z' ); 
    surfaceFitGnd = fit([sliceGroundX, sliceGroundY],sliceGroundZ,planeModel,...
    'StartPoint',centerSlice);
    surfaceFitSphere = fit([sliceSphereX, sliceSphereY],sliceSphereZ,planeModel,...
    'StartPoint',centerSlice);
end
if strcmp(sliceOpt,'coronal') == 1
    for step = 1:length(gndTruthLocY)
        % extract out center sagittal plane data
        if gndTruthLocY(step) == centerRow
            sliceCount = sliceCount + 1;
            sliceSphereX(sliceCount,1) = sphereLocX(step);
            sliceSphereY(sliceCount,1) = sphereLocY(step);
            sliceSphereZ(sliceCount,1) = sphereLocZ(step);
            sliceGroundX(sliceCount,1) = gndTruthLocX(step);
            sliceGroundY(sliceCount,1) = gndTruthLocY(step);
            sliceGroundZ(sliceCount,1) = gndTruthLocZ(step);
        end
    end
    planeModel = fittype( @(c,x,z) 0*(x) + z.*0 + c, ...
    'independent', {'x', 'z'}, ...
    'dependent', 'y' ); 
    surfaceFitGnd = fit([sliceGroundX, sliceGroundZ],sliceGroundY,planeModel,...
    'StartPoint',centerRow);
    surfaceFitSphere = fit([sliceSphereX, sliceSphereZ],sliceSphereY,planeModel,...
    'StartPoint',centerRow);
end
if strcmp(sliceOpt,'sagittal') == 1
    for step = 1:length(gndTruthLocX)
        % extract out center sagittal plane data
        if gndTruthLocX(step) == centerCol
            sliceCount = sliceCount + 1;
            sliceSphereX(sliceCount,1) = sphereLocX(step);
            sliceSphereY(sliceCount,1) = sphereLocY(step);
            sliceSphereZ(sliceCount,1) = sphereLocZ(step);
            sliceGroundX(sliceCount,1) = gndTruthLocX(step);
            sliceGroundY(sliceCount,1) = gndTruthLocY(step);
            sliceGroundZ(sliceCount,1) = gndTruthLocZ(step);
        end
    end
    planeModel = fittype( @(c,y,z) 0*(y) + z.*0 + c, ...
    'independent', {'y', 'z'}, ...
    'dependent', 'x' ); 
    surfaceFitGnd = fit([sliceGroundY, sliceGroundZ],sliceGroundX,planeModel,...
    'StartPoint',centerCol);
    surfaceFitSphere = fit([sliceSphereY, sliceSphereZ],sliceSphereX,planeModel,...
    'StartPoint',centerCol);
end
% calculate the shift required
coeffSph = coeffvalues(surfaceFitSphere);
coeffGnd = coeffvalues(surfaceFitGnd);
sliceShift = coeffGnd - coeffSph;

if plotPlaneFit == 1
    % calculate the values of the locations after the shift
    correctedX = zeros(length(sliceSphereX),1);
    correctedY = zeros(length(sliceSphereX),1);
    correctedZ = zeros(length(sliceSphereX),1);
    for step = 1:length(sliceSphereX)
        if strcmp(sliceOpt,'transverse') == 1
            correctedX(step) = sliceSphereX(step);
            correctedY(step) = sliceSphereY(step);
            correctedZ(step) = sliceSphereZ(step) + sliceShift;
        elseif strcmp(sliceOpt,'coronal') == 1
            correctedX(step) = sliceSphereX(step);
            correctedY(step) = sliceSphereY(step) + sliceShift;
            correctedZ(step) = sliceSphereZ(step);
        elseif strcmp(sliceOpt,'sagittal') == 1
            correctedX(step) = sliceSphereX(step) + sliceShift;
            correctedY(step) = sliceSphereY(step);
            correctedZ(step) = sliceSphereZ(step);
        end
    end
    
    % plot the slice that is extracted
    figure;
    scatter3(sliceSphereX,sliceSphereY,sliceSphereZ)
    hold on
    scatter3(sliceGroundX,sliceGroundY,sliceGroundZ,'r+')
    hold on
    scatter3(correctedX,correctedY,correctedZ,'k*')
    plotTitle = strcat([sliceOpt,' center slice, shift: ',num2str(sliceShift),' (index)']);
    % xlim([min(sphereLocX),max(sphereLocX)])
    % ylim([min(sphereLocY),max(sphereLocY)])
    % zlim([min(sphereLocZ) max(sphereLocZ)]);
    legend('Sphere Fit','Ground Truth','Corrected Locations')
    xlabel('x-axis (index)','FontSize',20)
    % switch to agree with clinical defintion
    ylabel('z-axis (index)','FontSize',20)
    zlabel('y-axis (index)','FontSize',20)
    title(plotTitle,'FontSize',20)
end





end