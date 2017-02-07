% Function to plot the deviations for the three cardinal axes for
% debugging. Note to user: the length of both the ground truth data and the
% sphere locations must be the same and they must correspond to one another
% (each datapoint in coordSphere must correspond to each datapoint in sphereGroundTruth)
%
% Input:
% finalSphereDataX The x-location of the found spheres after the entire analysis
% finalSphereDataY The y-location of the found spheres after the entire analysis
% finalSphereDataZ The z-location of the found spheres after the entire analysis
% xGndTruthFinal The x-location of the ground truth positions after the entire analysis
% yGndTruthFinal The y-location of the ground truth positions after the entire analysis
% zGndTruthFinal The z-location of the ground truth positions after the entire analysis
% centerCol The centermost ground truth sphere x-location
% centerSlice The centermost ground truth sphere z-location
% centerRow The centermost ground truth sphere y-location
% voxelHeight (mm) y-dimension (anterior-posterior)
% voxelWidth (mm) x-dimension (medial-lateral)
% voxelLength (mm) z-dimension (inferior-superior)
% plotPlaneFit Whether or not to plot the fit
% correctionCount The number of the correction
%
% Output:
%
% John Ginn
% Created: 10/26/16
% Modified: 12/13/16

function [] = cardAxisPlot(finalSphereDataX,finalSphereDataY,finalSphereDataZ,...
    xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,centerCol,centerSlice,centerRow,...
    voxelHeight,voxelWidth,voxelLength,plotPlaneFit,correctionCount)

    % function to calculate the deviation between the sphere and ground
    % truth locations
    %
    % Input:
    % sphX The sphere x-position
    % sphY The sphere y-position
    % sphZ The sphere z-position
    % gndX The ground-truth x-position
    % gndY The ground-truth y-position
    % gndZ The ground-truth z-position
    % voxHeight (mm) y-dimension (anterior-posterior)
    % voxWidth (mm) x-dimension (medial-lateral)
    % voxLeng (mm) z-dimension (inferior-superior)
    % centCol The centermost ground truth sphere x-location
    % centSlice The centermost ground truth sphere z-location
    % centRow The centermost ground truth sphere y-location
    %
    % Output:
    % dev2D (mm) 2D deviation in the axial plane
    % dev3D (mm) 3D deviation
    % devX (mm) 1D deviation along x
    % devY (mm) 1D deviation along y
    % devZ (mm) 1D deviation along z
    % radius (mm) The disance from isocenter
    %
    function [dev2D, dev3D, devX, devY, devZ, radius] = calcDev(sphX,sphY,sphZ,gndX,gndY,gndZ,...
            voxHeight,voxWidth,voxLeng,centCol,centSlice,centRow)
        % initialize array
        dev2D = zeros(1,length(sphX));
        dev3D = zeros(1,length(sphX));
        devX = zeros(1,length(sphX));
        devY = zeros(1,length(sphX));
        devZ = zeros(1,length(sphX));
        radius = zeros(1,length(sphX));
        for stepDev = 1:length(sphX)
            % 2D deviation in axial plane
            dev2D(stepDev) = sqrt((voxWidth*(sphX(stepDev) - gndX(stepDev))).^2 + ...
                (voxHeight*(sphY(stepDev) - gndY(stepDev))).^2);
            % 3D deviation 
            dev3D(stepDev) = sqrt((voxWidth*(sphX(stepDev) - gndX(stepDev))).^2 + ...
                (voxHeight*(sphY(stepDev) - gndY(stepDev))).^2 + ...
                (voxLeng*(sphZ(stepDev) - gndZ(stepDev))).^2);
            % distance from isocenter
            radius(stepDev) = sqrt((voxWidth*(sphX(stepDev) - centCol)).^2 + ...
                (voxHeight*(sphY(stepDev) - centRow)).^2 + ...
                (voxLeng*(sphZ(stepDev) - centSlice)).^2);
            
            devX(stepDev) = voxWidth.*(sphX(stepDev) - gndX(stepDev));
            devY(stepDev) = voxHeight.*(sphY(stepDev) - gndY(stepDev));
            devZ(stepDev) = voxLeng.*(sphZ(stepDev) - gndZ(stepDev));
        end
    end

    % function to plot the deviations for each slice
    function [] = plotDev(sphX,sphY,sphZ,dev3D,devX,devY,devZ,radius,sliceStr)

        % plot the individual deviations along each direction
        figure;
        subplot(2,2,1)
        plot(sphX,devX,'.')
        xlabel('x-position (ind)','FontSize',16)
        ylabel('x-deviation (mm)','FontSize',16)
        title(sliceStr,'FontSize',16)
        subplot(2,2,2)
        plot(sphY,devY,'.')
        xlabel('y-position (ind)','FontSize',16)
        ylabel('y-deviation (mm)','FontSize',16)
        title(sliceStr,'FontSize',16)
        subplot(2,2,3)
        plot(sphZ,devZ,'.')
        xlabel('z-position (ind)','FontSize',16)
        ylabel('z-deviation (mm)','FontSize',16)
        title(sliceStr,'FontSize',16)
        subplot(2,2,4)
        plot(radius,dev3D,'.')
        xlabel('radius (mm)','FontSize',16)
        ylabel('3D-deviation (mm)','FontSize',16)
        title(sliceStr,'FontSize',16)
        
        % plot the 3D deviation
        figure;
        subplot(2,2,1)
        plot(sphX,dev3D,'.')
        xlabel('x-position (ind)','FontSize',16)
        ylabel('3D-deviation (mm)','FontSize',16)
        title(sliceStr,'FontSize',16)
        subplot(2,2,2)
        plot(sphY,dev3D,'.')
        xlabel('y-position (ind)','FontSize',16)
        ylabel('3D-deviation (mm)','FontSize',16)
        title(sliceStr,'FontSize',16)
        subplot(2,2,3)
        plot(sphZ,dev3D,'.')
        xlabel('z-position (ind)','FontSize',16)
        ylabel('3D-deviation (mm)','FontSize',16)
        title(sliceStr,'FontSize',16)
        subplot(2,2,4)
        plot(radius,dev3D,'.')
        xlabel('radius (mm)','FontSize',16)
        ylabel('3D-deviation (mm)','FontSize',16)
        title(sliceStr,'FontSize',16)
    end

% extract the individual axes
[sliceShiftSagDebug, indAnalysisSagSphereX,indAnalysisSagSphereY,indAnalysisSagSphereZ,...
    indAnalysisSagGroundX,indAnalysisSagGroundY,indAnalysisSagGroundZ]  = ...
    extractSlice(finalSphereDataX,finalSphereDataY,finalSphereDataZ,...
    xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,centerCol,centerSlice,centerRow,plotPlaneFit,'sagittal');
[sliceShiftCorDebug, indAnalysisCorSphereX,indAnalysisCorSphereY,indAnalysisCorSphereZ,...
    indAnalysisCorGroundX,indAnalysisCorGroundY,indAnalysisCorGroundZ]  = ...
    extractSlice(finalSphereDataX,finalSphereDataY,finalSphereDataZ,...
    xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,centerCol,centerSlice,centerRow,plotPlaneFit,'coronal');
[sliceShiftTransDebug, indAnalysisTransSphereX,indAnalysisTransSphereY,indAnalysisTransSphereZ,...
    indAnalysisTransGroundX,indAnalysisTransGroundY,indAnalysisTransGroundZ] = ...
    extractSlice(finalSphereDataX,finalSphereDataY,finalSphereDataZ,...
    xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,centerCol,centerSlice,centerRow,plotPlaneFit,'transverse');


% calculate the deviations for each of the cardinal slices
[dev2DSag, dev3DSag,  devXSag, devYSag, devZSag,radiusSag] = ...
    calcDev(indAnalysisSagSphereX,indAnalysisSagSphereY,indAnalysisSagSphereZ,...
    indAnalysisSagGroundX,indAnalysisSagGroundY,indAnalysisSagGroundZ,...
    voxelHeight,voxelWidth,voxelLength,centerCol,centerSlice,centerRow);
[dev2DCor, dev3DCor, devXCor, devYCor, devZCor, radiusCor] = ...
    calcDev(indAnalysisCorSphereX,indAnalysisCorSphereY,indAnalysisCorSphereZ,...
    indAnalysisCorGroundX,indAnalysisCorGroundY,indAnalysisCorGroundZ,...
    voxelHeight,voxelWidth,voxelLength,centerCol,centerSlice,centerRow);
[dev2DTrans, dev3DTrans,  devXTrans, devYTrans, devZTrans,radiusTrans] = ...
    calcDev(indAnalysisTransSphereX,indAnalysisTransSphereY,indAnalysisTransSphereZ,...
    indAnalysisTransGroundX,indAnalysisTransGroundY,indAnalysisTransGroundZ,...
    voxelHeight,voxelWidth,voxelLength,centerCol,centerSlice,centerRow);
% entire volume
[dev2DAll, dev3DAll,  devXAll, devYAll, devZAll,radiusAll] = ...
    calcDev(finalSphereDataX,finalSphereDataY,finalSphereDataZ,...
    xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,...
    voxelHeight,voxelWidth,voxelLength,centerCol,centerSlice,centerRow);


%% plot the results
plotDev(indAnalysisSagSphereX,indAnalysisSagSphereY,indAnalysisSagSphereZ,...
    dev3DSag,devXSag,devYSag,devZSag,radiusSag,'sagittal')
plotDev(indAnalysisCorSphereX,indAnalysisCorSphereY,indAnalysisCorSphereZ,...
    dev3DCor,devXCor,devYCor,devZCor,radiusCor,'coronal')
plotDev(indAnalysisTransSphereX,indAnalysisTransSphereY,indAnalysisTransSphereZ,...
    dev3DTrans,devXTrans,devYTrans,devZTrans,radiusTrans,'transverse')
plotString = strcat(['Rot Corrections = ',num2str(correctionCount)]);

plotDev(finalSphereDataX,finalSphereDataY,finalSphereDataZ,...
    dev3DAll,devXAll,devYAll,devZAll,radiusAll,plotString)


end