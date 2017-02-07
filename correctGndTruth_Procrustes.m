% Function to rotate and translate the ground truth data to make it better
% align with the sphere locations in the image. Applies the reverse
% corrections optained from Procrustes method
%
% Input:
% proT The transform returned by procrustes metheod (see doc procrustes)
% proB The scaling returned by procrustes metheod (see doc procrustes, this will be 1)
% proC The translation returned by procrustes metheod (see doc procrustes)
% centerRow The center row determined by the user
% centerCol The center column determined by the user
% centerSlice The center slice determined by the user
% xSphere The found sphere x-locations before corrections
% ySphere The found sphere y-locations before corrections
% zSphere The found sphere z-locations before corrections
% proX The found sphere x-locations after corrections
% proY The found sphere y-locations after corrections
% proZ The found sphere z-locations after corrections
% xGnd The ground truth x-locations without corrections
% yGnd The ground truth y-locations without corrections
% zGnd The ground truth z-locations without corrections
% voxelWidth (mm) The width of each voxel
% voxelHeight (mm) The height of each voxel
% voxelLength (mm) The length of each voxel
% radiusSearch1 (mm) The distance from isocenter defining the first analysis region
% radiusSearch2 (mm) The distance from isocenter defining the second analysis region
% radiusSearch3 (mm) The distance from isocenter defining the third analysis region
% radiusThreshold1 (mm) The deviation threshold for the first analysis region
% radiusThreshold2 (mm) The deviation threshold for the second analysis region
% radiusThreshold3 (mm) The deviation threshold for the third analysis region
% centSliceImg The center slice image
% volData The volumetric dataset
% saveDeformationImg (y=1,n=0) Whether or not to save the deformation field image
%
% Output:
% volPass2D Volume of markers showing passing sphere locations for 2D deviation, in the three analysis regions 
% volFail2D Volume of markers showing failed sphere locations for 2D deviation, in the three analysis regions 
% volPass3D Volume of markers showing passing sphere locations for 3D deviation, in the three analysis regions 
% volFail3D Volume of markers showing failed sphere locations for 3D deviation, in the three analysis regions 
% volGndTruth Volume of ground truth locations corrected to match original found sphere location
% volAnalysisRegions Volume showing the different analysis regions
%
% John Ginn
% Created: 1/25/17
% Modified: 1/26/17

function[volPass2D,volFail2D,volPass3D,volFail3D,volGndTruth,...
    volAnalysisRegions] = ...
    correctGndTruth_Procrustes(proT,proB,proC,centerRow,centerCol,centerSlice,xSphere,ySphere,zSphere,...
    proX,proY,proZ,xGnd,yGnd,zGnd,...
    voxelWidth,voxelHeight,voxelLength,radiusSearch1,radiusSearch2,radiusSearch3,...
    radiusThreshold1,radiusThreshold2,radiusThreshold3,centSliceImg,volData,saveDeformationImg)
% factor for upscaling the image
scaleFactor = 3;
plotNonRotated = 0; % for debugging, whether or not to plot the non-rotated 

% extract the centermost slice
count = 0;
for step = 1:length(xGnd)
    if (round(zGnd(step)) == centerSlice)
        count = count + 1;
        xGndSmall(count) = xGnd(step);
        yGndSmall(count) = yGnd(step);
        zGndSmall(count) = zGnd(step);
        xSphereSmall(count) = xSphere(step);
        ySphereSmall(count) = ySphere(step);
        zSphereSmall(count) = zSphere(step);
        xSphereCorrSmall(count) = proX(step);
        ySphereCorrSmall(count) = proY(step);
        zSphereCorrSmall(count) = proZ(step);
    end
end
xGndCorr = xGndSmall;
yGndCorr = yGndSmall;
zGndCorr = zGndSmall;

gndCorrected = zeros(length(xGndCorr),3);
devGndCorr2D = zeros(1,length(xGndCorr));
devGndCorr3D = zeros(1,length(xGndCorr));
devSphCorr2D = zeros(1,length(xGndCorr));
devSphCorr3D = zeros(1,length(xGndCorr));
devNoCorr = zeros(1,length(xGndCorr));
radiusGndCorr = zeros(1,length(xGndCorr));
% step through the series of rotation and translation corrections
    
% for making image shapes
sphereFitPtsPass = zeros(scaleFactor*length(centSliceImg(:,1)),scaleFactor*length(centSliceImg(1,:)));
sphereFitPtsFail = sphereFitPtsPass;
grndTruthPts = sphereFitPtsPass;
grndTruthPtsNonRot = sphereFitPtsPass;

% apply the reverse of procrustes --> Y = (Z - C)T^-1
groundBeforeCorr = [xGndSmall',yGndSmall',zGndSmall'];
proCtemp = proC(1,:);
xProCsmall = proCtemp(1,1).*ones(length(xGndSmall),1);
yProCsmall = proCtemp(1,2).*ones(length(xGndSmall),1);
zProCsmall = proCtemp(1,3).*ones(length(xGndSmall),1);
proCsmall = [xProCsmall,yProCsmall,zProCsmall];
gndCorrected = (groundBeforeCorr - proCsmall)/(proT); % this does the inverse 

% test the inverse of the transform
% reverseTest = [proX,proY,proZ];
% reverseResult = (reverseTest - proC)/(proT); % this does the inverse 
% reverseX = reverseResult(:,1);
% reverseY = reverseResult(:,2);
% reverseZ = reverseResult(:,3);
% totDiff = sum([sum(xSphere - reverseX),sum(ySphere - reverseY),sum(zSphere - reverseZ)]);
for step = 1:length(xGndSmall)
    xGndCorr(step) = gndCorrected(step,1);
    yGndCorr(step) = gndCorrected(step,2);
    zGndCorr(step) = gndCorrected(step,3);
    % compare the deviation between the corrected ground truth data and the
    % original sphere location, with the corrected sphere location and
    % original ground truth
    devGndCorr2D(step) = sqrt((voxelWidth*(xGndCorr(step) - xSphereSmall(step)))^2 + ...
        (voxelHeight*(yGndCorr(step) - ySphereSmall(step)))^2);
    devGndCorr3D(step) = sqrt((voxelWidth*(xGndCorr(step) - xSphereSmall(step)))^2 + ...
        (voxelHeight*(yGndCorr(step) - ySphereSmall(step)))^2 + ...
        (voxelLength*(zGndCorr(step) - zSphereSmall(step)))^2);
    devSphCorr2D(step) = sqrt((voxelWidth*(xGndSmall(step) - xSphereCorrSmall(step)))^2 + ...
        (voxelHeight*(yGndSmall(step) - ySphereCorrSmall(step)))^2);
    devSphCorr3D(step) = sqrt((voxelWidth*(xGndSmall(step) - xSphereCorrSmall(step)))^2 + ...
        (voxelHeight*(yGndSmall(step) - ySphereCorrSmall(step)))^2 + ...
        (voxelLength*(zGndSmall(step) - zSphereCorrSmall(step)))^2);
    devNoCorr(step) = sqrt((voxelWidth*(xGndSmall(step) - xSphereSmall(step)))^2 + ...
        (voxelHeight*(yGndSmall(step) - ySphereSmall(step)))^2 + ...
        (voxelLength*(zGndSmall(step) - zSphereSmall(step)))^2);
    
    % determine the distance from isocenter
    radiusGndCorr(step) = sqrt((voxelWidth*(centerCol - xSphereSmall(step)))^2 + ...
        (voxelHeight*(centerRow - ySphereSmall(step)))^2 + ...
        (voxelLength*(centerSlice - zSphereSmall(step)))^2);
    
    
    
    % check different analysis regions and plot the results
    data = round([ySphereSmall(step) xSphereSmall(step) zSphereSmall(step)].*scaleFactor);
    shapeExtraSize = 3;
    xMarkerLoc = (data(1)-scaleFactor - shapeExtraSize):1:(data(1)+scaleFactor + shapeExtraSize); % make a symbol on image
    yMarkerLoc = (data(2)-scaleFactor - shapeExtraSize):1:(data(2)+scaleFactor + shapeExtraSize);% make a symbol on image
    if (radiusGndCorr(step) <= radiusSearch1)
        if devSphCorr2D(step)< radiusThreshold1
            % the sphere passes
            sphereFitPtsPass(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'square','thick'); % make a square
        else
            % the sphere fails
            sphereFitPtsFail(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'x','thick'); % make a x
        end
    end
    if (radiusGndCorr(step) <= radiusSearch2)&&(radiusGndCorr(step) > radiusSearch1)
        if devSphCorr2D(step)< radiusThreshold2
            % the sphere passes
            sphereFitPtsPass(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'square','thick'); % make a square
        else
            % the sphere fails
            sphereFitPtsFail(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'x','thick'); % make a x
        end
    end
    if (radiusGndCorr(step) <= radiusSearch3)&&(radiusGndCorr(step) > radiusSearch2)
        if devSphCorr2D(step)< radiusThreshold3
            % the sphere passes
            sphereFitPtsPass(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'square','thick'); % make a square
        else
            % the sphere fails
            sphereFitPtsFail(xMarkerLoc,yMarkerLoc) = makeShape(length(xMarkerLoc),'x','thick'); % make a x
        end
    end
    % Marker for fusing on the image (add on data in each slice)
    % marker{sliceNumber} = [current data; new data]
    % store location of ground truth for spheres
    if (radiusGndCorr(step) <= radiusSearch3)
        groundExtraSize = 2;
        data = round([yGndCorr(step) xGndCorr(step) zGndCorr(step)].*scaleFactor);
        xMarkerLoc = (data(1)-scaleFactor-groundExtraSize):1:(data(1)+scaleFactor+groundExtraSize); % make + symbol on image
        yMarkerLoc = (data(2)-scaleFactor-groundExtraSize):1:(data(2)+scaleFactor+groundExtraSize);% make + symbol on image
        grndTruthPts(xMarkerLoc,data(2)) = 1; % - portion of + marker for ground truth of spheres
        grndTruthPts(data(1),yMarkerLoc) = 1; % - portion of + marker for ground truth of spheres
        grndTruthPts(xMarkerLoc,data(2)-1) = 1; % - portion of + marker for ground truth of spheres
        grndTruthPts(data(1)-1,yMarkerLoc) = 1; % - portion of + marker for ground truth of spheres
        grndTruthPts(xMarkerLoc,data(2)+1) = 1; % - portion of + marker for ground truth of spheres
        grndTruthPts(data(1)+1,yMarkerLoc) = 1; % - portion of + marker for ground truth of spheres
        
        data = round([yGndSmall(step) xGndSmall(step) zGndSmall(step)].*scaleFactor);
        xMarkerLoc = (data(1)-scaleFactor-groundExtraSize):1:(data(1)+scaleFactor+groundExtraSize); % make + symbol on image
        yMarkerLoc = (data(2)-scaleFactor-groundExtraSize):1:(data(2)+scaleFactor+groundExtraSize);% make + symbol on image
        grndTruthPtsNonRot(xMarkerLoc,data(2)) = 1; % - portion of + marker for ground truth of spheres
        grndTruthPtsNonRot(data(1),yMarkerLoc) = 1; % - portion of + marker for ground truth of spheres
        grndTruthPtsNonRot(xMarkerLoc,data(2)-1) = 1; % - portion of + marker for ground truth of spheres
        grndTruthPtsNonRot(data(1)-1,yMarkerLoc) = 1; % - portion of + marker for ground truth of spheres
        grndTruthPtsNonRot(xMarkerLoc,data(2)+1) = 1; % - portion of + marker for ground truth of spheres
        grndTruthPtsNonRot(data(1)+1,yMarkerLoc) = 1; % - portion of + marker for ground truth of spheres
    end
end
%% Now do the plotting

% For debugging and checking the rotation
% avgDevGndCorr = sum(devGndCorr)/length(devGndCorr);
% avgDevSphCorr = sum(devSphCorr)/length(devSphCorr);
% avgDevNoCorr = sum(devNoCorr)/length(devNoCorr);
% 
% disp(strcat(['Average deviation correcting ground-truth: ',num2str(avgDevGndCorr)]))
% disp(strcat(['Average deviation correcting sphere locations: ',num2str(avgDevSphCorr)]))
% disp(strcat(['Average deviation without corrections: ',num2str(avgDevNoCorr)]))
% figure;
% plot(devGndCorr-devSphCorr,'.')


% upscale the image for plotting
% original information
origX = 1:length(centSliceImg(1,:));
origY = 1:length(centSliceImg(:,1));
[origMeshX,origMeshY] = meshgrid(origX,origY);
% new data
volUpscaleX = 1/scaleFactor:1/scaleFactor:length(centSliceImg(1,:));
volUpscaleY = 1/scaleFactor:1/scaleFactor:length(centSliceImg(:,1));
[upMeshX, upMeshY] = meshgrid(volUpscaleX,volUpscaleY);
% the upscaled data
SigDataUpscale = interp2(origMeshX,origMeshY,centSliceImg,...
    upMeshX,upMeshY,'cubic');



% Plot the deformation field
xSphDef = xSphereSmall.*scaleFactor;
ySphDef = ySphereSmall.*scaleFactor;
zSphDef = zSphereSmall.*scaleFactor;
xGndDef = xGndSmall.*scaleFactor;
yGndDef = yGndSmall.*scaleFactor;
zGndDef = zGndSmall.*scaleFactor;
plotDeformationField(xSphDef,ySphDef,zSphDef,xGndDef,yGndDef,zGndDef,devSphCorr3D,SigDataUpscale,...
    voxelWidth,voxelHeight,voxelLength,scaleFactor,saveDeformationImg)

sliceShift = 0;
upCenterRow = centerRow*scaleFactor;
upCenterCol = centerCol*scaleFactor;
upVoxWidth = voxelWidth/scaleFactor;
upVoxHeight = voxelHeight/scaleFactor;

[circleImg1] = imageCircle(SigDataUpscale,radiusSearch1,upVoxHeight,upVoxWidth,...
    upCenterRow,upCenterCol,sliceShift);
[circleImg2] = imageCircle(SigDataUpscale,radiusSearch2,upVoxHeight,upVoxWidth,...
    upCenterRow,upCenterCol,sliceShift);
[circleImg3] = imageCircle(SigDataUpscale,radiusSearch3,upVoxHeight,upVoxWidth,...
    upCenterRow,upCenterCol,sliceShift);




figure;
% fit marker plotting
imshow(SigDataUpscale,[])
for step = 1:3
    if step == 1
        currentImg = circleImg1;
    elseif step == 2
        currentImg = circleImg2;
    else
        currentImg = circleImg3;
    end
    rgb= [0,183,229]/255;
    blue = cat(3, rgb(1).*ones(size(currentImg)),rgb(2).*ones(size(currentImg)),rgb(3).*ones(size(currentImg)));
    hold on
    hBlue = imshow(blue);
    hold off
    set(hBlue, 'AlphaData', currentImg) % make color sheet only show markers
end
% cat(3,r,g,b)
green = cat(3, zeros(size(sphereFitPtsPass)),ones(size(sphereFitPtsPass)), zeros(size(sphereFitPtsPass)));
% bright = 0.9; 
% white version
% green = bright.*cat(3, ones(size(sphereFitPtsPass)),ones(size(sphereFitPtsPass)), ones(size(sphereFitPtsPass)));
hold on 
hGreen = imshow(green); 
hold off
set(hGreen, 'AlphaData', sphereFitPtsPass) % make color sheet only show markers
% the spheres that failed
red = cat(3, ones(size(sphereFitPtsFail)),zeros(size(sphereFitPtsFail)), zeros(size(sphereFitPtsFail)));
% white version
% red = bright.*cat(3, ones(size(sphereFitPtsFail)),ones(size(sphereFitPtsFail)), ones(size(sphereFitPtsFail)));
hold on 
hRed = imshow(red); 
hold off
set(hRed, 'AlphaData', sphereFitPtsFail) % make color sheet only show markers

% rotated ground truth marker plotting
% cat(3,r,g,b)
rgb= [64,224,208]/255;
white = cat(3, ones(size(grndTruthPts)),ones(size(grndTruthPts)),ones(size(grndTruthPts)));
hold on 
hBlue = imshow(white); 
hold off
set(hBlue, 'AlphaData', grndTruthPts) % make color sheet only show markers
% title('Spatial Integrity Phantom Center Slice','FontSize',20)
xlabel('x-position (mm)','FontSize',20)
ylabel('z-position (mm)','FontSize',20)

% non-rotated ground truth marker plotting
% cat(3,r,g,b)
if plotNonRotated == 1
    rgb= [255,128,0]/255;
    white = cat(3, rgb(1).*ones(size(grndTruthPtsNonRot)),rgb(2).*ones(size(grndTruthPtsNonRot)),rgb(3).*ones(size(grndTruthPtsNonRot)));
    hold on
    hBlue = imshow(white);
    hold off
    set(hBlue, 'AlphaData', grndTruthPtsNonRot) % make color sheet only show markers
    % title('Spatial Integrity Phantom Center Slice','FontSize',20)
    xlabel('x-position (mm)','FontSize',20)
    ylabel('z-position (mm)','FontSize',20)
end

% custom axes to show distances (the locations of the axes)
numOfTicks = 10;
xAxisSpacing = floor(length(SigDataUpscale(1,:))/numOfTicks);
yAxisSpacing = floor(length(SigDataUpscale(:,1))/numOfTicks);

% xAxisLocation = linspace(1,length(SigDataUpscale(1,:)),numOfTicks);
% yAxisLocation = linspace(1,length(SigDataUpscale(:,1)),numOfTicks);
xAxisLocation = 1:xAxisSpacing:(xAxisSpacing*numOfTicks+1);
yAxisLocation = 1:yAxisSpacing:(yAxisSpacing*numOfTicks+1);
% the values on the axes
xAxisValue = voxelWidth/scaleFactor.*(xAxisLocation - round(median(xAxisLocation)));
yAxisValue = voxelHeight/scaleFactor.*(yAxisLocation - round(median(yAxisLocation)));
yAxisLocation(length(yAxisLocation)) = length(SigDataUpscale(:,1)); % special case for plotting for paper
axis on
axisHandle = gca;
set(gca,'XTickMode','manual')
set(gca,'YTickMode','manual')
set(gca,'XTick',xAxisLocation)
set(gca,'YTick',yAxisLocation)
set(gca,'XTickLabel',xAxisValue)
set(gca,'YTickLabel',yAxisValue)
set(gca,'FontSize',16)




%% Repeat the deviation calculation for the entire volume for plotting
devSphCorrVol2D = zeros(1,length(xGnd));
devSphCorrVol3D = devSphCorrVol2D;
radiusGndCorrVol = devSphCorrVol2D;
gndCorrectedAll = zeros(length(xGnd),3);


xGndCorrAll = xGnd;
yGndCorrAll = yGnd;
zGndCorrAll = zGnd;

% make sure array is correctly oriented for correction
if size(proX,1) == 1;
    proX = proX';
    proY = proY';
    proZ = proZ';
end
groundBeforeCorr = [proX,proY,proZ];
xProClarge = proCtemp(1,1).*ones(length(proX),1);
yProClarge = proCtemp(1,2).*ones(length(proX),1);
zProClarge = proCtemp(1,3).*ones(length(proX),1);
proClarge = [xProClarge,yProClarge,zProClarge];
gndCorrectedAll = (groundBeforeCorr - proClarge)/(proT); % this does the inverse 


volPass2D = zeros(length(volData(:,1,1)),length(volData(1,:,1)),length(volData(1,1,:)));
volFail2D = volPass2D;
volPass3D = volPass2D;
volFail3D = volPass2D;
volGndTruth = volPass2D;
volAnalysisRegions = volPass2D;
for step = 1:length(xGnd)    
    xGndCorrAll(step) = gndCorrectedAll(step,1);
    yGndCorrAll(step) = gndCorrectedAll(step,2);
    zGndCorrAll(step) = gndCorrectedAll(step,3);
    
    % compare the deviation between the corrected ground truth data and the
    % original sphere location, with the corrected sphere location and
    % original ground truth
    devSphCorrVol2D(step) = sqrt((voxelWidth*(xGnd(step) - proX(step)))^2 + ...
        (voxelHeight*(yGnd(step) - proY(step)))^2);
    devSphCorrVol3D(step) = sqrt((voxelWidth*(xGnd(step) - proX(step)))^2 + ...
        (voxelHeight*(yGnd(step) - proY(step)))^2 + ...
        (voxelLength*(zGnd(step) - proZ(step)))^2);
    
    % determine the distance from isocenter
    radiusGndCorrVol(step) = sqrt((voxelWidth*(centerCol - proX(step)))^2 + ...
        (voxelHeight*(centerRow - proY(step)))^2 + ...
        (voxelLength*(centerSlice - proZ(step)))^2);
    

    % check different analysis regions and plot the results
    data = round([ySphere(step) xSphere(step) zSphere(step)]);
    xMarkerLoc = (data(1)-1 - 1):1:(data(1)+1 + 1); % make a symbol on image
    yMarkerLoc = (data(2)-1 - 1):1:(data(2)+1 + 1);% make a symbol on image
    if (radiusGndCorrVol(step) <= radiusSearch1)
        % 2D volume
        if devSphCorrVol2D(step)< radiusThreshold1
            % the sphere passes
            volPass2D(xMarkerLoc,yMarkerLoc,data(3)) =...
                makeShape(length(xMarkerLoc),'square','thin'); % make a square
        else
            % the sphere fails
            volFail2D(xMarkerLoc,yMarkerLoc,data(3)) = ...
                makeShape(length(xMarkerLoc),'x','thin'); % make a x
        end
        
        % 3D volume
        if devSphCorrVol3D(step)< radiusThreshold1
            % the sphere passes
            volPass3D(xMarkerLoc,yMarkerLoc,data(3)) =...
                makeShape(length(xMarkerLoc),'square','thin'); % make a square
        else
            % the sphere fails
            volFail3D(xMarkerLoc,yMarkerLoc,data(3)) = ...
                makeShape(length(xMarkerLoc),'x','thin'); % make a x
        end
    end
    if (radiusGndCorrVol(step) <= radiusSearch2)&&(radiusGndCorrVol(step) > radiusSearch1)
        % 2D volume
        if devSphCorrVol2D(step)< radiusThreshold2
            % the sphere passes
            volPass2D(xMarkerLoc,yMarkerLoc,data(3)) =...
                makeShape(length(xMarkerLoc),'square','thin'); % make a square
        else
            % the sphere fails
            volFail2D(xMarkerLoc,yMarkerLoc,data(3)) = ...
                makeShape(length(xMarkerLoc),'x','thin'); % make a x
        end
        
        % 3D volume
        if devSphCorrVol3D(step)< radiusThreshold2
            % the sphere passes
            volPass3D(xMarkerLoc,yMarkerLoc,data(3)) =...
                makeShape(length(xMarkerLoc),'square','thin'); % make a square
        else
            % the sphere fails
            volFail3D(xMarkerLoc,yMarkerLoc,data(3)) = ...
                makeShape(length(xMarkerLoc),'x','thin'); % make a x
        end
    end
    if (radiusGndCorrVol(step) <= radiusSearch3)&&(radiusGndCorrVol(step) > radiusSearch2)
        % 2D volume
        if devSphCorrVol2D(step)< radiusThreshold3
            % the sphere passes
            volPass2D(xMarkerLoc,yMarkerLoc,data(3)) =...
                makeShape(length(xMarkerLoc),'square','thin'); % make a square
        else
            % the sphere fails
            volFail2D(xMarkerLoc,yMarkerLoc,data(3)) = ...
                makeShape(length(xMarkerLoc),'x','thin'); % make a x
        end
        
        % 3D volume
        if devSphCorrVol3D(step)< radiusThreshold3
            % the sphere passes
            volPass3D(xMarkerLoc,yMarkerLoc,data(3)) =...
                makeShape(length(xMarkerLoc),'square','thin'); % make a square
        else
            % the sphere fails
            volFail3D(xMarkerLoc,yMarkerLoc,data(3)) = ...
                makeShape(length(xMarkerLoc),'x','thin'); % make a x
        end
    end
    % Marker for fusing on the image (add on data in each slice)
    % marker{sliceNumber} = [current data; new data]
    % store location of ground truth for spheres
    data = round([yGndCorrAll(step) xGndCorrAll(step) zGndCorrAll(step)]);
    xMarkerLoc = (data(1)-1):1:(data(1)+1); % make + symbol on image
    yMarkerLoc = (data(2)-1):1:(data(2)+1);% make + symbol on image
    volGndTruth(xMarkerLoc,data(2),data(3)) = 1; % - portion of + marker for ground truth of spheres
    volGndTruth(data(1),yMarkerLoc,data(3)) = 1; % - portion of + marker for ground truth of spheres
    % these make the '+' thicker
%     volGndTruth(xMarkerLoc,data(2)-1,data(3)) = 1; % - portion of + marker for ground truth of spheres
%     volGndTruth(data(1)-1,yMarkerLoc,data(3)) = 1; % - portion of + marker for ground truth of spheres
%     volGndTruth(xMarkerLoc,data(2)+1,data(3)) = 1; % - portion of + marker for ground truth of spheres
%     volGndTruth(data(1)+1,yMarkerLoc,data(3)) = 1; % - portion of + marker for ground truth of spheres
    
end

% make the volume of circles here...
for step = 1:length(volAnalysisRegions(1,1,:))
    % calculate the distance to the slice from isocenter
    sliceDist = voxelLength*(centerSlice - step);
    imageData = volAnalysisRegions(:,:,step);
    % radiusSearch1 image
    [circleImg1] = imageCircle(imageData,radiusSearch1,voxelHeight,voxelWidth,...
    centerRow,centerCol,sliceDist);
    volAnalysisRegions(:,:,step) = circleImg1;
    % radiusSearch2 image
    [circleImg2] = imageCircle(imageData,radiusSearch2,voxelHeight,voxelWidth,...
    centerRow,centerCol,sliceDist);
    volAnalysisRegions(:,:,step) = volAnalysisRegions(:,:,step) + circleImg2;
    % radiusSearch3 image
    [circleImg3] = imageCircle(imageData,radiusSearch3,voxelHeight,voxelWidth,...
    centerRow,centerCol,sliceDist);
    volAnalysisRegions(:,:,step) =  volAnalysisRegions(:,:,step) + circleImg3;
    
end

end



