% Master script to find the location of the spheres in the spatial integrity
% phantom produced by Integrated Medical Technologies. This script can analyze
% either a CT or MRI scan of the spatial integrity phantom. The parameters 
% below need to be changed manually to suit individual scan.
%
% NOTE: In general, x refers to the in-plane left-right location in the image.
% (typically medial-lateral) y refers to the in-plane up-down location in the 
% image (typically anterior-posterior) and z refers to the out-of-plane location
% (typically inferior-superior). The labels in the plots however are often
% reversed in the plots in order to agree with the clinical definition of the axes.
% The easiest way to decipher which definition is being used is to use the
% asymmetries of the phantom.
%
% All the dicom images of the spatial integrity phantom must be in a folder
% without any additional files in that folder
%
% The user must adjust all parameters according to their acquisition parameters
% in order for the analysis code to work appropriately
% 
% John Ginn
% Created: 6/23/16
% Modified: 2/7/17

clear all
close all

%% Parameters
gantryAngle = 0; % (degrees) angle the gantry was at during scanning or string 'NA' etc.
voxelHeight = 1.5; % (mm) y-dimension (anterior-posterior)
voxelWidth = 1.5; % (mm) x-dimension (medial-lateral)
voxelLength = 3; % (mm) z-dimension (inferior-superior)
corrThreshold = 0.70; % the minimum correlation threshold
radiusSearch1 = 50;  % (mm) distance for searching for spheres within 100 mm of isocenter
radiusThreshold1 = 1; % (mm) the passing criteria for spheres within 100 mm of isocenter
radiusSearch2 = 100;  % (mm) distance for searching for spheres within 175 mm of isocenter
radiusThreshold2 = 1; % (mm) the passing criteria for spheres within 175 mm of isocenter
radiusSearch3 = 175;  % (mm) distance for searching for spheres within 175 mm of isocenter
radiusThreshold3 = 2; % (mm) the passing criteria for spheres within 175 mm of isocenter
radiusSearch4 = 175;  % (mm) distance for searching for spheres within 175 mm of isocenter
radiusThreshold4 = 1; % (mm) the passing criteria for spheres within 175 mm of isocenter
allThreshold = 2; % (mm) threshold for passing in entire volume
MRorCT = 'mr'; % ('mr' or 'ct') Determines whether the sphere template contrast exhibits a CT or MR scan 
rotCT = 'y'; % ('y' or 'n') Whether or not to rotate the CT image depending on scanning orientation
% if orientation was the same as in the MR, do not rotate (only affects analysis when MRorCT == 'ct')
ctImageExport = 'sagittal'; % ('sagittal' or 'original') Whether the CT images
% were exported from MIM using the sagittal option or just the original
% orientation (Use this parameter only when the phantom was not oriented in
% the same manner as on the MRI scanner, and when a rotation is applied)
% (only affects analysis when (MRorCT == 'ct')&&(rotCT == 'y'))

% Options
% do you need to obtain the center slice, center column and center row location?
obtainLocations = 1; % y = 1,n = 0
computeInParallel = 1; % y = 1,n = 0 (whether or not to compute the locations of the spheres using the parallel toolbox)
% WARNING: this will overwrite images saved in this direcory with the same name
saveExcelFiles = 0; % y = 1, n = 0 save the data to an excel file
saveImages = 0; % y = 1, n = 0 save the images at a high resolution 
interpolateSlice = 0;  % y = 1, n = 0 interpolate single image set slice thickness
oldSliceThickness = voxelLength; % (mm) if super interpolateSlice option is on, the old slice thickness
newSliceThickness = 1.5; % if interpolateSlice option is on, the new slice thickness
analysisType = 'weighted';  % determines the data used in the analysis
% 'non-upscale' for non upscaled correlation coefficient location
% 'upscale' for upscaled correlation coefficient location
% 'weighted' method that calculates sphere locations based on the
% weighted sum of the signal.

% Other Options and Parameters Rarely Changed 
searchType = 'upscale'; 
% determines the method of finding the spheres
% 'upscale' if you want to upscale the data for the analysis
% 'non-upscale' locations based off the non-upscaled correlation coefficient
% the starting location for the weighted sum will be determined by the
% location found using either the non-upscaled or upscaled correlation coefficient
if strcmp(searchType,'upscale'); 
    % y = 1, n = 0 construct an upscaled dataset and compute the correlation coefficent to find the location of the sphere
    upscale = 1;
else
    upscale = 0;
end

% parameters
searchDist = 6; % (mm) distance in each direction to search for the sphere (orig corr coeff)
searchDistWeight = 6; % (mm) distance in each direction to use to calulate weighted sum location
binSize = 20; % (mm) The size of the bins for plotting the deviation
% spacing between the spheres (FOUND USING CT OF PHANTOM)
xSpacing = 16.01; % (mm) new x-direction spacing 
ySpacing = 16.05; % (mm) new y-direction spacing 
zSpacing = 16.08; % (mm) new z-direction spacing

% options
phantomSim = 0; % y = 1, n = 0 Phantom Sim Debugging? This will overwrite any data obtained from image
phantomSignalAmp = 130;
phantomSimNoise = 0;% y = 1, n = 0 whether or not to add noise to the phantom simulation
phantomSimNoiseMean = 0.00; % mean of gaussian noise added to phantom simulation
phantomSimVariance = 0.1;
sphereBoundary = 1; % y = 1,n = 0 make sphere model have a boundary of zeros around it
plotPlaneFit = 0; % y = 1, n = 0 plot the fit of the planes for translational correction
% old method and option below, see "Import data from MR images" section below or
% user manual for more
% superRes = 0; % y = 1, n = 0 whether or not you need to create a super-resolution volume of two image datasets
%
% change if you already know what these values are,
% otherwise they will be overwritten when you run the code
if obtainLocations == 0
    centerSlice = 146; % the index for the center slice of spheres
    centerCol = 168; % location of center column of spheres
    centerRow = 152; % index of the center row
end
%
%
% Phantom Simulation Parameters
if phantomSim == 1;
    % NOTE: use an odd number of spheres so there is a "center slice" along
    % each of the principle axes
    nSpheres = 5; % the number of spheres in each direction
    obtainLocations = 0; % the locations are calculated, don't obtain them
    phantomSimRotationTrans = 0.2; % (degrees) rotation of data
    phantomSimRotationCoronal = 0.1; % (degrees) rotation of data 
    phantomSimRotationSagittal = 0.1; % (degrees) rotation of data
    phantomSimRotDirTrans = 0; % clockwise = 1, counter-clockwise = 0
    phantomSimRotDirCoronal = 1; % clockwise = 1, counter-clockwise = 0
    phantomSimRotDirSagittal = 1; % clockwise = 1, counter-clockwise = 0
    phantomSimShiftX = 0.1; % shift the data
    phantomSimShiftY = -0.2; % shift the data
    phantomSimShiftZ = 0.3; % shift the data
    MRorCT = 'mr';
end

if (strcmp(MRorCT,'ct') && strcmp(rotCT,'y'))
    % rotate the CT image if necessary depending on scanning orientation
    voxelHeightTemp = voxelHeight; % (mm) y-dimension (anterior-posterior)
    voxelWidthTemp = voxelWidth; % (mm) x-dimension (medial-lateral)
    voxelLengthTemp = voxelLength; % (mm) z-dimension (inferior-superior)
    % change the voxel dimensions
    voxelHeight = voxelHeightTemp; % (mm) y-dimension (anterior-posterior)
    voxelWidth = voxelLengthTemp; % (mm) x-dimension (medial-lateral)
    voxelLength = voxelWidthTemp; % (mm) z-dimension (inferior-superior)
end


% phantom information
sphereDiameter = 8; % (mm)
distBtwnSpheres = 16; % (mm) The distance between the spheres in each direction (from phantom specifications)
xIndPerSph = xSpacing/voxelWidth; % number of pixels in x-direction between the spheres
yIndPerSph = ySpacing/voxelHeight; % number of pixels in y-direction between the spheres
zIndPerSph = zSpacing/voxelLength; % number of pixels in z-direction between the spheres
% number of spheres in each row from bottom to top
spheresInRow = [29 29 29 29 29 29 29 29 27 27 27 27 25 25 23 23 21 19 17 13 7]; 
centerRowNum = round(median(1:length(spheresInRow)));
nRows = 21; % there are 21 rows of spheres
nSlices = 9; % there are 9 "slices" of spheres
nTotData = nSlices.*sum(spheresInRow); % total number of spheres

%% Import data from MR images
if phantomSim ~= 1
    % obtain the data from the images
    currentDir = pwd;
        disp('Please select folder containing the dicom images')
        dataDir = uigetdir;
        linesToSkip = 3;
        [fileNames, fileData, fileMap, fileInfo] = loadImages(dataDir,currentDir,linesToSkip);
        % Construct the volumetric signal
        volData = zeros(size(fileData{1},1),size(fileData{1},2),length(fileData));
        for step = 1:length(fileData);
            volData(:,:,step) = fileData{step};
        end
        % interpolate a single image set to a new slice thickness
        if interpolateSlice == 1;
            xDim = length(volData(1,:,1));
            yDim = length(volData(:,1,1));
            zDim = length(volData(1,1,:));
            upscaleFactor = oldSliceThickness/newSliceThickness; % upscale the data
            [X_ns,Y_ns,Z_ns]=meshgrid((1:xDim)*voxelWidth,(1:yDim)*voxelHeight,(1:zDim)*oldSliceThickness);
            [Xq,Yq,Zq]=meshgrid((1:xDim)*voxelWidth,(1:yDim)*voxelHeight,(1:zDim*upscaleFactor)*newSliceThickness);
            volDataTemp = interp3(X_ns,Y_ns,Z_ns,single(volData),Xq,Yq,Zq);    
            % remove bounds containing NaN values
            volData = zeros(length(volDataTemp(:,1,1)),...
                length(volDataTemp(1,:,1)),length(volDataTemp(1,1,:))); % reset volData
            fileMap = cell(1,length(volDataTemp(1,1,:)));
            for step = 2:(length(volDataTemp(1,1,:)) - 1)
                volData(:,:,step) = volDataTemp(:,:,step); 
            end
            % update voxel length
            voxelLength = newSliceThickness;
            zIndPerSph = zSpacing/voxelLength; % number of pixels in z-direction between the spheres
        end
    if (strcmp(MRorCT,'ct') && strcmp(rotCT,'y'))
        % rotate the CT image if necessary depending on scanning orientation
        nonInvertData = volData; % store for checking later
        dimX = length(nonInvertData(1,:,1));
        dimY = length(nonInvertData(:,1,1));
        dimZ = length(nonInvertData(1,1,:));
        % volData = zeros(dimY,dimZ,dimX); % reset the volume data rotate original scan images
        volData = zeros(dimX,dimY,dimZ); % reset the volume data rotate sagittal exported images from mim
        fileMap = cell(dimX,1); % change the dimension of the file map for GUI
        for stepX = 1:dimX
            for stepY = 1:dimY
                for stepZ = 1:dimZ
                    % rotate the data
                    if strcmp(ctImageExport,'sagittal')                       
                        % rotate sagittal exported images from mim
                        volData(stepX,stepY,stepZ) = nonInvertData(stepY,stepX,stepZ);
                    else
                        % rotate original scan images
                        volData(stepY,stepZ,stepX) = volDataTemp(stepY,stepX,stepZ);
                    end
                end
            end
        end
        % invert the image
        [volData] = invertCT(volData);
    end
    disp('MRorCT has been changed to mr for sphere template')
    MRorCT = 'mr';
    %
    % Plot the volume extracted from images or simulation and obtain locations
    % of the center slice, row and column if necessary
    % store variables for GUI
    guiData{1} = volData;
    guiData{2} = fileMap;
    % display the GUI
    UserInputGUI2(guiData);
    if obtainLocations == 1;
        disp('Please find the center of the sphere')
        centerSlice = input('Enter center slice index: '); % the index for the center slice of spheres
        centerCol = input('Enter center column index [x]: '); % location of center column of spheres
        centerRow = input('Enter center row index, the phantom has 21 rows [y]: '); % index of the center row
    end
end

% the calculated sphere slice position
if phantomSim ~= 1;
    sliceCalc = zeros(1,nSlices);
    sliceGnd = sliceCalc;
    for step = 1:nSlices
        if step < 5;
            % move downward inverted because of matlab indexing in plot
            sliceCalc(step) = centerSlice - round(zIndPerSph*step);
            sliceGnd(step) = centerSlice - zIndPerSph*step;
        elseif step == 5; % the center slice of spheres
            sliceCalc(step) = centerSlice;
            sliceGnd(step) = centerSlice;
        else
            % move upward, inverted because of indexing in plot
            % the count should start at 1 and move up (thus subtract 5)
            slicePos = (step - 5);
            sliceCalc(step) = centerSlice + round(zIndPerSph*slicePos);
            sliceGnd(step) = centerSlice + zIndPerSph*slicePos;
        end
    end
    sliceCalc = sort(sliceCalc);
    sliceGnd = sort(sliceGnd);
else
    % phantom simulation testing
    distBtwnSpheres = 16; % (mm) The distance between the spheres in each direction
    xSpacing = distBtwnSpheres; % (mm) new x-direction spacing
    ySpacing = distBtwnSpheres; % (mm) new y-direction spacing
    zSpacing = distBtwnSpheres; % (mm) new z-direction spacing
    xIndPerSph = xSpacing/voxelWidth; % number of pixels in x-direction between the spheres
    yIndPerSph = ySpacing/voxelHeight; % number of pixels in y-direction between the spheres
    zIndPerSph = zSpacing/voxelLength; % number of pixels in z-direction between the spheres

    [volData, sliceCalc, sliceGnd] = ...
        makePhantom(voxelHeight, voxelWidth, voxelLength, sphereDiameter,nSpheres,sphereBoundary);
    fileMap = cell(1,length(volData));
    % add noise to the simulation
    if phantomSimNoise == 1;
        for stepPhantom = 1:length(volData(1,1,:));
            phantomImg = volData(:,:,stepPhantom); % extract individual images
            % add the noise
            volData(:,:,stepPhantom) = ...
                imnoise(phantomImg,'gaussian',phantomSimNoiseMean,phantomSimVariance);
        end
    end
    centerSlice = median(sliceCalc);
    centerCol = median(1:length(volData(1,:,1)));
    centerRow = median(1:length(volData(:,1,1)));
    nSlices = nSpheres;
    % simulated phantom
    spheresInRow = ones(1,nSpheres).*nSpheres;
    nRows = nSpheres; % there are 21 rows of spheres
    centerRowNum = round(median(1:nRows)); % the center row number
    sliceCalc = sort(sliceCalc);
    sliceGnd = sort(sliceGnd);
    centerCol = round(centerCol);
    centerRow = round(centerRow);
    centerSlice = round(centerSlice);
    guiData{1} = volData;
    guiData{2} = fileMap;
    % display the GUI
    UserInputGUI2(guiData);
end


% Correlation fitting to find the spheres

% make sphere for correlation fitting
sphereVol = makeSphere(voxelHeight, voxelWidth, voxelLength, sphereDiameter,sphereBoundary,'mr');
% check the volume of the sphere template for debugging
% sphereTempGUI{1} = sphereVol;
% emptyArray = cell(1,length(sphereVol(1,1,:)));
% sphereTempGUI{2} = emptyArray; % the map for the image scale
% UserInputGUI2(sphereTempGUI)

%% Scan the entire phantom for the sphere locations
tic
if computeInParallel == 1
    [coordSphereData, coordSpherePlot, optCorrelation, sphereGroundTruth,...
        weightCoord, weightCoordPlot, coordSphereDataUp,coordSpherePlotUp, optCorrelationUp,...
        volSearchRegion,volAllGndTruth,volWeightSum] = ...
        PhantomScanParallel(volData, sphereVol, voxelHeight,voxelWidth, voxelLength,sliceCalc,...
        centerCol,centerRow,centerSlice,sliceGnd,nRows,nSlices,spheresInRow,centerRowNum,searchDist,searchDistWeight,upscale,sphereBoundary,MRorCT,...
        xSpacing,ySpacing,zSpacing);
else
    [coordSphereData, coordSpherePlot, optCorrelation, sphereGroundTruth,...
        weightCoord, weightCoordPlot, coordSphereDataUp,coordSpherePlotUp, optCorrelationUp,...
        volSearchRegion,volAllGndTruth,volWeightSum] = ...
        PhantomScan(volData, sphereVol, voxelHeight,voxelWidth, voxelLength,sliceCalc,...
        centerCol,centerRow,centerSlice,sliceGnd,nRows,nSlices,spheresInRow,centerRowNum,searchDist,searchDistWeight,upscale,sphereBoundary,MRorCT,...
        xSpacing,ySpacing,zSpacing);
end

finalTime = toc;
disp(strcat(['Analysis Search Duration: ',num2str(finalTime/60),' (min)']))


% check the search regions for all the spheres
guiDataSearch{1} = volData;
guiDataSearch{2} = fileMap;
guiDataSearch{3} = volWeightSum;
guiDataSearch{4} = volAllGndTruth;
CheckSearchRegion(guiDataSearch);

if phantomSim == 1
    % the location of the ground truth for the simulated data may vary
    % depending on the voxel size because the center of each sphere is
    % assigned an integer value. Round the data to determine the ground
    % truth locations of the spheres.
    for step = 1:length(sphereGroundTruth)
        sphereGroundTruth{step} = round(sphereGroundTruth{step});
    end
end
%% construct a volume of the sphere locations for the image and remove data 

% sort the data into individual slices (required to fit a plane to the data
% to obtain the rotation of the phantom setup). Extract out the approprate
% analysis data.
if strcmp(analysisType,'upscale') == 1
    if strcmp(searchType,'upscale')
        % upscale analysis has occurred
        [sortedSphere, sortedGndTruth] = sortSlices(coordSphereDataUp, sphereGroundTruth);
        corrData = optCorrelationUp;
    else
       error('The upscaled search was not performed, change the analysis type') 
    end
elseif strcmp(analysisType,'non-upscale') == 1
    % non-upscaled correlation coefficient results
    [sortedSphere, sortedGndTruth] = sortSlices(coordSpherePlot, sphereGroundTruth);
    corrData = optCorrelation;
else
    % the weighted sum method
    [sortedSphere, sortedGndTruth] = sortSlices(weightCoordPlot, sphereGroundTruth);
    corrData = optCorrelation;
end
corrDataArray = zeros(1,length(corrData));
for step = 1:length(corrData)
   corrDataArray(step) = corrData{step}; 
end

% initialize arrays
sphereImg = zeros(length(volData(:,1,1)),length(volData(1,:,1)),length(volData(1,1,:)));
sphereFitPts = sphereImg; % store just the location of the spheres
grndTruthPts = sphereImg; % the location of the "ground truth" points
sphereVolX = length(sphereVol(1,:,1)); % x-dimension of the sphere
sphereVolY = length(sphereVol(:,1,1)); % y-dimension of the sphere
sphereVolZ = length(sphereVol(1,1,:)); % z-dimension of the sphere
% count the data that actually passes the correlation coeff. threshold
countPass = 0;
% count data in each slice that passes
countSlicePass = 0;
stepCount = 1; % total step counter
countFail = 0; % data that fail
% step through the slices
for stepSlice = 1:size(sortedSphere,2)
    % step through the spheres in each slice
    for stepSphere = 1:size(sortedSphere,1)
        if corrData{stepCount} < corrThreshold
            SphereLocData = coordSpherePlot{stepCount};
            GroundTruthFinal = sphereGroundTruth{stepCount};
            % if correlation is less than the minimum correlation threshold,
            % do not plot the sphere
            countFail = countFail + 1;
            xSphereFail(countFail) = SphereLocData(2); 
            ySphereFail(countFail) = SphereLocData(1); 
            zSphereFail(countFail) = SphereLocData(3); 
            xGndTruthFail(countFail) = GroundTruthFinal(2); 
            yGndTruthFail(countFail) = GroundTruthFinal(1);
            zGndTruthFail(countFail) = GroundTruthFinal(3);
        else
            % count the data that actually passes the correlation coeff.
            % threshold
            countPass = countPass + 1;
            % correlation threshold met
            % store data for plotting
            if mod(sphereVolX,2) == 0; % x-dimension of sphere is even
                sphRangeX = sphereVolX/2;
            else % x-dimension of sphere is odd
                sphRangeX = (sphereVolX-1)/2;
            end
            if mod(sphereVolY,2) == 0; % y-dimension of sphere is even
                sphRangeY = sphereVolY/2;
            else % y-dimension of sphere is odd
                sphRangeY = (sphereVolY-1)/2;
            end
            if mod(sphereVolZ,2) == 0; % z-dimension of sphere is even
                sphRangeZ = sphereVolZ/2;
            else % z-dimension of sphere is odd
                sphRangeZ = (sphereVolZ-1)/2;
            end
            data = coordSpherePlot{stepCount};
            xLoc = (data(2)-sphRangeX):1:(data(2)+sphRangeX); 
            yLoc = (data(1)-sphRangeY):1:(data(1)+sphRangeY); 
            zLoc = (data(3)-sphRangeZ):1:(data(3)+sphRangeZ); 
            xMarkerLoc = (data(1)-1):1:(data(1)+1); % make + symbol on image
            yMarkerLoc = (data(2)-1):1:(data(2)+1);% make + symbol on image
            sphereImg(yLoc,xLoc,zLoc) = sphereVol; % store the sphere volume
            sphereFitPts(xMarkerLoc,data(2),data(3)) = 1; % - portion of + marker for fit of spheres
            sphereFitPts(data(1),yMarkerLoc,data(3)) = 1; % | portion of + marker for fit of spheres
            % store final data without spheres that did not meet correlation
            % requirement
            % data for sphere agreement calculation
            if strcmp(analysisType,'upscale') == 1 % upscale analysis has occurred
                SphereLocData = coordSphereDataUp{stepCount};
                optCorrelationFinal{countPass} = optCorrelationUp{stepCount};
            elseif strcmp(analysisType,'non-upscale') == 1
                % non-upscaled correlation coefficient results
                SphereLocData = coordSpherePlot{stepCount};
                optCorrelationFinal{countPass} = corrData{stepCount};
            else
                % the weighted sum results (uses correlation coefficient
                % from non-upscaled data to remove spheres with a correlation 
                % coefficient less than the specified threshold)
                SphereLocData = weightCoordPlot{stepCount};
                optCorrelationFinal{countPass} = corrData{stepCount};
            end
            xSphere(countPass) = SphereLocData(2); % sphere location x-component
            ySphere(countPass) = SphereLocData(1); % sphere location y-component
            zSphere(countPass) = SphereLocData(3); % sphere location z-component
            % store location of ground truth for spheres
            data = sphereGroundTruth{stepCount};
            xMarkerLoc = round(data(1)-1):1:(data(1)+1); % make + symbol on image
            yMarkerLoc = round(data(2)-1):1:(data(2)+1);% make + symbol on image
            grndTruthPts(xMarkerLoc,round(data(2)),round(data(3))) = 1; % - portion of + marker for ground truth of spheres
            grndTruthPts(round(data(1)),yMarkerLoc,round(data(3))) = 1; % - portion of + marker for ground truth of spheres
            % ground truth points that did not meet correlation requirement
            GroundTruthFinal = sphereGroundTruth{stepCount};
            xGndTruthFinal(countPass) = GroundTruthFinal(2);
            yGndTruthFinal(countPass) = GroundTruthFinal(1);
            zGndTruthFinal(countPass) = GroundTruthFinal(3);
            % store the slice data that passes
            countSlicePass = countSlicePass + 1;
            sortedSpherePass{countSlicePass,stepSlice} = ...
                sortedSphere{stepSphere,stepSlice};
            sortedGndTruthPass{countSlicePass,stepSlice}  = ...
                sortedGndTruth{stepSphere,stepSlice};
        end
        stepCount = stepCount + 1; % add one step
    end
    % the number of rows in each slice (used later)
    sliceSpheres(stepSlice) = countSlicePass;
    countSlicePass = 0; % reset the count in the slice of spheres that pass
end

% 3d plotting of the ground truth to find rotation of phantom
% extract out a single slice for fitting to correct the rotation
fittingSlice = round(median(1:nSlices)); % this is the middle slice
for step = 1:sliceSpheres(fittingSlice)
    GndTruthData = sortedGndTruthPass{step,fittingSlice};
    xGndTruthSlice(step,1) = GndTruthData(2); % ground truth x-component
    yGndTruthSlice(step,1) = GndTruthData(1); % ground truth y-component
    zGndTruthSlice(step,1) = GndTruthData(3); % ground truth z-component
    SphereLocData = sortedSpherePass{step,fittingSlice};
    xSphereRotSlice(step,1) = SphereLocData(2); % sphere location x-component
    ySphereRotSlice(step,1) = SphereLocData(1); % sphere location y-component
    zSphereRotSlice(step,1) = SphereLocData(3); % sphere location z-component
end

% phantom sim data debugging, add shift and roation to data
if phantomSim == 1 % y = 1, n = 0;
    % shift and rotate the data that is used for fitting 
    [xSphereRotSlice, ySphereRotSlice, zSphereRotSlice] = ...
        rotatePhantomSim(xSphereRotSlice,ySphereRotSlice,zSphereRotSlice,...
    phantomSimRotationTrans,phantomSimRotationCoronal,phantomSimRotationSagittal,...
    phantomSimRotDirTrans,phantomSimRotDirCoronal,...
    phantomSimRotDirSagittal,centerRow,centerCol,centerSlice,...
    phantomSimShiftX,phantomSimShiftY,phantomSimShiftZ,0);
    % ensure that it is in correct formatting for the fitting
    xSphereRotSlice = xSphereRotSlice';
    ySphereRotSlice = ySphereRotSlice';
    zSphereRotSlice = zSphereRotSlice';
    % shift and rotate all the data
    [xSphere, ySphere, zSphere] = ...
        rotatePhantomSim(xSphere,ySphere,zSphere,...
    phantomSimRotationTrans,phantomSimRotationCoronal,phantomSimRotationSagittal,...
    phantomSimRotDirTrans,phantomSimRotDirCoronal,...
    phantomSimRotDirSagittal,centerRow,centerCol,centerSlice,...
    phantomSimShiftX,phantomSimShiftY,phantomSimShiftZ,0);
    % plot after the rotation has been applied
    plotSym3D(xSphere,ySphere,zSphere,xGndTruthFinal,yGndTruthFinal,zGndTruthFinal)
end

% used for plottingcoma
for step = 1:length(sphereGroundTruth)
    % ground-truth data [y, x, z] (stored this way for plotting in imshow) 
    GndTruthData = sphereGroundTruth{step};
    xGndTruth(step) = GndTruthData(2); % ground truth x-component
    yGndTruth(step) = GndTruthData(1); % ground truth y-component
    zGndTruth(step) = GndTruthData(3); % ground truth z-component
end


% select out the cardinal planes
[xSpherePro,ySpherePro,zSpherePro,xGndPro,yGndPro,zGndPro] =...
    selectCardinalPlanes(centerRow,centerCol,centerSlice,...
    xSphere,ySphere,zSphere,xGndTruthFinal,yGndTruthFinal,zGndTruthFinal);

% apply procrustes method.
plotPlaneFitCard = 0;
proSphereBefore = [xSpherePro',ySpherePro',zSpherePro'];
proGround = [xGndPro',yGndPro',zGndPro'];
[d, proSphereAfter,transform] = procrustes(proGround,proSphereBefore,'scaling',false);
% apply the procrustes transform to the entire volume now
% see doc procrustes for more
proT = transform.T; % rotation matrix from the transform
proB = transform.b;
proCtemp = transform.c;
% the translation component is the same for all spheres. Make an array that
% is of the same dimensions as xSphere, ySphere, and zSphere
xProC = proCtemp(1,1).*ones(length(xSphere),1);
yProC = proCtemp(1,2).*ones(length(xSphere),1);
zProC = proCtemp(1,3).*ones(length(xSphere),1);
proC = [xProC,yProC,zProC];
tempTransDim = size(xSphere);
if (tempTransDim(1) > 1)
    transformAll = [xSphere,ySphere,zSphere];
else
    transformAll = [xSphere',ySphere',zSphere'];
end
proSphereAll = (proB*transformAll*proT + proC);


% use this code to verify that the transform was correcly performed
% newLoc = (proB*proSphereBefore*proT + proCtemp);
% totDiff = 0;
% for step = 1:size(newLoc,1)
%    totDiff =  totDiff + sqrt((newLoc(step,1) - proSphereAfter(step,1))^2 + ...
%        (newLoc(step,2) - proSphereAfter(step,2))^2 + ...
%        (newLoc(step,3) - proSphereAfter(step,3))^2);
% end
    
proX = proSphereAll(:,1);
proY = proSphereAll(:,2);
proZ = proSphereAll(:,3);


finalSphereDataX = proX';
finalSphereDataY = proY';
finalSphereDataZ = proZ';
% dummy variables so code will run but doesnt apply any additional shifts
shiftCorrection = [0,0,0];
shiftCorrectionReCalc  = [0,0,0];
shiftForRotTrans = [0, 0, 0];
shiftForRotCoronal = [0, 0, 0];
shiftForRotSagittal = [0, 0, 0];
RotationAngleCoronal = 0;
RotationAngleTrans = 0;
RotationAngleSagittal = 0;
netTransRot = 0;
netCoronalRot = 0;
netSagittalRot = 0;
shiftGndSagittal = 0;
shiftGndCoronal = 0;
shiftGndTransverse = 0;
netShiftSagittal = 0;
netShiftCoronal = 0;
netShiftTransverse = 0;
rotClockTrans = @(x)[1,0,0;0,1,0;0,0,1]; % identity matrix
rotCountTrans = @(x)[1,0,0;0,1,0;0,0,1];
% NOTE: this function does not apply any shifts, just extracts the slice.
% This function was used in an older version of the code to apply
% corrections, but corrections are no longer implemented because of the
% values above.
[sliceShiftTransverseFinal, indAnalysisSphereX,indAnalysisSphereY,indAnalysisSphereZ,...
    indAnalysisGroundX,indAnalysisGroundY,indAnalysisGroundZ] = extractSlice(finalSphereDataX,finalSphereDataY,finalSphereDataZ,...
    xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,centerCol,centerSlice,centerRow,plotPlaneFit,'transverse');

% calcualte the deviation for the center transverse slice to compare to the
% ViewRay analysis software
% Note: if you remove an outlier, this analysis will be redone to ensure
% that if the outlier was in the central slice, it is not included in the
% final result
centerSliceImg = volData(:,:,centerSlice);
[vrCompareDeviationXY,vrCompareDeviationXYZ,vrAvgDevCompViewRayXY,vrAvgDevCompViewRayXYZ,...
    vrPercentPassXY1,vrPercentPassXYZ1,vrPercentPassXY2,vrPercentPassXYZ2,vrPercentPassXY3,vrPercentPassXYZ3,...
    vrAvgDevRadius1XY,vrAvgDevRadius1XYZ,vrAvgDevRadius2XY,vrAvgDevRadius2XYZ,vrAvgDevRadius3XY,vrAvgDevRadius3XYZ,...
    vrMaxDevRadius1XY,vrMaxDevRadius1XYZ,vrMaxDevRadius2XY,vrMaxDevRadius2XYZ,vrMaxDevRadius3XY,vrMaxDevRadius3XYZ,...
    vrAvgDevRadius4XY,vrAvgDevRadius4XYZ,vrPercentPassXY4,vrPercentPassXYZ4,vrMaxDevRadius4XY,vrMaxDevRadius4XYZ]  = ...
    compViewRay(indAnalysisSphereX,indAnalysisSphereY,indAnalysisSphereZ,...
    indAnalysisGroundX,indAnalysisGroundY,indAnalysisGroundZ,voxelHeight,...
    voxelWidth,voxelLength,centerSliceImg,shiftCorrectionReCalc,radiusThreshold1,...
    radiusThreshold2,radiusThreshold3,radiusThreshold4,centerCol,centerSlice,...
    centerRow,radiusSearch1,radiusSearch2,radiusSearch3,radiusSearch4);


% find the number of spheres failing in the different regions specified by
% the user
numFailRadiusSearch1 = 0;
numFailRadiusSearch2 = 0;
numFailRadiusSearch3 = 0;
numFailRadiusSearch4 = 0;
if exist('xSphereFail')
    radiusFail = zeros(1,length(xSphereFail));
    for step = 1:length(xSphereFail)
        radiusFail(step) = sqrt((voxelWidth*(xSphereFail(step) - centerCol))^2 + ...
            (voxelHeight*(ySphereFail(step) - centerRow))^2 + ...
            (voxelLength*(zSphereFail(step) - centerSlice))^2);
        if radiusFail(step) < radiusSearch1
            numFailRadiusSearch1 = numFailRadiusSearch1 + 1;
        end
        if radiusFail(step) < radiusSearch2
            numFailRadiusSearch2 = numFailRadiusSearch2 + 1;
        end
        if radiusFail(step) < radiusSearch3
            numFailRadiusSearch3 = numFailRadiusSearch3 + 1;
        end
        if radiusFail(step) < radiusSearch4
            numFailRadiusSearch4 = numFailRadiusSearch4 + 1;
        end
    end
end


radiusAll = zeros(1,length(finalSphereDataX)); % in 3D, X,Y,Z
for step = 1:length(finalSphereDataX)
    % in 3D
    radiusAll(step) = sqrt((voxelWidth*(finalSphereDataX(step) - centerCol)).^2 +...
        (voxelHeight*(finalSphereDataY(step) - centerRow)).^2 +...
        (voxelLength*(finalSphereDataZ(step) - centerSlice)).^2);
end

% if you know you need to remove an outlier because an incorrect point,
% like one of the rods that passes through the phantom, has been found
% instead of a sphere. A GUI will allow you to cycle through all the
% images and select the spheres you want to omit from the analysis
outliersRemoved = 0; %
lastCalc = 0; % allows the computation of one more calculation following removing the data
while outliersRemoved == 0
    % calculate the distance from isocenter for all the spheres and determine
    % the percent that pass at the two tolerances
    spheres3DThresh1 = 0; % number of spheres within range of first search radius
    spheres3DThresh2 = 0; % number of spheres within range of second search radius
    spheres3DThresh3 = 0; % number of spheres within range of third search radius
    spheres3DThresh4 = 0; % number of spheres within range of fourth search radius
    spheres3DThresh1Pass = 0; % number of spheres passing tolerance within range of first search radius
    spheres3DThresh2Pass = 0; % number of spheres passing tolerance within range of second search radius
    spheres3DThresh3Pass = 0; % number of spheres passing tolerance within range of third search radius
    spheres3DThresh4Pass = 0; % number of spheres passing tolerance within range of fourth search radius
    spheres2DThresh1 = 0; % number of spheres within range of first search distance
    spheres2DThresh2 = 0; % number of spheres within range of second search distance
    spheres2DThresh3 = 0; % number of spheres within range of third search distance
    spheres2DThresh4 = 0; % number of spheres within range of fourth search distance
    spheres2DThresh1Pass = 0; % number of spheres passing tolerance within range of first search distance
    spheres2DThresh2Pass = 0; % number of spheres passing tolerance within range of second search distance
    spheres2DThresh3Pass = 0; % number of spheres passing tolerance within range of third search distance
    spheres2DThresh4Pass = 0; % number of spheres passing tolerance within range of fourth search distance
    % initialize more arrays
    diffXY = zeros(1,length(finalSphereDataX));
    diffX = zeros(1,length(finalSphereDataX));
    diffY = zeros(1,length(finalSphereDataX));
    diffZ = zeros(1,length(finalSphereDataX));
    distZ = zeros(1,length(finalSphereDataX));
    finalDiffCorr = zeros(1,length(finalSphereDataX));
    % radius 1 = radiusSearch1, radius 2 = radiusSearch2
    avgDiff3DRadius1 = 0; % average deviation between ground truth and sphere location
    avgDiff3DRadius2 = 0; % average deviation between ground truth and sphere location
    avgDiff3DRadius3 = 0; % average deviation between ground truth and sphere location
    avgDiff3DRadius4 = 0; % average deviation between ground truth and sphere location
    avgDiff2DRadius1 = 0; % average deviation between ground truth and sphere location
    avgDiff2DRadius2 = 0; % average deviation between ground truth and sphere location
    avgDiff2DRadius3 = 0; % average deviation between ground truth and sphere location
    avgDiff2DRadius4 = 0; % average deviation between ground truth and sphere location
    avgDiff3DAll = 0; % average deviation for all the spheres
    avgDiff2DAll =0; % average deviation for all the spheres
    % deviation specific pass rate
    passAll2D = 0;
    passAll3D = 0;
    diff3DThresh1 = 0;
    diff2DThresh1 = 0;
    diff3DThresh2 = 0;
    diff2DThresh2 = 0;
    diff3DThresh3 = 0;
    diff2DThresh3 = 0;
    diff3DThresh4 = 0;
    diff2DThresh4 = 0;
    avgDiff2DCentralAxis1 = 0;
    avgDiff2DCentralAxis2 = 0;
    avgDiff2DCentralAxis3 = 0;
    avgDiff2DCentralAxis4 = 0;
    countCentAxisRadius1 = 0;
    countCentAxisRadius2 = 0;
    countCentAxisRadius3 = 0;
    countCentAxisRadius4 = 0;
    radiusAllXY = zeros(1,length(radiusAll)); % in 2D, X,Y only
    for step = 1:length(finalSphereDataX)
        radiusAllXY(step) = sqrt((voxelWidth*(finalSphereDataX(step) - centerCol)).^2 +...
            (voxelHeight*(finalSphereDataY(step) - centerRow)).^2);
        diffXY(step) = sqrt((voxelWidth*(finalSphereDataX(step) - xGndTruthFinal(step))).^2 + ...
            (voxelHeight*(finalSphereDataY(step) - yGndTruthFinal(step))).^2);
        finalDiffCorr(step) = sqrt((voxelWidth*(finalSphereDataX(step) - xGndTruthFinal(step))).^2 + ...
            (voxelHeight*(finalSphereDataY(step) - yGndTruthFinal(step))).^2 + ...
            (voxelLength*(finalSphereDataZ(step) - zGndTruthFinal(step))).^2);
        diffX(step) = abs(voxelWidth*(finalSphereDataX(step) - xGndTruthFinal(step)));
        diffY(step) = abs(voxelHeight*(finalSphereDataY(step) - yGndTruthFinal(step)));
        diffZ(step) = abs(voxelLength*(finalSphereDataZ(step) - zGndTruthFinal(step)));
        % used for plotting deviationin x-y as a function of slice location
        distZ(step) = sqrt((voxelLength*(finalSphereDataZ(step) - centerSlice)).^2);
        avgDiff3DAll = avgDiff3DAll + finalDiffCorr(step);
        avgDiff2DAll = avgDiff2DAll + diffXY(step);
        % test different thresholds
        % first threshold
        if radiusAll(step) < radiusSearch1
            spheres3DThresh1 = spheres3DThresh1 + 1;
            avgDiff3DRadius1 = avgDiff3DRadius1 + finalDiffCorr(step);
            avgDiff2DRadius1 = avgDiff2DRadius1 + diffXY(step);
            % 3D case
            if finalDiffCorr(step) < radiusThreshold1
                spheres3DThresh1Pass = spheres3DThresh1Pass + 1;
            end
            % if they are within the distance in 3D they will also be within the
            % distance in 2D case
            spheres2DThresh1 = spheres2DThresh1 + 1;
            if diffXY(step) < radiusThreshold1
                spheres2DThresh1Pass = spheres2DThresh1Pass + 1;
            end
            % store the deviations
            diff3DThresh1(spheres3DThresh1) = finalDiffCorr(step);
            diff2DThresh1(spheres2DThresh1) = diffXY(step);
        end
        % second threshold
        if radiusAll(step) < radiusSearch2
            spheres3DThresh2 = spheres3DThresh2 + 1;
            avgDiff3DRadius2 = avgDiff3DRadius2 + finalDiffCorr(step);
            avgDiff2DRadius2 = avgDiff2DRadius2 + diffXY(step);
            % 3D case
            if finalDiffCorr(step) < radiusThreshold2
                spheres3DThresh2Pass = spheres3DThresh2Pass + 1;
            end
            % if they are within the distance in 3D they will also be within the
            % distance in 2D case
            spheres2DThresh2 = spheres2DThresh2 + 1;
            if diffXY(step) < radiusThreshold2
                spheres2DThresh2Pass = spheres2DThresh2Pass + 1;
            end
            % store the deviations
            diff3DThresh2(spheres3DThresh2) = finalDiffCorr(step);
            diff2DThresh2(spheres2DThresh2) = diffXY(step);
        end
        % third threshold
        if radiusAll(step) < radiusSearch3
            spheres3DThresh3 = spheres3DThresh3 + 1;
            avgDiff3DRadius3 = avgDiff3DRadius3 + finalDiffCorr(step);
            avgDiff2DRadius3 = avgDiff2DRadius3 + diffXY(step);
            % 3D case
            if finalDiffCorr(step) < radiusThreshold3
                spheres3DThresh3Pass = spheres3DThresh3Pass + 1;
            end
            % if they are within the distance in 3D they will also be within the
            % distance in 2D case
            spheres2DThresh3 = spheres2DThresh3 + 1;
            if diffXY(step) < radiusThreshold3
                spheres2DThresh3Pass = spheres2DThresh3Pass + 1;
            end
            % store the deviations
            diff3DThresh3(spheres3DThresh3) = finalDiffCorr(step);
            diff2DThresh3(spheres2DThresh3) = diffXY(step);
        end
        % fourth threshold
        if radiusAll(step) < radiusSearch4
            spheres3DThresh4 = spheres3DThresh4 + 1;
            avgDiff3DRadius4 = avgDiff3DRadius4 + finalDiffCorr(step);
            avgDiff2DRadius4 = avgDiff2DRadius4 + diffXY(step);
            % 3D case
            if finalDiffCorr(step) < radiusThreshold4
                spheres3DThresh4Pass = spheres3DThresh4Pass + 1;
            end
            % if they are within the distance in 3D they will also be within the
            % distance in 2D case
            spheres2DThresh4 = spheres2DThresh4 + 1;
            if diffXY(step) < radiusThreshold4
                spheres2DThresh4Pass = spheres2DThresh4Pass + 1;
            end
            % store the deviations
            diff3DThresh4(spheres3DThresh4) = finalDiffCorr(step);
            diff2DThresh4(spheres2DThresh4) = diffXY(step);
        end
        % determine if deviation thresholds have been exceeded
        % 2 mm deviation
        if diffXY(step)<= allThreshold
            passAll2D = passAll2D + 1;
        end
        % 3 mm deviation
        if finalDiffCorr(step)<= allThreshold
            passAll3D = passAll3D + 1;
        end
        % Distance from central axis calcuation
        if radiusAllXY(step) < radiusSearch1
            % average deviation in 2D for spheres within radius1 of central axis
            countCentAxisRadius1 = countCentAxisRadius1 + 1;
            avgDiff2DCentralAxis1 = avgDiff2DCentralAxis1 + diffXY(step);
        end
        if radiusAllXY(step) < radiusSearch2
            % average deviation in 2D for spheres within radius2 of central axis
            countCentAxisRadius2 = countCentAxisRadius2 + 1;
            avgDiff2DCentralAxis2 = avgDiff2DCentralAxis2 + diffXY(step);
        end
        if radiusAllXY(step) < radiusSearch3
            % average deviation in 2D for spheres within radius2 of central axis
            countCentAxisRadius3 = countCentAxisRadius3 + 1;
            avgDiff2DCentralAxis3 = avgDiff2DCentralAxis3 + diffXY(step);
        end
        if radiusAllXY(step) < radiusSearch4
            % average deviation in 2D for spheres within radius2 of central axis
            countCentAxisRadius4 = countCentAxisRadius4 + 1;
            avgDiff2DCentralAxis4 = avgDiff2DCentralAxis4 + diffXY(step);
        end
    end
    
    % the percentage that pass
    pass3DThresh1 = 100*spheres3DThresh1Pass/spheres3DThresh1;
    pass3DThresh2 = 100*spheres3DThresh2Pass/spheres3DThresh2;
    pass3DThresh3 = 100*spheres3DThresh3Pass/spheres3DThresh3;
    pass3DThresh4 = 100*spheres3DThresh4Pass/spheres3DThresh4;
    pass2DThresh1 = 100*spheres2DThresh1Pass/spheres2DThresh1;
    pass2DThresh2 = 100*spheres2DThresh2Pass/spheres2DThresh2;
    pass2DThresh3 = 100*spheres2DThresh3Pass/spheres2DThresh3;
    pass2DThresh4 = 100*spheres2DThresh4Pass/spheres2DThresh4;
    % average deviation between ground truth and sphere location
    avgDiff3DAll = avgDiff3DAll/length(finalSphereDataX);
    avgDiff2DAll = avgDiff2DAll/length(finalSphereDataX);
    avgDiff3DRadius1 = avgDiff3DRadius1/spheres3DThresh1;
    avgDiff3DRadius2 = avgDiff3DRadius2/spheres3DThresh2;
    avgDiff3DRadius3 = avgDiff3DRadius3/spheres3DThresh3;
    avgDiff3DRadius4 = avgDiff3DRadius4/spheres3DThresh4;
    avgDiff2DRadius1 = avgDiff2DRadius1/spheres2DThresh1;
    avgDiff2DRadius2 = avgDiff2DRadius2/spheres2DThresh2;
    avgDiff2DRadius3 = avgDiff2DRadius3/spheres2DThresh3;
    avgDiff2DRadius4 = avgDiff2DRadius4/spheres2DThresh4;
    avgDiff2DCentralAxis1 = avgDiff2DCentralAxis1/countCentAxisRadius1;
    avgDiff2DCentralAxis2 = avgDiff2DCentralAxis2/countCentAxisRadius2;
    avgDiff2DCentralAxis3 = avgDiff2DCentralAxis3/countCentAxisRadius3;
    avgDiff2DCentralAxis4 = avgDiff2DCentralAxis4/countCentAxisRadius4;
    % the percentage that were found (passed correlation coefficient tolerance)
    spheresFoundThresh1 = 100*spheres3DThresh1/(numFailRadiusSearch1 + spheres3DThresh1);
    spheresFoundThresh2 = 100*spheres3DThresh2/(numFailRadiusSearch2 + spheres3DThresh2);
    spheresFoundThresh3 = 100*spheres3DThresh3/(numFailRadiusSearch3 + spheres3DThresh3);
    spheresFoundThresh4 = 100*spheres3DThresh4/(numFailRadiusSearch4 + spheres3DThresh4);
    
    % calculate the maximum deviation
    finalMaxDeviation2D = max(diffXY);
    finalMaxDeviation3D = max(finalDiffCorr);
    maxDiff3DThresh1 = max(diff3DThresh1);
    maxDiff2DThresh1 = max(diff2DThresh1);
    maxDiff3DThresh2 = max(diff3DThresh2);
    maxDiff2DThresh2 = max(diff2DThresh2);
    maxDiff3DThresh3 = max(diff3DThresh3);
    maxDiff2DThresh3 = max(diff2DThresh3);
    maxDiff3DThresh4 = max(diff3DThresh4);
    maxDiff2DThresh4 = max(diff2DThresh4);
    
    passAll2DPercent = 100*passAll2D/length(finalSphereDataX);
    passAll3DPercent = 100*passAll3D/length(finalSphereDataX);    
    % NOTE: this if statement is only evaluated if the user decides to remove outliers.
    if lastCalc == 1; % this is the last iteration of calculations
        outliersRemoved = 1;
        [sliceShiftTransverse, indAnalysisSphereX,indAnalysisSphereY,indAnalysisSphereZ,...
            indAnalysisGroundX,indAnalysisGroundY,indAnalysisGroundZ] = extractSlice(finalSphereDataX,finalSphereDataY,finalSphereDataZ,...
            xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,centerCol,centerSlice,centerRow,0,'transverse');
        shiftCorrectionReCalc = [0,0,0]; % do not apply any corrections
        
        % calcualte the deviation for the center transverse slice to compare to the
        % ViewRay analysis software now that data has been removed
        centerSliceImg = volData(:,:,centerSlice);
        [vrCompareDeviationXY,vrCompareDeviationXYZ,vrAvgDevCompViewRayXY,vrAvgDevCompViewRayXYZ,...
            vrPercentPassXY1,vrPercentPassXYZ1,vrPercentPassXY2,vrPercentPassXYZ2,vrPercentPassXY3,vrPercentPassXYZ3,...
            vrAvgDevRadius1XY,vrAvgDevRadius1XYZ,vrAvgDevRadius2XY,vrAvgDevRadius2XYZ,vrAvgDevRadius3XY,vrAvgDevRadius3XYZ,...
            vrMaxDevRadius1XY,vrMaxDevRadius1XYZ,vrMaxDevRadius2XY,vrMaxDevRadius2XYZ,vrMaxDevRadius3XY,vrMaxDevRadius3XYZ,...
            vrAvgDevRadius4XY,vrAvgDevRadius4XYZ,vrPercentPassXY4,vrPercentPassXYZ4,vrMaxDevRadius4XY,vrMaxDevRadius4XYZ]  = ...
            compViewRay(indAnalysisSphereX,indAnalysisSphereY,indAnalysisSphereZ,...
            indAnalysisGroundX,indAnalysisGroundY,indAnalysisGroundZ,voxelHeight,...
            voxelWidth,voxelLength,centerSliceImg,shiftCorrectionReCalc,radiusThreshold1,...
            radiusThreshold2,radiusThreshold3,radiusThreshold4,centerCol,centerSlice,...
            centerRow,radiusSearch1,radiusSearch2,radiusSearch3,radiusSearch4);
    end
    if outliersRemoved == 0;
        figure
        plot(radiusAll,finalDiffCorr,'b.');
        xlabel('Distance from Isocenter (mm)','FontSize',22)
        ylabel('Deviation (mm)','FontSize',22)
        title('Entire Volume 3D Deviation from Ground Truth','FontSize',22)
        if(input('Do you need to check for spheres in artifact regions? (y = 1, n = 0):'))
            disp('Select any spheres that should be removed from the analysis')
            [numRemoved,finalSphereDataX,finalSphereDataY,finalSphereDataZ,...
                xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,radiusAll,finalDiffCorr,...
                diffX,diffY,diffZ,xSphere,ySphere,zSphere] = ...
                removeSpheresWithGUI(finalSphereDataX,finalSphereDataY,finalSphereDataZ,diffX,diffY,diffZ,...
                xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,radiusAll,finalDiffCorr,...
                volData,xSphere,ySphere,zSphere);
            figure
            plot(radiusAll,finalDiffCorr,'b.');
            xlabel('Distance from Isocenter (mm)','FontSize',22)
            ylabel('Deviation (mm)','FontSize',22)
            lastCalc = input('Are all the spheres in artifact regions removed? (y = 1, n = 0): ');
            if (lastCalc ~=1)
                radiusRemove(1) = input('Input min radius from isocenter to be highlighted:');
                radiusRemove(2) = input('Input max radius from isocenter to be highlighted:');
                rangeDevRemove(1) = input('Input min deviation to be highlighted:');
                rangeDevRemove(2) = input('Input max deviation to be highlighted:');
                plotOutliers(finalSphereDataX,finalSphereDataY,finalSphereDataZ,diffX,diffY,diffZ,...
                    xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,radiusAll,finalDiffCorr,...
                    radiusThreshold1,volData,sphRangeX,sphRangeY,sphRangeZ,sphereVol,radiusRemove,rangeDevRemove,0)
            end
            % if some spheres were removed, re-calculate procrustes
            if numRemoved > 0
                [xSpherePro,ySpherePro,zSpherePro,xGndPro,yGndPro,zGndPro] =...
                    selectCardinalPlanes(centerRow,centerCol,centerSlice,...
                    finalSphereDataX,finalSphereDataY,finalSphereDataZ,xGndTruthFinal,yGndTruthFinal,zGndTruthFinal);
                proSphereBefore = [xSpherePro',ySpherePro',zSpherePro'];
                proGround = [xGndPro',yGndPro',zGndPro'];
                [d, proSphereAfter,transform] = procrustes(proGround,proSphereBefore,'scaling',false);
                % apply the procrustes transform to the entire volume now
                % see doc procrustes for more
                proT = transform.T; % rotation matrix from the transform
                proB = transform.b;
                proCtemp = transform.c;
                xProC = proCtemp(1,1).*ones(length(finalSphereDataX),1);
                yProC = proCtemp(1,2).*ones(length(finalSphereDataX),1);
                zProC = proCtemp(1,3).*ones(length(finalSphereDataX),1);
                proC = [xProC,yProC,zProC];
                tempTransDim = size(finalSphereDataX);
                if (tempTransDim(1) > 1)
                    transformAll = [finalSphereDataX,finalSphereDataY,finalSphereDataZ];
                else
                    transformAll = [finalSphereDataX',finalSphereDataY',finalSphereDataZ'];
                end
                proSphereAll = (proB*transformAll*proT + proC);
                proX = proSphereAll(:,1);
                proY = proSphereAll(:,2);
                proZ = proSphereAll(:,3);
                finalSphereDataX = proX';
                finalSphereDataY = proY';
                finalSphereDataZ = proZ';
            end
        else
            % we don't want to remove any outliers
            outliersRemoved = 1;
            lastCalc = 1;
        end
    end
end

% the total number of spheres found after
spheresFoundAll = 100*length(radiusAll)/nTotData;




% calculate the combined uncertainty in finding the location of the sphere
% and the manufacturing of the sphere location
[devFromExpected,avgDevFromExpected,numSpheresCalDevExp, calcSurLoc,foundSurSphere,rmsDevFromExp,rmsDevFromGnd] =...
    calcUncertainty(centerRow,centerCol,centerSlice,voxelHeight,voxelWidth,voxelLength,finalSphereDataX,...
    finalSphereDataY,finalSphereDataZ,xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,phantomSim,xSpacing,ySpacing,zSpacing);

% compute error in a single measurement, current metric only measures error
% in one dimension (multiply by sqrt(3)), and contains the error for two
% spheres (divide by sqrt(2))
rmsDevFromExp = rmsDevFromExp*sqrt(3)/sqrt(2);
rmsDevFromGnd = rmsDevFromGnd*sqrt(3)/sqrt(2);
% calculate the SNR in the image
[SNR, avgSig, bgStandardDev] = calcSNR(volData,voxelHeight,voxelWidth,voxelLength,centerSlice,...
centerCol,centerRow,sphereDiameter,xGndTruthFinal,yGndTruthFinal,zGndTruthFinal,...
xGndTruthFinal,yGndTruthFinal,zGndTruthFinal);

% Here we plot the center slice of the ViewRay analysis and save the image
% if necessary.
%
% For plotting we want to use the original locations of the spheres found
% by the software, along with the ground truth locations rotated and
% shifted by the opposite correction applied to the found locations. This
% makes the plots of the found and ground truth locations agree better with
% the images.
% 
% using a single net correction
saveDeformationImg = 0;
radiusThreshold3 = 1;
[volPass2D,volFail2D,volPass3D,volFail3D,volGndTruth,...
    volAnalysisRegions] = ...
    correctGndTruth_Procrustes(proT,proB,proC,centerRow,centerCol,centerSlice,xSphere,ySphere,zSphere,...
    finalSphereDataX,finalSphereDataY,finalSphereDataZ,xGndTruthFinal,yGndTruthFinal,...
    zGndTruthFinal,voxelWidth,voxelHeight,voxelLength,...
    radiusSearch1,radiusSearch2,radiusSearch3,radiusThreshold1,...
    radiusThreshold2,radiusThreshold3,centerSliceImg,volData,saveDeformationImg);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('Phantom3D_ImageOnly_0deg172s','-dtiff','-r300')

viewRayPlotTitle = 'CompareViewRay';
if saveImages == 1
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(viewRayPlotTitle,'-dtiff','-r300')
end

%     fig = gcf;
%     fig.PaperPositionMode = 'auto';
%     print('DeformationFieldMap','-dtiff','-r300')
disp(' ') % add some spaces 
disp(' ') % add some spaces  
disp(' ') % add some spaces 
disp(strcat(['Deviation from expected distance RMS for the ',num2str(numSpheresCalDevExp),' spheres surrounding isocenter: ',num2str(rmsDevFromExp),' (mm)']))
disp(strcat(['Deviation from ground truth RMS for the ',num2str(numSpheresCalDevExp+1),' spheres nearest to isocenter: ',num2str(rmsDevFromGnd),' (mm)']))
disp(strcat(['Spheres found with corr. threshold of ',num2str(corrThreshold),' in a ',num2str(radiusSearch1),' (mm) radius: ',num2str(spheresFoundThresh1),'%'])) 
disp(strcat(['Spheres found with corr. threshold of ',num2str(corrThreshold),' in a ',num2str(radiusSearch2),' (mm) radius: ',num2str(spheresFoundThresh2),'%'])) 
disp(strcat(['Spheres found with corr. threshold of ',num2str(corrThreshold),' in a ',num2str(radiusSearch3),' (mm) radius: ',num2str(spheresFoundThresh3),'%'])) 
disp(strcat(['Spheres found  with corr. threshold of ',num2str(corrThreshold),' in entire volume: ',num2str(spheresFoundAll),'%'])) 
disp(strcat([num2str(radiusSearch1),' (mm) central axis dist, avg. 2D deviation: ', num2str(avgDiff2DCentralAxis1), ' (mm)']))
disp(strcat([num2str(radiusSearch2),' (mm) central axis dist, avg. 2D deviation: ', num2str(avgDiff2DCentralAxis2), ' (mm)']))
disp(strcat([num2str(radiusSearch3),' (mm) central axis dist, avg. 2D deviation: ', num2str(avgDiff2DCentralAxis3), ' (mm)']))
disp(strcat([num2str(radiusSearch1),' (mm) radius volume avg. 3D deviation: ', num2str(avgDiff3DRadius1), ' (mm)']))
disp(strcat([num2str(radiusSearch1),' (mm) radius volume avg. 2D deviation: ', num2str(avgDiff2DRadius1), ' (mm)']))
disp(strcat([num2str(radiusSearch2),' (mm) radius volume avg. 3D deviation: ', num2str(avgDiff3DRadius2), ' (mm)']))
disp(strcat([num2str(radiusSearch2),' (mm) radius volume avg. 2D deviation: ', num2str(avgDiff2DRadius2), ' (mm)']))
disp(strcat([num2str(radiusSearch3),' (mm) radius volume avg. 3D deviation: ', num2str(avgDiff3DRadius3), ' (mm)']))
disp(strcat([num2str(radiusSearch3),' (mm) radius volume avg. 2D deviation: ', num2str(avgDiff2DRadius3), ' (mm)']))
disp(strcat(['Whole volume avg. 3D deviation: ', num2str(avgDiff3DAll), ' (mm)']))
disp(strcat(['Whole volume avg. 2D deviation: ', num2str(avgDiff2DAll), ' (mm)']))
disp(strcat(['3D deviation, percent passing in ',num2str(radiusSearch1),' (mm) radius volume, tolerance of ',num2str(radiusThreshold1),' (mm): ', num2str(pass3DThresh1),'%']))
disp(strcat(['2D deviation, percent passing in ',num2str(radiusSearch1),' (mm) radius volume, tolerance of ',num2str(radiusThreshold1),' (mm): ', num2str(pass2DThresh1),'%']))
disp(strcat(['3D deviation, percent passing in ',num2str(radiusSearch2),' (mm) radius volume, tolerance of ',num2str(radiusThreshold2),' (mm): ', num2str(pass3DThresh2),'%']))
disp(strcat(['2D deviation, percent passing in ',num2str(radiusSearch2),' (mm) radius volume, tolerance of ',num2str(radiusThreshold2),' (mm): ', num2str(pass2DThresh2),'%']))
disp(strcat(['3D deviation, percent passing in ',num2str(radiusSearch3),' (mm) radius volume, tolerance of ',num2str(radiusThreshold3),' (mm): ', num2str(pass3DThresh3),'%']))
disp(strcat(['2D deviation, percent passing in ',num2str(radiusSearch3),' (mm) radius volume, tolerance of ',num2str(radiusThreshold3),' (mm): ', num2str(pass2DThresh3),'%']))
disp(strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch1),' (mm) dev <= ',...
    num2str(radiusThreshold1),' (mm), all 3D: ',num2str(vrPercentPassXYZ1),'%']))
disp(strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch1),' (mm) dev <= ',...
    num2str(radiusThreshold1),' (mm), just x,y: ',num2str(vrPercentPassXY1),'%']))
disp(strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch1),' (mm)'...
    ' avg 2D deviation: ',num2str(vrAvgDevRadius1XY),' (mm)']))
disp(strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch2),' (mm) dev <= ',...
    num2str(radiusThreshold2),' (mm), all 3D: ',num2str(vrPercentPassXYZ2),'%']))
disp(strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch2),' (mm) dev <= ',...
    num2str(radiusThreshold2),' (mm), just x,y: ',num2str(vrPercentPassXY2),'%']))
disp(strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch2),' (mm)'...
    ' avg 2D deviation: ',num2str(vrAvgDevRadius2XY),' (mm)']))
disp(strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch3),' (mm) dev <= ',...
    num2str(radiusThreshold3),' (mm), all 3D: ',num2str(vrPercentPassXYZ3),'%']))
disp(strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch3),' (mm) dev <= ',...
    num2str(radiusThreshold3),' (mm), just x,y: ',num2str(vrPercentPassXY3),'%']))
disp(strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch3),' (mm)'...
    ' avg 2D deviation: ',num2str(vrAvgDevRadius3XY),' (mm)']))
disp(strcat(['Percent of spheres in entire volume with less than ',num2str(allThreshold),' mm 2D deviation: ',num2str(passAll2DPercent),'%'])) 
disp(strcat(['Percent of spheres in entire volume with less than ',num2str(allThreshold),' mm 3D deviation: ',num2str(passAll3DPercent),'%'])) 
disp(strcat(['The maximum 3D deviation of a sphere in the entire volume: ',num2str(finalMaxDeviation3D),' (mm)']))
disp(strcat(['The maximum 2D deviation of a sphere in the entire volume: ',num2str(finalMaxDeviation2D),' (mm)']))
disp(strcat(['The maximum 3D deviation in ',num2str(radiusSearch1),' (mm) radius volume: ',num2str(maxDiff3DThresh1),' (mm)']))
disp(strcat(['The maximum 2D deviation in ',num2str(radiusSearch1),' (mm) radius volume: ',num2str(maxDiff2DThresh1),' (mm)']))
disp(strcat(['The maximum 3D deviation in ',num2str(radiusSearch2),' (mm) radius volume: ',num2str(maxDiff3DThresh2),' (mm)']))
disp(strcat(['The maximum 2D deviation in ',num2str(radiusSearch2),' (mm) radius volume: ',num2str(maxDiff2DThresh2),' (mm)']))
disp(strcat(['The maximum 3D deviation in ',num2str(radiusSearch3),' (mm) radius volume: ',num2str(maxDiff3DThresh3),' (mm)']))
disp(strcat(['The maximum 2D deviation in ',num2str(radiusSearch3),' (mm) radius volume: ',num2str(maxDiff2DThresh3),' (mm)']))
disp(strcat(['Center Sphere SNR: ',num2str(SNR),', Average Signal: ',num2str(avgSig),', Bg. Standard Dev: ',num2str(bgStandardDev)]))
% Plot the deviation in just two dimensions as a function of their distance
% from the y-axis
% plotDeviation(radiusAllXY,diffXY,binSize);
% xlabel('Distance from Y-Axis (mm)','FontSize',22)
% ylabel('XZ-Deviation (mm)','FontSize',22)
% title('Entire Volume 2D Deviation from Ground Truth','FontSize',22)
%
% Plot the deviation in three dimensions as a function of their distance
% from ioscenter
plotDeviation(radiusAll,finalDiffCorr,binSize);
xlabel('Distance from isocenter (mm)','FontSize',22)
ylabel('3D Deviation (mm)','FontSize',22)
% title('Entire Volume 3D Deviation from Ground Truth','FontSize',22)
if saveImages == 1
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print('DeviationPlot','-dtiff','-r300')
end
% Plot the deviation in two dimensions as a function of their distance
% from the central axis
plotDeviation(radiusAllXY,diffXY,binSize)
xlabel('Distance from central axis (mm)','FontSize',22)
ylabel('2D Deviation (mm)','FontSize',22)
% title('Entire Volume 2D Deviation from Ground Truth','FontSize',22)
if saveImages == 1
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print('CentralAxisDeviationPlot','-dtiff','-r300')
end


% save the data to an excel file
if saveExcelFiles == 1
    saveExcel(gantryAngle,dataDir,numSpheresCalDevExp,rmsDevFromExp,...
    rmsDevFromGnd,corrThreshold,spheresFoundAll,radiusSearch1,spheresFoundThresh1,radiusSearch2,...
    spheresFoundThresh2,radiusSearch3,spheresFoundThresh3,avgDiff2DCentralAxis1,avgDiff2DCentralAxis2,avgDiff2DCentralAxis3,avgDiff3DRadius1,...
    avgDiff2DRadius1,avgDiff3DRadius2,avgDiff2DRadius2,avgDiff3DRadius3,avgDiff2DRadius3,avgDiff3DAll,avgDiff2DAll,...
    radiusThreshold1,radiusThreshold2,radiusThreshold3,pass3DThresh1,pass2DThresh1,pass3DThresh2,pass2DThresh2,pass3DThresh3,pass2DThresh3,...
    vrPercentPassXYZ1,vrPercentPassXY1,vrAvgDevRadius1XY,vrPercentPassXYZ2,vrPercentPassXY2,...
    vrAvgDevRadius2XY,vrPercentPassXYZ3,vrPercentPassXY3,vrAvgDevRadius3XY,passAll2DPercent,passAll3DPercent,finalMaxDeviation3D,finalMaxDeviation2D,...
    maxDiff3DThresh1,maxDiff2DThresh1,maxDiff3DThresh2,maxDiff2DThresh2,maxDiff3DThresh3,maxDiff2DThresh3,SNR,avgSig,...
    bgStandardDev,radiusSearch4,radiusThreshold4,spheresFoundThresh4,avgDiff2DCentralAxis4,...
    avgDiff3DRadius4,avgDiff2DRadius4,pass3DThresh4,pass2DThresh4,vrPercentPassXYZ4,...
    vrPercentPassXY4,vrAvgDevRadius4XY,maxDiff3DThresh4,maxDiff2DThresh4,allThreshold);
end


%% Plot the final sphere locations with the image
% detrmine opposite rotation matrix
guiData{1} = volData; % phantom image data
guiData{2} = volAnalysisRegions; % analysis regions
guiData{3} = volPass3D; % markers of sphere locations for fusing to the image
guiData{4} = volFail3D; % markers of sphere locations for fusing to the image
guiData{5} = volGndTruth; % ground truth data
thisGUI = CheckCorrGUI2(guiData);