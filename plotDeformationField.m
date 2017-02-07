% Function to plot the deformation field map for a slice in the phantom.
% NOTE: the indices specified below MUST correspond to the locations in the
% image. Additinoally, the ground truth locations must not be corrected
% by a transform
%
% Input:
% xSph (index) The sphere found x-location
% ySph (index) The sphere found y-location
% zSph (index) The sphere found z-location
% xGnd (index) The ground truth x-location
% yGnd (index) The ground truth y-location
% zGnd (index) The ground truth z-location
% deviation (mm) The deviation between the sphere and ground truth locations
% sliceImg The image of the slice for these spheres 
% voxelWidth (mm) The width of the voxels
% voxelHeight (mm) The height of the voxels
% voxelLength (mm) The length of the voxels
% scaleFactor The factor by which you want to upscale the image
% saveDeformationImg (y=1,n=0) Whether or not to save the deformation field image
%
% Output:
%
% John Ginn
% Created: 1/25/17
% Modified: 1/26/17

function [] = plotDeformationField(xSph,ySph,zSph,xGnd,yGnd,zGnd,deviation,sliceImg,...
    voxelWidth,voxelHeight,voxelLength,scaleFactor,saveDeformationImg)
% sets the scale for the colorbar plot
scaleMin = 0;
scaleMax = 3;
extraDist = 0; % 3 mm extra plotting

% find the maximum number of spheres in each row
spheresInRow = [29 29 29 29 29 29 29 29 27 27 27 27 25 25 23 23 21 19 17 13 7]; 
maxSpheresDim = max(spheresInRow) + 2;
% extra distance for plotting
% find bounds of region for deformation map
minX = round(min(xSph) - extraDist/voxelWidth*scaleFactor);
maxX = round(max(xSph) + extraDist/voxelWidth*scaleFactor);
minY = round(min(ySph) - extraDist/voxelHeight*scaleFactor);
maxY = round(max(ySph) + extraDist/voxelHeight*scaleFactor);
minZ = round(min(zSph) - extraDist/voxelLength*scaleFactor);
maxZ = round(max(zSph) + extraDist/voxelLength*scaleFactor);

xArray = minX:maxX;
yArray = minY:maxY;
zArray = minZ:maxZ;
[meshX,meshY] = meshgrid(xArray,yArray);

% make data monotonically increase
% [xSphSort, xSortInd] = sort(xSph);
% % sort y and z accordingly
% yTemp = ySph; % don't overwrite the data you need to access!
% zTemp = zSph; % don't overwrite the data you need to access!
% for stepOther = 1:length(xSortInd)
%     ySph(stepOther) = yTemp(xSortInd(stepOther));
%     zSph(stepOther) = zTemp(xSortInd(stepOther));
% end
minGndX = (min(xGnd));
maxGndX = (max(xGnd));
minGndY = (min(xGnd));
maxGndY = (max(ySph));
minGndZ = (min(zSph));
maxGndZ = (max(zSph));

xGndUnique = unique(xGnd);
yGndUnique = unique(yGnd);
zGndUnique = unique(zGnd);

[meshGndX,meshGndY] = meshgrid(xGndUnique,yGndUnique);
devMesh = zeros(length(meshGndX(:,1)),length(meshGndX(1,:)));
for step = 1:length(xSph)
    % find the location of this value in the grid
    thisX = find(xGnd(step) == xGndUnique);
    thisY = find(yGnd(step) == yGndUnique);
    devMesh(thisY,thisX) = deviation(step);
end
X = meshGndX;
Y = meshGndY;
V = devMesh;
Xq = meshX;
Yq= meshY;
interpVal = interp2(X,Y,V,Xq,Yq,'linear',0);
% remove any deformation less than 0
for stepX = 1:length(interpVal(1,:))
    for stepY = 1:length(interpVal(:,1))
        if interpVal(stepY,stepX) < scaleMin;
            interpVal(stepY,stepX) = 0;
        end
    end
end
% create the interpolated image. interpVal is not necessarily the same dimension as
% the image
interpDefImg = zeros(length(sliceImg(:,1)),length(sliceImg(1,:)));
interpDefImg(yArray,xArray) = interpVal;
% create the mask for this image
interpMask = zeros(length(sliceImg(:,1)),length(sliceImg(1,:)));
for stepX = 1:length(interpDefImg(1,:))
    for stepY = 1:length(interpDefImg(:,1))
        if interpDefImg(stepY,stepX) > scaleMin;
            interpMask(stepY,stepX) = 0.5;
        end
    end
end
% normalize the image scale of the phantom to the deformation map
maxDef = max(max(interpDefImg));
minDef = 0;
% remove any NaN values from the image and set to zero (around the
% boundaries)
for stepX = 1:length(sliceImg(1,:))
    for stepY = 1:length(sliceImg(:,1))
        TF = isnan(sliceImg(stepY,stepX));
        if TF == 1
           sliceImg(stepY,stepX) = 0; 
        end
    end
end
normImg = sliceImg;
minVal = min(min(normImg));
% make sure the smallest value is at least zero
if minVal < 0
    normImg(:) = normImg(:) + abs(minVal);
end
% normalize the image scale to the deformation scale
maxImg = max(max(normImg));
normFact = maxDef/maxImg;
normImg = normImg.*normFact; 

figure;
h = imagesc(normImg);
colormap('gray')
% Now make an RGB image that matches display from IMAGESC:
C = colormap;  % Get the figure's colormap.
L = size(C,1);
% Scale the matrix to the range of the map.
Gs = round(interp1(linspace(min(normImg(:)),max(normImg(:)),L),1:L,normImg));
colorPhantomImg = reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.
close gcf

% here is the plot that actually gets saved
figure;
imshow(colorPhantomImg,[]);
hold on
h = imshow(interpDefImg,[scaleMin scaleMax]);
hold off
colormap('jet')
barObj = colorbar('eastoutside');
% alpha 0.3
set(h, 'AlphaData', interpMask);
% change the colorbar ticks to include mm
currTicks = barObj.Ticks;
newTicks = cell(1,length(currTicks));
for step = 1:length(currTicks)
    newTicks{step} = strcat([num2str(currTicks(step)),' (mm)']);
end
newTicks{step} = strcat('>',[num2str(currTicks(step)),' (mm)']);
newTicks{1} = strcat([num2str(currTicks(1)),' (mm)']);
colorbar('Ticks',currTicks,...
         'TickLabels',newTicks)
% % This works with the exception of the colorbar
% figure;
% imshow(colorInterpImg);
% hold on
% h = imshow(normImg,[]);
% hold off
% % alpha 0.3
% set(h, 'AlphaData', interpMask);
% set(h,'colorbar',barObjTemp)
% barObj = colorbar('eastoutside');

% custom axes to show distances (the locations of the axes)
numOfTicks = 10;
xAxisSpacing = floor(length(sliceImg(1,:))/numOfTicks);
yAxisSpacing = floor(length(sliceImg(:,1))/numOfTicks);

% xAxisLocation = linspace(1,length(SigDataUpscale(1,:)),numOfTicks);
% yAxisLocation = linspace(1,length(SigDataUpscale(:,1)),numOfTicks);
xAxisLocation = 1:xAxisSpacing:(xAxisSpacing*numOfTicks+1);
yAxisLocation = 1:yAxisSpacing:(yAxisSpacing*numOfTicks+1);
% the values on the axes
xAxisValue = voxelWidth/scaleFactor.*(xAxisLocation - round(median(xAxisLocation)));
yAxisValue = voxelHeight/scaleFactor.*(yAxisLocation - round(median(yAxisLocation)));
yAxisLocation(length(yAxisLocation)) = length(sliceImg(:,1)); % special case for plotting for paper
axis on
axisHandle = gca;
set(gca,'XTickMode','manual')
set(gca,'YTickMode','manual')
set(gca,'XTick',xAxisLocation)
set(gca,'YTick',yAxisLocation)
set(gca,'XTickLabel',xAxisValue)
set(gca,'YTickLabel',yAxisValue)
set(gca,'FontSize',16)
xlabel('x-position (mm)','FontSize',20)
ylabel('z-position (mm)','FontSize',20)

if saveDeformationImg
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print('DeformationFieldMap','-dtiff','-r300')
end
end