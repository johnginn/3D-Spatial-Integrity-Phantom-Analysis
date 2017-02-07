% Function to create a cirlce for overlaying on top of the current image
%
% Input:
% imageData The image of the phantom
% radiusSearch (mm) The distance from isocenter defining this region
% voxelHeight (mm) The height of each voxel in imageData
% voxelWidth (mm) The width of each voxel in imageData
% centerRow The row index of isocenter
% centerCol The column index of isocenter
% sliceShift (mm) The shift of the phantom away from isocenter
%
% Output:
% circleImg The image of the circle for plotting (ones and zeros)
%
% John Ginn
% Created: 11/4/16
% Modified: 11/4/16
function [circleImg] = imageCircle(imageData,radiusSearch,voxelHeight,voxelWidth,...
    centerRow,centerCol,sliceShift)

% number of pixels in the image
nPixels = length(imageData(:,1))*length(imageData(1,:));
% initialize circle image
circleImg = zeros(length(imageData(:,1)),length(imageData(1,:)));

angleStep = 360/nPixels; % for polar coordinates

% take the absolute value of the slice shift
sliceShift = abs(sliceShift);
if sliceShift >= radiusSearch
    circleRadius = 0;
else
    circleRadius = sqrt(radiusSearch^2 - sliceShift^2); % account for slice shift
end
% step through nPixel times to ensure the circle is continuous (no breaks
% in the line)
xDim = length(imageData(1,:));
yDim = length(imageData(:,1));
for step = 1:nPixels;
    % if the circle radius is zero don't plot anything
    if circleRadius > 0
        xPos = round(circleRadius*cosd(angleStep*step)/voxelWidth + centerCol);
        yPos = round(circleRadius*sind(angleStep*step)/voxelHeight + centerRow);
        % make sure radius does not exceeed image
        if (xPos > 0)&&(xPos <=xDim)&&(yPos > 0)&&(yPos <=yDim)
            circleImg(yPos,xPos) = 1; % store a 1 because this is the boundary of the circle
        end
    end
end
%% For debugging
% figure;
% imshow(circleImg,[])
% title('The circle boundary','FontSize',16)

end