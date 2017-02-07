% A function to invert the contrast of the CT images to make them work with
% the weighted sum method of finding the spheres
%
% Input:
% volData The volume data from the CT images
%
% Output:
% invertedVolData The volume data from the CT images after inverting the contrast
%
% John Ginn
% Created: 10/24/16
% Modified: 10/24/16

function [invertedVolData] = invertCT(volData)

xDim = length(volData(1,:,1));
yDim = length(volData(:,1,1));
zDim = length(volData(1,1,:));

dataArray = zeros(1,xDim*yDim*zDim);
count = 0;
for x = 1:xDim
    for y = 1:yDim
        for z = 1:zDim
            count = count + 1;
            % store the data as a single array
            dataArray(count) = volData(y,x,z);
        end
    end
end

% make a histogram to determine threshold for different contrast regions
nBins = 1000; % number of bins in the histogram
figure
hist(dataArray,nBins);
xlabel('attenuation level','FontSize',16)
ylabel('number of voxels','FontSize',16)
title('Voxel Attenuation Histogram','FontSize',16)

% select the plastic signal peak
disp('Select window for attenuation values. (all outside regions will be set to zero)')
disp('Select minimum window position')
[xMin,yMin] = ginput(1);
disp('Select maximum window position')
[xMax,yMax] = ginput(1);
xlim([xMin,xMax])
disp('Select plastic peak')
[xPlastic,yPlastic] = ginput(1);
close(gcf);



% initialize
invertedVolData = zeros(yDim,xDim,zDim);
for x = 1:xDim
    for y = 1:yDim
        for z = 1:zDim
            count = count + 1;
            currentData = volData(y,x,z); % the current datapoint
            % store signal values only within the predefined window
            if (currentData >= xMin)&&(currentData <= xMax)
            	% invert the image
                invertedVolData(y,x,z) = -1*(currentData - xPlastic);
            end
        end
    end
end
% % for debugging 
% fileMap = cell(1,length(invertedVolData(1,1,:)));
% guiData{1} = invertedVolData;
% guiData{2} = fileMap;
% % display the GUI
% UserInputGUI2(guiData);


%% Old method
% % select a window of attenuation values and invert
% disp('Select window for attenuation values. (all outside regions will be set to zero)')
% disp('Select minimum window position')
% [xMin,yMin] = ginput(1);
% disp('Select maximum window position')
% [xMax,yMax] = ginput(1);
% close(gcf);
% 
% 
% 
% % initialize
% invertedVolData = zeros(yDim,xDim,zDim);
% for x = 1:xDim
%     for y = 1:yDim
%         for z = 1:zDim
%             count = count + 1;
%             currentData = volData(y,x,z); % the current datapoint
%             % store signal values only within the predefined window
%             if (currentData >= xMin)&&(currentData <= xMax)
%             	% invert the image
%                 invertedVolData(y,x,z) = abs(currentData - xMax);
%                 % invert the image, normalize the signal to range of 0 to 1
% %                 invertedVolData(y,x,z) = abs(currentData - xMax)/(xMax - xMin);
%             end
%         end
%     end
% end

% % for debugging 
% fileMap = cell(1,length(invertedVolData(1,1,:)));
% guiData{1} = invertedVolData;
% guiData{2} = fileMap;
% % display the GUI
% UserInputGUI2(guiData);


end