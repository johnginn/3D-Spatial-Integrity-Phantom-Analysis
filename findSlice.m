% function to create an array of the slice locations. Note to user: this
% function cannot be used after the rotation has been applied, otherwise
% the datapoints will no longer lie in the same plane
%
% Input:
% zLocation The z-position
%
% Output: 
% sliceNumber The slice number in the phantom
%
% John Ginn
% Created: 7/1/16
% Modified: 7/1/16
function [sliceNumber] = findSlice(zLocation)
slicePosition = zLocation(1); % the first slice
[slicePosition, sliceNumber] = sort(slicePosition);
    for x = 1:length(zLocation)
        for searchSlices = 1:length(slicePosition)
             if zLocation(x) == slicePosition(searchSlices)
                 sliceFound = 1;
             end
        end 
    end
        
end