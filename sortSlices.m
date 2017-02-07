% Function to sort the data into individual slices so that a plane can be
% fitted to each slice individually. Note to user: the length of both the
% ground truth data and the sphere locations must be the same and they must
% correspond to one another (each datapoint in coordSphere must correspond
% to each datapoint in sphereGroundTruth)
% 
% Input:
% coordSphere [yInd,xInd,zInd] The coordinates of the sphere determined by the program
% sphereGoundTruth [yInd,xInd,zInd] The known, ground truth location of the sphere
%
% Output:
% sortedSphere{sphere,slice} [yInd,xInd,zInd] The coordinates of the
% spheres sorted into each slice
% sortedGndTruth{sphere,slice} [yInd,xInd,zInd] The coordinates of the
% groundTruth sorted into each slice
%
% John Ginn
% Created: 6/28/16
% Modified: 8/16/16
function [sortedSphere, sortedGndTruth] = sortSlices(coordSphere, sphereGroundTruth)

% all the GroundTruth points will be in the same plane
initialData = sphereGroundTruth{1};
% the start in slice for compairison
thisSlice = initialData(3);

sliceCount = 1; % counts the slices
sphereCount = 1; % counts spheres in each plane
sliceArray(sliceCount) = thisSlice; % the index of each slice in the volume
numSpheresInSlice(sliceCount) = sphereCount; % number of spheres in each slice
for step = 1:length(sphereGroundTruth)
    currentLocation = sphereGroundTruth{step};
    % in the same plane, add onto current array
    if thisSlice == currentLocation(3)
        sortedSphere{sphereCount,sliceCount} = coordSphere{step};
        sortedGndTruth{sphereCount,sliceCount} = sphereGroundTruth{step};
        numSpheresInSlice(sliceCount) = sphereCount; 
        sphereCount = sphereCount + 1;
    else
        % check to see if this slice already exists
        sliceExists = 0; % does the slice already exist y = 1, n = 0;
        for sliceStep = 1:length(sliceArray)
            % the slice already exists
            if sliceArray(sliceStep) == currentLocation(3)
                sliceExists = 1; 
                sliceCount = sliceStep; % the index of the slice
            end
        end
        % the slice already exists
        if sliceExists == 1
            % the number of spheres already in the slice
            sphereCount = numSpheresInSlice(sliceCount);
            % add on the new sphere to the array
            sphereCount = sphereCount + 1;
            sortedSphere{sphereCount,sliceCount} = coordSphere{step};
            sortedGndTruth{sphereCount,sliceCount} = sphereGroundTruth{step};
        else
            % in different plane, start new row in dataset
            sphereCount = 1; % restart sphere count
            sliceCount = sliceCount + 1;
            sortedSphere{sphereCount,sliceCount} = coordSphere{step};
            sortedGndTruth{sphereCount,sliceCount} = sphereGroundTruth{step};
        end
        % count another sphere so the one just added does not get deleted
        sphereCount = sphereCount + 1; 
    end
    % store slice to compare to see if the next sphere is in the same plane
    thisSlice = currentLocation(3);
end

end