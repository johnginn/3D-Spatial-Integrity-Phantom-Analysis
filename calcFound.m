% Function to calculate the number of spheres failing the correlation
% correlation threshold in the different regions after the shift correction
% is applied to their location.
%
% Input:
% sliceSphereX The x-location of the spheres in the slice extracted
% sliceSphereY The y-location of the spheres in the slice extracted
% sliceSphereZ The z-location of the spheres in the slice extacted
% sliceGroundX The x-location of the ground truth in the slices
% sliceGroundY The y-location of the ground truth in the slices
% sliceGroundZ The z-location of the ground truth in the slices
% voxelHeight (mm) The height of the voxels
% voxelWidth (mm) The width of the voxels
% voxelLength (mm) The length of the voxels
% shiftCorrection [sagittal,coronal,transverse] The correction required by fitting each center plane
% shiftForRot [column, row, slice] The shift to center the rotation (from plane fit typically)
% RotationAngle The angle of the rotation
% rotationMatrix The rotation matrix
% centerCol The index of the center column of the phantom
% centerSlice The index of the center slice of the phantom
% centerRow The index of the center row of the phantom
% radiusSearch1 (mm) The radius from isocenter that defines the subvolume for the first threshold
% radiusSearch2 (mm) The radius from isocenter that defines the subvolume for the second threshold
% radiusSearch3 (mm) The radius from isocenter that defines the subvolume for the third threshold
% radiusSearch4 (mm) The radius from isocenter that defines the subvolume for the fourth threshold
% rotTransFinal The transverse rotational matrix that minimizes deviation from ground-truth locations 
% rotCoronalFinal The coronal rotational matrix that minimizes deviation from ground-truth locations 
% rotSagittalFinal The sagittal rotational matrix that minimizes deviation from ground-truth locations 
% shiftForRotTrans The shift required to center the locations at isocenter before the transverse rotation
% shiftForRotCoronal The shift required to center the locations at isocenter before the coronal rotation
% shiftForRotSagittal The shift required to center the locations at isocenter before the sagittal rotation
%
% Output:
% failSphRadiusSearch1 The number of spheres in RadiusSearch1 failing the correlation threshold
% failSphRadiusSearch2 The number of spheres in RadiusSearch2 failing the correlation threshold
% failSphRadiusSearch3 The number of spheres in RadiusSearch3 failing the correlation threshold
% failSphRadiusSearch4 The number of spheres in RadiusSearch4 failing the correlation threshold
%
% John Ginn
% Created: 7/21/16
% Modified: 1/24/17
function [failSphRadiusSearch1, failSphRadiusSearch2, failSphRadiusSearch3,failSphRadiusSearch4] = ...
    calcFound(xSphereFail,ySphereFail,zSphereFail,xGndTruthFail,yGndTruthFail,zGndTruthFail,...
    voxelHeight,voxelWidth,voxelLength,centerCol,centerSlice,centerRow,shiftCorrection,radiusSearch1,radiusSearch2,radiusSearch3,...
    radiusSearch4,rotTransFinal,rotCoronalFinal,rotSagittalFinal,shiftForRotTrans,shiftForRotCoronal,shiftForRotSagittal,...
    RotationAngleTrans,RotationAngleCoronal,RotationAngleSagittal)

% the rotation matricies 
rotationTrans = rotTransFinal(RotationAngleTrans); 
rotationCoronal = rotCoronalFinal(RotationAngleCoronal);
rotationSagittal = rotSagittalFinal(RotationAngleSagittal);
% shift the vector prior to rotation and after rotation to determine axis
% of rotation
shiftVectorX = @(shift,vector) [vector(1) + shift,vector(2),vector(3)];
shiftVectorY = @(shift,vector) [vector(1),vector(2) + shift,vector(3)];
shiftVectorZ = @(shift,vector) [vector(1),vector(2),vector(3) + shift];


% calcualte the deviation for the center transverse slice to compare to the
% ViewRay analysis software
radius = zeros(1,length(xSphereFail));
radiusGnd = zeros(1,length(xSphereFail));

% number of spheres that did not meet correlation coefficient tolerance
failSphRadiusSearch1 = 0; % sphere location within first radius
failSphRadiusSearch2 = 0; % sphere location within second radius
failSphRadiusSearch3 = 0; % sphere location within third radius
failSphRadiusSearch4 = 0; % sphere location within fourth radius

failGndRadiusSearch1 = 0; % gnd-truth location within first radius
failGndRadiusSearch2 = 0; % gnd-truth location within second radius
failGndRadiusSearch3 = 0; % gnd-truth location within third radius
failGndRadiusSearch4 = 0; % gnd-truth location within fourth radius
% relative to isocenter
for step = 1:length(xSphereFail)
    % rotate and shift the data with the corrections calculated from the
    % passing dataset
    % transverse rotation correction
    SphereFailComp(step,:) = [xSphereFail(step), ySphereFail(step), zSphereFail(step)];
    % shift the vector so the ceter column is at the origin
    SphereFailComp(step,:)= shiftVectorX(-shiftForRotTrans(1),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorY(-shiftForRotTrans(2),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorZ(-shiftForRotTrans(3),SphereFailComp(step,:));
    % now at origin, apply rotation
    SphereFailComp(step,:) = rotationTrans*SphereFailComp(step,:)';
    % shift back
    SphereFailComp(step,:)= shiftVectorX(shiftForRotTrans(1),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorY(shiftForRotTrans(2),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorZ(shiftForRotTrans(3),SphereFailComp(step,:));
    
    % coronal rotation correction
    SphereFailComp(step,:)= shiftVectorX(-shiftForRotCoronal(1),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorY(-shiftForRotCoronal(2),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorZ(-shiftForRotCoronal(3),SphereFailComp(step,:));
    SphereFailComp(step,:) = rotationCoronal*SphereFailComp(step,:)';
    SphereFailComp(step,:)= shiftVectorX(shiftForRotCoronal(1),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorY(shiftForRotCoronal(2),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorZ(shiftForRotCoronal(3),SphereFailComp(step,:));
    
    % sagittal rotation correction
    SphereFailComp(step,:)= shiftVectorX(-shiftForRotSagittal(1),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorY(-shiftForRotSagittal(2),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorZ(-shiftForRotSagittal(3),SphereFailComp(step,:));
    SphereFailComp(step,:) = rotationSagittal*SphereFailComp(step,:)';
    SphereFailComp(step,:)= shiftVectorX(shiftForRotSagittal(1),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorY(shiftForRotSagittal(2),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorZ(shiftForRotSagittal(3),SphereFailComp(step,:));
    
    
    % apply the shift correction from translational correction
    SphereFailComp(step,:) = shiftVectorX(shiftCorrection(1),SphereFailComp(step,:));
    SphereFailComp(step,:) = shiftVectorY(shiftCorrection(2),SphereFailComp(step,:));
    SphereFailComp(step,:)= shiftVectorZ(shiftCorrection(3),SphereFailComp(step,:));

    xSphereFailComp(step) = SphereFailComp(1);
    ySphereFailComp(step) = SphereFailComp(2);
    zSphereFailComp(step) = SphereFailComp(3);
    
    % current distance from isocenter
    radius(step) = sqrt((voxelWidth*(xSphereFailComp(step) - centerCol)).^2 + ...
        (voxelHeight*(ySphereFailComp(step) - centerRow)).^2 +...
        (voxelLength*(zSphereFailComp(step) - centerSlice)).^2);
    
    radiusGnd(step) = sqrt((voxelWidth*(xGndTruthFail(step) - centerCol)).^2 + ...
        (voxelHeight*(yGndTruthFail(step) - centerRow)).^2 +...
        (voxelLength*(zGndTruthFail(step) - centerSlice)).^2);
    
    
    % count the number of spheres that fail within the specified radius of
    % isocenter
    if radius(step) < radiusSearch1
        failSphRadiusSearch1 = failSphRadiusSearch1 + 1;
    end
    if radius(step) < radiusSearch2
        failSphRadiusSearch2 = failSphRadiusSearch2 + 1;
    end
    if radius(step) < radiusSearch3
        failSphRadiusSearch3 = failSphRadiusSearch3 + 1;
    end
    if radius(step) < radiusSearch4
        failSphRadiusSearch4 = failSphRadiusSearch4 + 1;
    end
    if radiusGnd(step) < radiusSearch1
        failGndRadiusSearch1 = failGndRadiusSearch1 + 1;
    end
    if radiusGnd(step) < radiusSearch2
        failGndRadiusSearch2 = failGndRadiusSearch2 + 1;
    end
    if radiusGnd(step) < radiusSearch3
        failGndRadiusSearch3 = failGndRadiusSearch3 + 1;
    end
    if radiusGnd(step) < radiusSearch4
        failGndRadiusSearch4 = failGndRadiusSearch4 + 1;
    end
end

end