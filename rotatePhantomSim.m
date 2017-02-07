% Function to rotate the datapoints from the phantom simulation to
% determine if the roation fitting is working appropriately
% 
% Input:
% xData The x-location data of the spheres to be rotated
% yData The y-location data of the spheres to be rotated
% zData The z-location data of the spheres to be rotated
% angleTrans (degrees) The angle to rotate the data in the transverse plane
% angleCoronal (degrees) The angle to rotate the data in the coronal plane
% angleSagittal (degrees) The angle to rotate the data in the sagittal plane
% directionTrans The direction to rotate the data in the transverse plane 1 = clockwise, 0 = counter-clokwise
% directionCoronal The direction to rotate the data in the coronal plane 1 = clockwise, 0 = counter-clokwise
% directionSagittal The direction to rotate the data in the sagittal plane 1 = clockwise, 0 = counter-clokwise
% centerRow The center row of the phantom
% centerCol The center column of the phantom
% centerSlice The center slice of the phantom
% xShift An additional x-shift you want to apply 
% yShift An additional y-shift you want to apply
% zShift An additional z-shift you want to apply 
% troubleshoot Plot for troubleshooting 1 = y, 0 = n
%
% Output:
% xRotData The x-component of the rotated data
% yRotData The y-component of the rotated data
% zRotData The z-component of the rotated data
%
% John Ginn
% Created: 7/5/16
% Modified: 8/16/16
function [xRotData,yRotData,zRotData] = rotatePhantomSim(xData,yData,zData,...
    angleTrans,angleCoronal,angleSagittal,directionTrans,directionCoronal,...
    directionSagittal,centerRow,centerCol,centerSlice,xShift,yShift,zShift,troubleshoot)


% matrices for rotating the data (select clockwise vs. counter-clockwise
% rotation)
% transverse
if directionTrans == 1;
    rotMatrixTrans = @(theta)[cosd(theta),  0,    -sind(theta);...
                                            0,         1,          0;
                                         sind(theta),  0,     cosd(theta)];
else
    rotMatrixTrans =  @(theta)[cosd(theta),  0,    sind(theta);...
                                         0,         1,          0;
                                     -sind(theta),  0,     cosd(theta)];
end
% coronal
if directionCoronal == 1;
    rotMatrixCoronal = @(theta)[1,            0,              0;...
                                      0,        cosd(theta),    -sind(theta);
                                      0,        sind(theta),     cosd(theta)];
else
    rotMatrixCoronal =  @(theta)[1,            0,              0;...
                                         0,        cosd(theta),    sind(theta);
                                         0,        -sind(theta),   cosd(theta)];
end
% sagittal
if directionSagittal == 1;
    rotMatrixSagittal = @(theta)[cosd(theta),  -sind(theta),    0;...
                                      sind(theta),   cosd(theta),    0;
                                          0,             0,          1];
else
    rotMatrixSagittal = @(theta)[cosd(theta),    sind(theta),    0;...
                                        -sind(theta),   cosd(theta),    0;
                                             0,             0,          1];
end
% the final rotation matrices
rotationTrans = rotMatrixTrans(angleTrans);
rotationCoronal = rotMatrixCoronal(angleCoronal);
rotationSagittal = rotMatrixSagittal(angleSagittal);


% shifts
shiftVectorX = @(shift,vector) [vector(1) + shift,vector(2),vector(3)];
shiftVectorY = @(shift,vector) [vector(1), vector(2) + shift, vector(3)];
shiftVectorZ = @(shift,vector) [vector(1),vector(2),vector(3) + shift];

for step = 1:length(xData)
    SphereRotated(step,:) = [xData(step), yData(step), zData(step)];
    % shift the vector so the ceter column is at the origin
    SphereRotated(step,:)= shiftVectorX(-centerCol,SphereRotated(step,:));
    SphereRotated(step,:)= shiftVectorY(-centerRow,SphereRotated(step,:));
    SphereRotated(step,:)= shiftVectorZ(-centerSlice,SphereRotated(step,:));
    % apply the rotations
    centeredData(step,:) = SphereRotated(step,:);
    SphereRotated(step,:) = rotationTrans*SphereRotated(step,:)';
    SphereRotated(step,:) = rotationCoronal*SphereRotated(step,:)';
    SphereRotated(step,:) = rotationSagittal*SphereRotated(step,:)';
    centeredDataRotated(step,:) = SphereRotated(step,:);
    % shift back and apply additional prescribed shift
    SphereRotated(step,:) = shiftVectorX(centerCol + xShift,SphereRotated(step,:));
    SphereRotated(step,:)= shiftVectorY(centerRow + yShift,SphereRotated(step,:));
    SphereRotated(step,:)= shiftVectorZ(centerSlice + zShift,SphereRotated(step,:));
    % does not use data that does not have high enough correlation coeff
    % location of spheres from fit [y, x, z]
    xRotData(step) = SphereRotated(step,1); % sphere location x-component
    yRotData(step) = SphereRotated(step,2); % sphere location y-component
    zRotData(step) = SphereRotated(step,3); % sphere location z-component
    
end

% for troubleshooting
if troubleshoot == 1;
    figure;
    scatter3(centeredData(:,1),centeredData(:,2),centeredData(:,3))
    hold on
    scatter3(centeredDataRotated(:,1),centeredDataRotated(:,2),centeredDataRotated(:,3),'r+')
    hold on
    scatter3(xRotData,yRotData,zRotData,'k.')
    title('Troubleshoot Sim Phantom Rotation','FontSize',20)
    legend('centered not-rotated','centered rotated','final')
    xlabel('x-axis','FontSize',20)
    ylabel('y-axis','FontSize',20)
    zlabel('z-axis','FontSize',20)
end

end