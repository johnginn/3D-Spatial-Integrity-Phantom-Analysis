% This function saves the different variables obtained from the analysis to
% excel files. WARNING: Be sure to transfer over data from a previous
% analysis prior to running this function, otherwise the data will be
% overwritten
%
% Input:
% gantryAngle The angle of the gantry during the can
% scan The title of the scan file for referencing scanning parameters
% numSpheresCalDevExp Number of spheres found surrounding centermost sphere used to calculate uncertainty
% rmsDevFromExp Root-mean square of the deviation from expected distance for spheres found surrounding centermost sphere
% rmsDevFromGnd Root-mean square of the deviation from ground-truth locations for spheres surrounding centermost sphere
% corrThreshold Correlation coefficient threshold to reject spheres in artifacts
% spheresFoundAll Percent of spheres passing correlation threshold
% radiusSearch1 Radius for user defined region 1 of analysis
% spheresFoundThresh1 Deviation threshold for spheres within radius1
% radiusSearch2 Radius for user defined region 1 of analysis
% spheresFoundThresh2 Deviation threshold for spheres within radius2
% radiusSearch3 Radius for user defined region 3 of analysis
% spheresFoundThresh3 Deviation threshold for spheres within radius3
% avgDiff2DCentralAxis1 Average 2D deviation within radius 1 of central axis
% avgDiff2DCentralAxis2 Average 2D deviation within radius 2 of central axis
% avgDiff2DCentralAxis3 Average 2D deviation within radius 3 of central axis
% avgDiff3DRadius1 Average 3D deviation for spheres within radius 1
% avgDiff2DRadius1 Average 2D deviation for spheres within radius 1
% avgDiff3DRadius2 Average 3D deviation for spheres within radius 2
% avgDiff2DRadius2 Average 2D deviation for spheres within radius 2
% avgDiff3DRadius3 Average 3D deviation for spheres within radius 3
% avgDiff2DRadius3 Average 2D deviation for spheres within radius 3
% avgDiff3DAll Average 3D deviation within the entire phantom
% avgDiff2DAll Average 2D deviation within the entire phantom
% radiusThreshold1 User-defined deviation threshold for spheres within radius1
% radiusThreshold2 User-defined deviation threshold for spheres within radius2
% radiusThreshold3 User-defined deviation threshold for spheres within radius3
% pass3DThresh1 Percent of spheres passing with a 3D deviation less than threshold1 in radius1
% pass2DThresh1 Percent of spheres passing with a 2D deviation less than threshold1 in radius1
% pass3DThresh2 Percent of spheres passing with a 3D deviation less than threshold2 in radius2
% pass2DThresh2 Percent of spheres passing with a 2D deviation less than threshold2 in radius2
% pass3DThresh3 Percent of spheres passing with a 3D deviation less than threshold2 in radius3
% pass2DThresh3 Percent of spheres passing with a 2D deviation less than threshold2 in radius3
% vrPercentPassXYZ1 Centermost axial-slice, spheres with a 3D deviation less than threshold1 in radius1
% vrPercentPassXY1 Centermost axial-slice, spheres with a 2D deviation less than threshold1 in radius1
% vrAvgDevRadius1XY Centermost axial-slice, average 2D deviation of spheres within radius1
% vrPercentPassXYZ2 Centermost axial-slice, spheres with a 3D deviation less than threshold2 in radius2
% vrPercentPassXY2 Centermost axial-slice, spheres with a 2D deviation less than threshold2 in radius2
% vrAvgDevRadius2XY  Centermost axial-slice, average 2D deviation of spheres within radius2
% vrPercentPassXYZ3 Centermost axial-slice, spheres with a 3D deviation less than threshold3 in radius3
% vrPercentPassXY3 Centermost axial-slice, spheres with a 2D deviation less than threshold3 in radius3
% vrAvgDevRadius3XY  Centermost axial-slice, average 2D deviation of spheres within radius3
% pass2mmDevPercent Percent of spheres in entire volume with less than 2mm deviation from ground truth
% pass3mmDevPercent Percent of spheres in entire volume with less than 3mm deviation from ground truth
% finalMaxDeviation3D Max 3D deviation for a single sphere in entire volume
% finalMaxDeviation2D Max 2D deviation for a single sphere in entire volume
% maxDiff3DThresh1 Max 3D deviation for a single sphere in radius1 of isocenter
% maxDiff2DThresh1 Max 2D deviation for a single sphere in radius1 of isocenter
% maxDiff3DThresh2 Max 3D deviation for a single sphere in radius2 of isocenter
% maxDiff2DThresh2 Max 2D deviation for a single sphere in radius2 of isocenter
% maxDiff3DThresh3 Max 3D deviation for a single sphere in radius3 of isocenter
% maxDiff2DThresh3 Max 2D deviation for a single sphere in radius3 of isocenter
% SNR Signal to noise ratio calculated for centermost sphere
% avgSig Average signal for centermost sphere
% bgStandardDev Standard deviation of the background signal
% allThreshold The passing threshold for all spheres in the volume
%
% Output:
%
% John Ginn
% Created: 9/2/16
% Modified: 9/6/16

function [] = saveExcel(gantryAngle,scan,numSpheresCalDevExp,rmsDevFromExp,...
    rmsDevFromGnd,corrThreshold,spheresFoundAll,radiusSearch1,spheresFoundThresh1,radiusSearch2,...
    spheresFoundThresh2,radiusSearch3,spheresFoundThresh3,avgDiff2DCentralAxis1,avgDiff2DCentralAxis2,avgDiff2DCentralAxis3,avgDiff3DRadius1,...
    avgDiff2DRadius1,avgDiff3DRadius2,avgDiff2DRadius2,avgDiff3DRadius3,avgDiff2DRadius3,avgDiff3DAll,avgDiff2DAll,...
    radiusThreshold1,radiusThreshold2,radiusThreshold3,pass3DThresh1,pass2DThresh1,pass3DThresh2,pass2DThresh2,pass3DThresh3,pass2DThresh3,...
    vrPercentPassXYZ1,vrPercentPassXY1,vrAvgDevRadius1XY,vrPercentPassXYZ2,vrPercentPassXY2,...
    vrAvgDevRadius2XY,vrPercentPassXYZ3,vrPercentPassXY3,vrAvgDevRadius3XY,passAll2DPercent,passAll3DPercent,finalMaxDeviation3D,finalMaxDeviation2D,...
    maxDiff3DThresh1,maxDiff2DThresh1,maxDiff3DThresh2,maxDiff2DThresh2,maxDiff3DThresh3,maxDiff2DThresh3,SNR,avgSig,...
    bgStandardDev,radiusSearch4,radiusThreshold4,spheresFoundThresh4,avgDiff2DCentralAxis4,avgDiff3DRadius4,avgDiff2DRadius4,...
    pass3DThresh4,pass2DThresh4,vrPercentPassXYZ4,vrPercentPassXY4,vrAvgDevRadius4XY,maxDiff3DThresh4,maxDiff2DThresh4,allThreshold)

% The old way of storing the data
oldData = {'Scan Title',scan,' ';...
    'Gantry Angle',gantryAngle,' (degrees)';
    strcat(['Deviation from expected distance RMS for the ',num2str(numSpheresCalDevExp),' spheres surrounding isocenter: ']),num2str(rmsDevFromExp),'(mm)';...
    strcat(['Deviation from ground truth RMS for the ',num2str(numSpheresCalDevExp+1),' spheres nearest to isocenter: ']),num2str(rmsDevFromGnd),'(mm)';...
strcat(['Spheres found with corr. threshold of ',num2str(corrThreshold),' in a ',num2str(radiusSearch1),' (mm) radius: ']),num2str(spheresFoundThresh1),'%';...
strcat(['Spheres found with corr. threshold of ',num2str(corrThreshold),' in a ',num2str(radiusSearch2),' (mm) radius: ']),num2str(spheresFoundThresh2),'%';...
strcat(['Spheres found with corr. threshold of ',num2str(corrThreshold),' in a ',num2str(radiusSearch4),' (mm) radius: ']),num2str(spheresFoundThresh4),'%';...
strcat(['Spheres found with corr. threshold of ',num2str(corrThreshold),' in a ',num2str(radiusSearch3),' (mm) radius: ']),num2str(spheresFoundThresh3),'%';...
strcat(['Spheres found  with corr. threshold of ',num2str(corrThreshold),' in entire volume: ']),num2str(spheresFoundAll),'%';...
strcat([num2str(radiusSearch1),' (mm) central axis dist, avg. 2D deviation: ']), num2str(avgDiff2DCentralAxis1), ' (mm)';...
strcat([num2str(radiusSearch2),' (mm) central axis dist, avg. 2D deviation: ']), num2str(avgDiff2DCentralAxis2), ' (mm)';...
strcat([num2str(radiusSearch4),' (mm) central axis dist, avg. 2D deviation: ']), num2str(avgDiff2DCentralAxis4), ' (mm)';...
strcat([num2str(radiusSearch3),' (mm) central axis dist, avg. 2D deviation: ']), num2str(avgDiff2DCentralAxis3), ' (mm)';...
strcat([num2str(radiusSearch1),' (mm) radius volume avg. 3D deviation: ']), num2str(avgDiff3DRadius1), ' (mm)';...
strcat([num2str(radiusSearch1),' (mm) radius volume avg. 2D deviation: ']), num2str(avgDiff2DRadius1), ' (mm)';...
strcat([num2str(radiusSearch2),' (mm) radius volume avg. 3D deviation: ']), num2str(avgDiff3DRadius2), ' (mm)';...
strcat([num2str(radiusSearch2),' (mm) radius volume avg. 2D deviation: ']), num2str(avgDiff2DRadius2), ' (mm)';...
strcat([num2str(radiusSearch4),' (mm) radius volume avg. 3D deviation: ']), num2str(avgDiff3DRadius4), ' (mm)';...
strcat([num2str(radiusSearch4),' (mm) radius volume avg. 2D deviation: ']), num2str(avgDiff2DRadius4), ' (mm)';...
strcat([num2str(radiusSearch3),' (mm) radius volume avg. 3D deviation: ']), num2str(avgDiff3DRadius3), ' (mm)';...
strcat([num2str(radiusSearch3),' (mm) radius volume avg. 2D deviation: ']), num2str(avgDiff2DRadius3), ' (mm)';...
strcat('Whole volume avg. 3D deviation: '), num2str(avgDiff3DAll), ' (mm)';...
strcat('Whole volume avg. 2D deviation: '), num2str(avgDiff2DAll), ' (mm)';...
strcat(['3D deviation, percent passing in ',num2str(radiusSearch1),' (mm) radius volume, tolerance of ',num2str(radiusThreshold1),' (mm): ']), num2str(pass3DThresh1),'%';...
strcat(['2D deviation, percent passing in ',num2str(radiusSearch1),' (mm) radius volume, tolerance of ',num2str(radiusThreshold1),' (mm): ']), num2str(pass2DThresh1),'%';...
strcat(['3D deviation, percent passing in ',num2str(radiusSearch2),' (mm) radius volume, tolerance of ',num2str(radiusThreshold2),' (mm): ']), num2str(pass3DThresh2),'%';...
strcat(['2D deviation, percent passing in ',num2str(radiusSearch2),' (mm) radius volume, tolerance of ',num2str(radiusThreshold2),' (mm): ']), num2str(pass2DThresh2),'%';...
strcat(['3D deviation, percent passing in ',num2str(radiusSearch4),' (mm) radius volume, tolerance of ',num2str(radiusThreshold4),' (mm): ']), num2str(pass3DThresh4),'%';...
strcat(['2D deviation, percent passing in ',num2str(radiusSearch4),' (mm) radius volume, tolerance of ',num2str(radiusThreshold4),' (mm): ']), num2str(pass2DThresh4),'%';...
strcat(['3D deviation, percent passing in ',num2str(radiusSearch3),' (mm) radius volume, tolerance of ',num2str(radiusThreshold3),' (mm): ']), num2str(pass3DThresh3),'%';...
strcat(['2D deviation, percent passing in ',num2str(radiusSearch3),' (mm) radius volume, tolerance of ',num2str(radiusThreshold3),' (mm): ']), num2str(pass2DThresh3),'%';...
strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch1),' (mm) dev <= ',num2str(radiusThreshold1),' (mm), all 3D: ']),num2str(vrPercentPassXYZ1),'%';...
strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch1),' (mm) dev <= ',num2str(radiusThreshold1),' (mm), just x,y: ']),num2str(vrPercentPassXY1),'%';...
strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch1),' (mm)',' avg 2D deviation: ']),num2str(vrAvgDevRadius1XY),' (mm)';...
strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch2),' (mm) dev <= ',num2str(radiusThreshold2),' (mm), all 3D: ']),num2str(vrPercentPassXYZ2),'%';...
strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch2),' (mm) dev <= ',num2str(radiusThreshold2),' (mm), just x,y: ']),num2str(vrPercentPassXY2),'%';...
strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch2),' (mm)',' avg 2D deviation: ']),num2str(vrAvgDevRadius2XY),' (mm)';...
strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch4),' (mm) dev <= ',num2str(radiusThreshold4),' (mm), all 3D: ']),num2str(vrPercentPassXYZ4),'%';...
strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch4),' (mm) dev <= ',num2str(radiusThreshold4),' (mm), just x,y: ']),num2str(vrPercentPassXY4),'%';...
strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch4),' (mm)',' avg 2D deviation: ']),num2str(vrAvgDevRadius4XY),' (mm)';...
strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch3),' (mm) dev <= ',num2str(radiusThreshold3),' (mm), all 3D: ']),num2str(vrPercentPassXYZ3),'%';...
strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch3),' (mm) dev <= ',num2str(radiusThreshold3),' (mm), just x,y: ']),num2str(vrPercentPassXY3),'%';...
strcat(['Single slice to compare to ViewRay w/in ',num2str(radiusSearch3),' (mm)',' avg 2D deviation: ']),num2str(vrAvgDevRadius3XY),' (mm)';...
strcat(['Percent of spheres in entire volume with less ',num2str(allThreshold),' mm 2D deviation: ']),num2str(passAll2DPercent),'%';...
strcat(['Percent of spheres in entire volume with less than ',num2str(allThreshold),' mm 3D deviation: ']),num2str(passAll3DPercent),'%';...
strcat(['The maximum 3D deviation of a sphere in the entire volume: ']),num2str(finalMaxDeviation3D),' (mm)';...
strcat(['The maximum 2D deviation of a sphere in the entire volume: ']),num2str(finalMaxDeviation2D),' (mm)';...
strcat(['The maximum 3D deviation in ',num2str(radiusSearch1),' (mm) radius volume: ']),num2str(maxDiff3DThresh1),' (mm)';...
strcat(['The maximum 2D deviation in ',num2str(radiusSearch1),' (mm) radius volume: ']),num2str(maxDiff2DThresh1),' (mm)';...
strcat(['The maximum 3D deviation in ',num2str(radiusSearch2),' (mm) radius volume: ']),num2str(maxDiff3DThresh2),' (mm)';...
strcat(['The maximum 2D deviation in ',num2str(radiusSearch2),' (mm) radius volume: ']),num2str(maxDiff2DThresh2),' (mm)';...
strcat(['The maximum 3D deviation in ',num2str(radiusSearch4),' (mm) radius volume: ']),num2str(maxDiff3DThresh4),' (mm)';...
strcat(['The maximum 2D deviation in ',num2str(radiusSearch4),' (mm) radius volume: ']),num2str(maxDiff2DThresh4),' (mm)';...
strcat(['The maximum 3D deviation in ',num2str(radiusSearch3),' (mm) radius volume: ']),num2str(maxDiff3DThresh3),' (mm)';...
strcat(['The maximum 2D deviation in ',num2str(radiusSearch3),' (mm) radius volume: ']),num2str(maxDiff2DThresh3),' (mm)';...
'Center Sphere SNR: ',num2str(SNR),' (au)';...
'Average Signal: ',num2str(avgSig),' (au)';...
'Bg. Signal Standard Dev: ',num2str(bgStandardDev),' (au)'};

% table containing the new data
newData = {'Scan Title','Gantry Angle (degrees)',strcat('Whole volume avg. 3D deviation: '),...
    strcat('Whole volume avg. 2D deviation: '),strcat([num2str(radiusSearch1),' (mm) radius volume avg. 3D deviation: ']),...
strcat([num2str(radiusSearch1),' (mm) radius volume avg. 2D deviation: ']),...
strcat([num2str(radiusSearch2),' (mm) radius volume avg. 3D deviation: ']),...
strcat([num2str(radiusSearch2),' (mm) radius volume avg. 2D deviation: ']),...
strcat([num2str(radiusSearch4),' (mm) radius volume avg. 3D deviation: ']),...
strcat([num2str(radiusSearch4),' (mm) radius volume avg. 2D deviation: ']),...
strcat([num2str(radiusSearch3),' (mm) radius volume avg. 3D deviation: ']),...
strcat([num2str(radiusSearch3),' (mm) radius volume avg. 2D deviation: ']);...
' ','(degrees)','(mm)','(mm)','(mm)','(mm)','(mm)','(mm)','(mm)','(mm)','(mm)','(mm)';...
scan,num2str(gantryAngle),num2str(avgDiff3DAll),num2str(avgDiff2DAll),...
num2str(avgDiff3DRadius1),num2str(avgDiff2DRadius1),num2str(avgDiff3DRadius2),...
num2str(avgDiff2DRadius2),num2str(avgDiff3DRadius4),num2str(avgDiff2DRadius4),...
num2str(avgDiff3DRadius3),num2str(avgDiff2DRadius3)};
% save the data
filename1 = 'oldDataReportStyle.xlsx';
xlswrite(filename1,oldData)

filename2 = 'newTable.xlsx';
xlswrite(filename2,newData);

end