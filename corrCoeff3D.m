% This function calculates the normalized cross correlation coefficient in
% 3D for finding the spheres in the IMT spatial integrity phantom.(see Phantom-based 
% characterization of distortion on a magnetic resonance imaging simulator 
% for radiation oncology. Huang et. al.) This function is computationally
% less expensive than MATLAB's correlation coeffienct function
%
% Input:
% imgVolume The image volume 
% tempVolume The template volume
%
% Output:
% gamma The normalized ccross correlation coefficient
%
% John Ginn
% Created: 11/22/16
% Modified: 11/22/16

function [gamma] = corrCoeff3D(imgVolume,tmpVolume)

imgArray(:) = imgVolume;
tmpArray(:) = tmpVolume;

% compute avg signal in the volume
sumImg = sum(imgArray);
sumTmp = sum(tmpArray);
nDataImg = length(imgArray);
nDataTmp = length(tmpArray);
% average signal values 
avgImg = sumImg/nDataImg;
avgTmp = sumTmp/nDataTmp;

% initialize the sums (see Phantom-based characterization of distortion on
% a magnetic resonance imaging simulator for radiation oncology. Huang et. al.)
numer = sum((imgArray - avgImg).*(tmpArray - avgTmp));
denom1 = sum((imgArray - avgImg).^2);
denom2 = sum((tmpArray - avgTmp).^2);

gamma = (numer/sqrt(denom1*denom2) + 1)/2;

end