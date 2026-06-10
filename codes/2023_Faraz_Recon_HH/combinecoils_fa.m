function [ImgArrayCombined,KspaceArrayCombined,b] = combinecoils_fa(ImgArray,KspaceArray,NoiseSize,b)
% Array In:         [ncol,nlin,nslices,nchannels,nrepetitions]
% Array Out:        [ncol,nlin,nslices,nrepetitions]

% check inputs
if nargin < 2, KspaceArray = []; end
if nargin < 3, NoiseSize = []; end
if nargin < 4, b = []; end
if isempty(KspaceArray), KspaceArrayCombined = []; end
if isempty(NoiseSize), NoiseSize = 6; end

% get dimensions
dim = size(ImgArray);
    
% fix for images with no repetitions
try 
    dim(5) = dim(5);
catch
    dim(5) = 1;
end

% define noise regions
NoiseRegionX = (dim(1)-NoiseSize):(dim(1)-1);
NoiseRegionY = (dim(2)-NoiseSize):(dim(2)-1);

% define NoiseRegionZ if 3D/iso
if dim(1) == dim(2) && dim(2) == dim(3)
    NoiseRegionZ = (dim(3)-NoiseSize):(dim(3)-1);
else
    NoiseRegionZ = 1:dim(3);
end
    
% normalize each channel by its noise
for ic = 1:dim(4)
    ImgArray(:, :, :, ic, :) = ImgArray(:, :, :, ic, :) / mean(abs(reshape(ImgArray(NoiseRegionX,NoiseRegionY,NoiseRegionZ,ic,:),[],1)));
end

if isempty(b)
    % sum all repetitions for calculating coil sensitivity
    b = sum(ImgArray, 5);
    
    % normalize by image norms
    b = b ./ sum(abs(ImgArray),5);
    
    % reshape
    b = permute(b, [4 1 2 3]);
end

% reshape
ImgArray = permute(ImgArray, [4 1 2 3 5]);

% combine optimally as per Bydder et al, MRM 47:539-458 (2002)
ImgArrayCombined = zeros(dim(1), dim(2), dim(3), dim(5));
for ix = 1:dim(1)
    for iy = 1:dim(2)
        for is = 1:dim(3)
            for ii = 1:dim(5)
                thisb = b(:, ix, iy, is);
                ImgArrayCombined(ix, iy, is, ii) = real((thisb' * ImgArray(:, ix, iy, is, ii)));
            end
        end
    end
end

% normalize by its noise
ImgArrayCombined = ImgArrayCombined / mean(reshape(ImgArrayCombined(NoiseRegionX,NoiseRegionY,NoiseRegionZ,:),[],1)); 

% zero negatives
% ImgArrayCombined(ImgArrayCombined < 0) = 0;

% just add up kspace channels
if ~isempty(KspaceArray)
    KspaceArrayCombined = sum(KspaceArray,4);
    KspaceArrayCombined = reshape(KspaceArrayCombined,size(KspaceArray,[1 2 3 5]));
end

end