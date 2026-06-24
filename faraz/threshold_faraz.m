function [output, ind] = threshold_faraz(input,threshold,dir)
% threshold: fraction of max value to be used for thresholding
%            OR 'hist' to threshold based off of histogram
% dir: whether to threshold greater or less than the provided value

% output: thresholded image
% ind: indices for the non-thresholded points (broken)

% does NOT handle negative values because abs is used (for complex values)
% does NOT modify input if threshold = 0

% check inputs
if nargin < 2, threshold = []; end
if nargin < 3, dir = []; end
if isempty(dir), dir = 'less'; end

% define initial output
output = input;

if ischar(threshold) || isempty(threshold)
    for j = 1:size(input,4)
        for i = 1:size(input,3)
            tmp = input(:,:,i,j);
            [~,edges] = histcounts(tmp);
            tmp(tmp < edges(2)) = 0;
            output(:,:,i,j) = tmp;
        end
    end
else
    if threshold ~= 0
        switch lower(dir)
            case 'less'
                for j = 1:size(input,4)
                    temp = output(:,:,:,j);
                    ind = find(abs(temp) >= threshold*max(abs(reshape(temp,[],1))));
                    temp(abs(temp) < threshold*max(abs(reshape(temp,[],1)))) = 0;
                    output(:,:,:,j) = temp;
                end
%                 ind = find(abs(output)>=threshold*max(abs(output(:))));
%                 output(abs(output)<threshold*max(abs(output(:)))) = 0;
            case 'greater'
                for j = 1:size(input,4)
                    temp = output(:,:,:,j);
                    ind = find(abs(temp) < threshold*max(abs(reshape(temp,[],1))));
                    temp(abs(temp) >= threshold*max(abs(reshape(temp,[],1)))) = 0;
                    output(:,:,:,j) = temp;
                end
%                 ind = find(abs(output)<threshold*max(abs(output(:))));
%                 output(abs(output)>=threshold*max(abs(output(:)))) = 0;
        end
    end
end
end