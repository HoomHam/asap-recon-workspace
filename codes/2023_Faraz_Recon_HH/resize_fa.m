function out = resize_fa(img,scale,method,threshold)
% img:              absolute value of image matrix
% scale:            interpolation factor
% method:           interpolation method
% threshold:        threshold for masking edges

% check inputs
if nargin < 3, method = []; end
if nargin < 4, threshold = []; end
if isempty(method), method = 'bilinear'; end
if isempty(threshold), threshold = 0.5; end

% do nothing if scale = 1
if scale == 1
    out = img;
else
    
% create mask
if any(isnan(img(:))) 
    % include 0s if background is NaNs 
    mask = double(abs(img) >= 0);
else
    % don't include 0s if background is 0s
    mask = double(abs(img) > 0);
end

% set nans to 0
img(isnan(img)) = 0;

% resize img and mask
out = imresize(img,scale,method);
maskout = imresize(mask,scale,method);

% threshold mask/img
maskout(maskout < threshold) = nan;
out = out .* maskout;

end
end