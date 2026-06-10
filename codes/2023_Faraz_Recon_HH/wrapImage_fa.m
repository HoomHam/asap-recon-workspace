function img_large = wrapImage_fa(img, rows, cols, Options)

% parse arguments
arguments 
   img 
   rows double
   cols double
   Options.hgap double = 10
   Options.vgap double = 10
   Options.scale double = 1
end

% get dimensions
dim = size(img);
if length(dim) == 3
    dim(4) = 1;
end

% check inputs and assign defaults
if nargin < 2, rows = []; end
if nargin < 3, cols = []; end
if isempty(rows), rows = ceil(dim(3)/3); end
if isempty(cols), cols = 3; end

% unpack options
hgap = Options.hgap;
vgap = Options.vgap;
scale = Options.scale;

% make gaps even
if mod(hgap,2) == 1
    hgap = ceil(hgap + 1);
end
if mod(vgap,2) == 1
    vgap = ceil(vgap + 1);
end

% initialize
sizex = dim(1)*scale;
sizey = dim(2)*scale;
img_large = zeros((sizex+vgap)*rows,(sizey+hgap)*cols,dim(4));

% pad and insert slices into img_large
for m = 1:dim(4)
    k = 1;
    for i = 1:rows
        for j = 1:cols
            if k <= dim(3)
                if hgap>=0 && vgap >=0
                    img_large(1+(i-1)*(sizex+vgap):i*(sizex+vgap),1+(j-1)*(sizey+hgap):j*(sizey+hgap),m) = padarray(resize_fa(img(:,:,k,m),scale,'bilinear'),abs([vgap/2 hgap/2]),0,'both');
                else
                    tmp = resize_fa(abs(img(:,:,k,m)),scale,'bilinear');
                    if hgap < 0 && vgap >=0
                        img_large(1+(i-1)*(sizex+vgap):i*(sizex+vgap),1+(j-1)*(sizey+hgap):j*(sizey+hgap),m) = tmp(:,abs(hgap/2)+1:end-abs(hgap/2));
                    elseif hgap < 0 && vgap < 0
                        img_large(1+(i-1)*(sizex+vgap):i*(sizex+vgap),1+(j-1)*(sizey+hgap):j*(sizey+hgap),m) = tmp(abs(vgap/2)+1:end-abs(vgap/2),abs(hgap/2)+1:end-abs(hgap/2));
                    elseif hgap >= 0 && vgap < 0
                        img_large(1+(i-1)*(sizex+vgap):i*(sizex+vgap),1+(j-1)*(sizey+hgap):j*(sizey+hgap),m) = tmp(abs(vgap/2)+1:end-abs(vgap/2),:);
                    end
                end
                k = k + 1;
            end
        end
    end
end
end
