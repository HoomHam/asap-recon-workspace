% This script demonstare a CS reconstruction from 
% Randomly undersmapled phase encodes of 2D FSE
% of a brain image.

close all
clc
clear

addpath(strcat(pwd,'/utils'));

if exist('FWT2_PO') <2
	error('must have Wavelab installed and in the path');
end

load brain512

img = fftshift(fft2(data));
img_mag = abs(img);
img_mag = img_mag/max(img_mag(:));
img_mag = flipud(fliplr(img_mag));
figure("Name","orginal image");
imshow(img_mag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(data); 	% image Size
DN = size(data); 	% data Size
TVWeight = 0.002; 	% Weight for TV penalty
xfmWeight = 0.005;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations


%generate Fourier sampling operator
FT = p2DFT(mask, N, 1, 2);

% scale data
im_dc = FT'*(data.*mask./pdf);
data = data/max(abs(im_dc(:)));
im_dc = im_dc/max(abs(im_dc(:)));

%generate transform operator
XFM = Wavelet('Daubechies',4,4);	% Wavelet

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
tic
for n=1:5
	res = fnlCg(res,param);
	im_res = XFM'*res;
	figure(100), imshow(abs(im_res),[]), drawnow
end
toc


figure, imshow(abs(cat(2,im_dc,im_res)),[]);
figure, imshow(abs(cat(2,im_dc(155:304,110:209), im_res(155:304,110:209))),[0,1],'InitialMagnification',200);
title(' zf-w/dc              l_1 Wavelet');

figure, imshow(abs(cat(2, img_mag(155:304,110:209), im_dc(155:304,110:209), im_res(155:304,110:209))),[0,1],'InitialMagnification',200);



