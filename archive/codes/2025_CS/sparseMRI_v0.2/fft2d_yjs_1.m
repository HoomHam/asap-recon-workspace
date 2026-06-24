close all
clc
clear

% generate image
img = phantom(256);
figure
imshow(img);

% full k-space
k = fftshift(fft2(img));
k_mag = log(abs(k));
figure
imshow(k_mag, []);

% undersampling: skip every other phase line
k_us = zeros(256,256);
k_us(1:2:end, :) = k(1:2:end, :);
img_us = ifft2(k_us);
figure
imshow(abs(img_us), [])

% variable density random sampling in phase encoding
DN = size(img);
P = 5;			% Variable density polymonial degree
pctg = [0.5];  	% undersampling factor
pdf = genPDF(256, P, pctg , 2 ,0.1, 1);	% generates the sampling PDF
rand_sampling = genSampling(pdf,10,60);		% generates a sampling pattern

% 
k_rand = k.*rand_sampling; % mulplication broadcast; 
figure
k_rand_mag = log(abs(k_rand));
imshow(k_rand_mag,[])

img_rand = ifft2(k_rand);
figure
imshow(abs(img_rand), [])

% CS
rand_sampling2d = zeros(size(k));
for i = 1:256
    rand_sampling2d(i,:) = rand_sampling;
end

FT = p2DFT(rand_sampling2d, size(k), 1, 2); % phase: 1

data = FT*img;

XFM = 1;				% Identity transform 	
TVWeight = 0.01; 	% Weight for TV penalty
xfmWeight = 0.00;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

res = XFM*img_rand;

% do iterations
tic
for n=1:8
	res = fnlCg(res,param);
	im_res = XFM'*res;
	figure(100+n), imshow(abs(im_res),[]), drawnow
end
toc

