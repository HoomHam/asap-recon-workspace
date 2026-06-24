% this is a script to demonstrate a TV recon from undersampled variable density spirals
%
% (c) Michael Lustig 2007

% note : to run this code, run setup.m in irt folder to setup paths
% for nufft toolbox; 
close all
clear
clc
addpath(strcat(pwd,'/utils'));

% add nufft toolbox path : C:\Users\P53-LOCAL\Desktop\AI\compressed_sensing\codes\irt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = [160,160]; 		% image Size
TVWeight = 0.01; 	% Weight for TV penalty
xfmWeight = 0.00;	% Weight for Transform L1 penalty
Itnlim = 25;		% Number of iterations

% generate a variable density spiral with gaussian density.
SIGMA = 3;
FOV = gausswin(40,SIGMA);
FOV = FOV(end/2+1:end);
FOV = FOV*10+7;
RADIUS = linspace(0,1,20);
RESOLUTION = 1; % in mm, such that image is 256x256
NITLV = 16; % number of spiral interleaves: 16
Gmax = 2 ; % maximum gradient in [G/CM];
Smax = 10; %Maximum slew-rate; [G/(CM*ms)]
T = 4e-3; % time sampling (in mS);

disp('design spiral')
gamma = 4.257;
[k,g,s,time] = vdSpiralDesign_old(NITLV, RESOLUTION,FOV,RADIUS,gamma,Gmax,Smax,T,'cubic');
k = k(2:end).'*exp(2*pi*i*[1:NITLV]/NITLV);
k = k(:)/max(abs(k(:)))/2; % scale to range [-0.5,0.5]

% plot k-trajectory in (kx. ky);
readout = length(k)/NITLV;
figure('Name', "k-trajectory (kx, ky)");
for shot=1:1
    index = ((shot-1)*readout+1):shot*readout; 
    plot(real(k(index)), imag(k(index)));
    hold on
end
axis equal
grid

% plot k-trajectory in time domain;
t = linspace(0,time,readout);
figure('Name', "K-trajectory: time domain");
for shot=1:1
    index = ((shot-1)*readout+1):shot*readout; 
    plot(t, real(k(index)));
    hold on
    plot(t, imag(k(index)));
    hold on
end
grid

% plot gradients in (gx, gy);
g = g*10; % G/cm -> mT/m;
figure('Name', "g-trajectory (gx, gy)");
for shot=1:1
    index = ((shot-1)*readout+1):shot*readout; 
    plot(real(g(index)), imag(g(index)));
    hold on
end
axis equal
grid

% plot gradients in time domain;
figure('Name', "K-trajectory : time domain");
for shot=1:1
    index = ((shot-1)*readout+1):shot*readout; 
    plot(t, real(g(index)));
    hold on
    plot(t, imag(g(index)));
    hold on
end
grid
xlabel('Time (ms)');
ylabel('Gradient (mT/m)');

% plot slewrate in (sx, sy);
s = s*10; % mT/(m*ms)
figure('Name', "s-trajectory (sx, sy)");
for shot=1:1
    index = ((shot-1)*readout+1):shot*readout; 
    plot(real(s(index)), imag(s(index)));
    hold on
end
axis equal
grid

% plot slewrate in time domain;
figure('Name', "s-trajectory : time domain");
for shot=1:1
    index = ((shot-1)*readout+1):shot*readout; 
    plot(t, real(s(index)));
    hold on
    plot(t, imag(s(index)));
    hold on
end
grid
xlabel('Time (ms)');
ylabel('Slew rate (mT/(m*ms))');


w = voronoidens(k); % calculate voronoi density compensation function
w = w/max(w(:));

% generate circular mask (spirals hav a circular FOV support
[xx,yy] = meshgrid(linspace(-1,1,N(1)));
ph = double(sqrt(xx.^2 + yy.^2)<1);

%generate image
im = (phantom(N(1)))  + randn(N)*0.01 + i*randn(N)*0.01;

%generate Fourier sampling operator
FT = NUFFT(k,1, ph, 0,N, 2);

% scale w
tmp=zeros(N); 
tmp(end/2+1,end/2+1)=1; 
tmp=FT'*(w.*(FT*tmp)); 
w = w/max(abs(tmp(:)));

data = FT*im;

%generate transform operator

%XFM = Wavelet('Daubechies',6,4);	% Wavelet
%XFM = TIDCT(8,4);			% DCT
XFM = 1;				% Identity transform 	

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

im_dc = FT'*(data.*w);	% init with zf-w/dc (zero-fill with density compensation)
figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
tic
for n=1:8
	res = fnlCg(res,param);
	im_res = XFM'*res;
% 	figure(100+n), imshow(abs(im_res),[]), drawnow
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
end
toc


im_full = ifft2c(zpad(fft2c(im),N(1),N(2)));
figure, imshow(abs(cat(2,im_full,im_dc,im_res)),[0,1]);
title('original                      zf-w/dc              TV');

figure, plot(1:N(1), abs(im_full(end/2,:)), 1:N(2), abs(im_dc(end/2,:)), 1:N(2), abs(im_res(end/2,:)),'LineWidth',2);
legend('original', 'zf-w/dc', 'TV');


