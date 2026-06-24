% note : to run this code, first run setup.m in irt folder to setup paths
% for nufft toolbox; then select mapVBVD-main folder and right click to add
% path with subfolders
% voronoidens : add sparseMRI_v0.2/utils path

% VD_yjs_1H.m and spiralVD_G_scanner.m : generate the Archemidia spiral
% with optimal time

% close all
clc
clear
MS_recon_desired = 100;

msnow = MS_recon_desired;
%
data_frames = load('ACR_test.mat');

ktrajs = data_frames.ktrajs;
kdatas = data_frames.kdatas;
kcomps = data_frames.kcomps;

%
frame = 1;
dkx = ktrajs(frame, 1, :);
dky = ktrajs(frame, 2, :);
dkz = ktrajs(frame, 3, :);

% [-0.5, 0.5]
dkx = dkx(:)/(2*pi);
dky = dky(:)/(2*pi);
dkz = dkz(:)/(2*pi);

k = [dkx, dky, dkz];

% generate circular mask (spirals hav a circular FOV support
% [xx,yy] = meshgrid(linspace(-1,1,N));
% ph = double(sqrt(xx.^2 + yy.^2)<1);
ph = 1;

slices = msnow;
imSize = [msnow, slices, msnow];
FT = NUFFT3D(k, 1, ph, 0, imSize, 2);

% w = voronoidens(k); % calculate voronoi density compensation function
% w = w/max(w(:));

w = kcomps(frame,:)';
w = w/max(w(:));

data = kdatas(frame,:)';
data = data / max(abs(data)); % scale kdata, NOT scaleing initial image

im_dc = FT'*(data.*w);

figure(1);
imshow(squeeze(abs(im_dc(:,:,45))),[])

img_gp = abs(im_dc);

for jj=1:size(img_gp,1)
    tmp = squeeze(img_gp(jj,:,:));
    img_gp(jj,:,:) = rot90(tmp);
end

img_gp = img_gp/max(img_gp(:));

% Version 2
nslcd = round(MS_recon_desired/2);
fasel = 3;

imagesc([fliplr(flipud(squeeze(img_gp(nslcd,:,:)))),rot90(squeeze(img_gp(:,nslcd+fasel-20,:))),rot90(squeeze(img_gp(:,nslcd+fasel*2-20+2,:))); ...
    rot90(squeeze(img_gp(:,nslcd+fasel*3-20+1,:))), rot90(squeeze(img_gp(:,nslcd+fasel*4-20+2,:))),rot90(squeeze(img_gp(:,nslcd+fasel*5-20+2,:))); ...
    rot90(squeeze(img_gp(:,nslcd+fasel*6-20+2,:))),rot90(squeeze(img_gp(:,nslcd+fasel*7-20+3,:))),rot90(squeeze(img_gp(:,nslcd+fasel*8-20+4,:))); ...
    rot90(squeeze(img_gp(:,nslcd+fasel*9-20+4,:))),rot90(squeeze(img_gp(:,nslcd+fasel*10-20+5,:))),rot90(squeeze(img_gp(:,nslcd+fasel*11-20+6,:)))]);

axis image off; colormap gray; %caxis([-2 22]);

drawnow()


% Wavelet operator: Daubechies-4, 4 levels
% XFM = Wavelet('Daubechies',4,4);   % db4, 4 levels (per-slice wavelets in many toolboxes)

%%%% CS
XFM = 1;				% Identity transform
TVWeight = 0.01; 	% Weight for TV penalty
xfmWeight = 0.01;	% Weight for Transform L1 penalty
Itnlim = 15;		% Number of iterations

% initialize Parameters for reconstruction
param = init;
param.FT = FT; % image to k-space transform
param.XFM = XFM;
param.TV = TVOP3D;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

% scale w?
% tmp=zeros(N); tmp(end/2+1,end/2+1)=1; tmp=FT'*(w.*(FT*tmp)); w = w/max(abs(tmp(:)));

% im_dc = FT'*(data.*w);	% init with zf-w/dc (zero-fill with density compensation)
% figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
tic
for n=1:15
    res = fnlCg(res,param);
    im_res = XFM'*res;
    % image = squeeze(abs(im_res(:,:,45)));
    figure(100+n), %imshow(image,[]), drawnow
    img_gp = abs(im_res);

    for jj=1:size(img_gp,1)
        tmp = squeeze(img_gp(jj,:,:));
        img_gp(jj,:,:) = rot90(tmp);
    end

    img_gp = img_gp/max(img_gp(:));

    %     for jj=1:size(img_gp,1)
    % tmp = squeeze(img_gp(jj,:,:));
    % img_gp(jj,:,:) = rot90(tmp);
    % end

    % Version 2
    nslcd = round(MS_recon_desired/2);
    fasel = 3;

    imagesc([fliplr(flipud(squeeze(img_gp(nslcd,:,:)))),rot90(squeeze(img_gp(:,nslcd+fasel-20,:))),rot90(squeeze(img_gp(:,nslcd+fasel*2-20+2,:))); ...
        rot90(squeeze(img_gp(:,nslcd+fasel*3-20+1,:))), rot90(squeeze(img_gp(:,nslcd+fasel*4-20+2,:))),rot90(squeeze(img_gp(:,nslcd+fasel*5-20+2,:))); ...
        rot90(squeeze(img_gp(:,nslcd+fasel*6-20+2,:))),rot90(squeeze(img_gp(:,nslcd+fasel*7-20+3,:))),rot90(squeeze(img_gp(:,nslcd+fasel*8-20+4,:))); ...
        rot90(squeeze(img_gp(:,nslcd+fasel*9-20+4,:))),rot90(squeeze(img_gp(:,nslcd+fasel*10-20+5,:))),rot90(squeeze(img_gp(:,nslcd+fasel*11-20+6,:)))]);

    axis image off; colormap gray; %caxis([-5 75]);
    drawnow()

    gas(n,:,:,:) = img_gp;

    %     set(gcf,'units','normalized','outerposition',[0 0 1 1])
end

save tv01g01 gas
toc