% note : to run this code, first run setup.m in irt folder to setup paths
% for nufft toolbox; then select mapVBVD-main folder and right click to add
% path with subfolders
% voronoidens : add sparseMRI_v0.2/utils path

% VD_yjs_1H.m and spiralVD_G_scanner.m : generate the Archemidia spiral
% with optimal time

% close all
clc
clear
fov=350;
dropoint = 1;
dropointst = 3;
dropointen = 5;

lengul = 16;
shiftg = 32;
reptur = 10;

MS_recon_desired = 128;            %%% CHANGED (was 90). Wavelet needs power-of-two size.
resul = fov/MS_recon_desired;

msnow = MS_recon_desired;
%
data_frames = load('ACR_test.mat');

ktrajs = data_frames.ktrajs;
kdatas = data_frames.kdatas;
kcomps = data_frames.kcomps;

%
frame = 1;
dkx = ktrajs(frame, 1, dropointst:end-dropointen);
dky = ktrajs(frame, 2, dropointst:end-dropointen);
dkz = ktrajs(frame, 3, dropointst:end-dropointen);

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

w = kcomps(frame,dropointst:end-dropointen)';
w = w/max(w(:));

data = kdatas(frame,dropointst:end-dropointen)';
% data = data / max(abs(data)); % scale kdata, NOT scaleing initial image

im_dc = FT'*(data.*w);

figure(1);
imshow(squeeze(abs(im_dc(:,:,45))),[])

img_gp = abs(im_dc);

for jj=1:size(img_gp,1)
    tmp = squeeze(img_gp(jj,:,:));
    img_gp(jj,:,:) = rot90(tmp);
end

% img_gp = img_gp/max(img_gp(:));

% Version 2
nslcd = round(MS_recon_desired/2);
fasel = 3;

imagesc([fliplr(flipud(squeeze(img_gp(nslcd,:,:)))),rot90(squeeze(img_gp(:,nslcd+fasel-20,:))),rot90(squeeze(img_gp(:,nslcd+fasel*2-20+2,:))); ...
    rot90(squeeze(img_gp(:,nslcd+fasel*3-20+1,:))), rot90(squeeze(img_gp(:,nslcd+fasel*4-20+2,:))),rot90(squeeze(img_gp(:,nslcd+fasel*5-20+2,:))); ...
    rot90(squeeze(img_gp(:,nslcd+fasel*6-20+2,:))),rot90(squeeze(img_gp(:,nslcd+fasel*7-20+3,:))),rot90(squeeze(img_gp(:,nslcd+fasel*8-20+4,:))); ...
    rot90(squeeze(img_gp(:,nslcd+fasel*9-20+4,:))),rot90(squeeze(img_gp(:,nslcd+fasel*10-20+5,:))),rot90(squeeze(img_gp(:,nslcd+fasel*11-20+6,:)))]);
%
% axis image off; colormap gray; caxis([-0.05 .7]);

% Version 240
% nslcd = round(MS_recon_desired/2);
% 
% 
% stepg = round(lengul/resul);
% 
% imagesc([fliplr(flipud(squeeze(img_gp(nslcd,:,:)))),rot90(squeeze(img_gp(:,nslcd-shiftg,:))),rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*1,:))); ...
%     rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*2,:))), rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*3,:))),rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*4,:))); ...
%     rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*5,:))),rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*6,:))),rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*7,:))); ...
%     rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*8,:))),rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*9,:))),rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*10,:)))]);

axis image off; colormap gray; %caxis([-.1 .8]);

drawnow()

% Wavelet operator: Daubechies-4, 4 levels
XFM = Wavelet('Daubechies',4,4);   % db4, 4 levels (per-slice wavelets in many toolboxes)

%%%% CS
% XFM = 1;				% Identity transform
TVWeight = 0.00; 	% Weight for TV penalty
xfmWeight = 0.01;	% Weight for Transform L1 penalty
Itnlim = 200;		% Number of iterations

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

res = XFM*im_dc;               % This now works because size is 128x128 per slice.

% do iterations
tic
for n=1:reptur
    res = fnlCg(res,param);
    im_res = XFM'*res;
    % image = squeeze(abs(im_res(:,:,45)));
    figure(500+n), %imshow(image,[]), drawnow
    % figure(800), %imshow(image,[]), drawnow
    img_gp = abs(im_res);

    for jj=1:size(img_gp,1)
        tmp = squeeze(img_gp(jj,:,:));
        img_gp(jj,:,:) = rot90(tmp);
    end

    % img_gp = img_gp/max(img_gp(:));

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
    %
    % axis image off; colormap gray; caxis([-0.05 .7]);

    % Version 240
    % Version 240
    % nslcd = round(MS_recon_desired/2);
    % 
    % stepg = round(lengul/resul);
    % 
    % imagesc([fliplr(flipud(squeeze(img_gp(nslcd,dropoint:end-dropointen,dropoint:end-dropointen)))),rot90(squeeze(img_gp(dropoint:end-dropointen,nslcd-shiftg,dropoint:end-dropointen))),rot90(squeeze(img_gp(dropoint:end-dropointen,nslcd-shiftg+stepg*1,dropoint:end-dropointen))); ...
    %     rot90(squeeze(img_gp(dropoint:end-dropointen,nslcd-shiftg+stepg*2,dropoint:end-dropointen))), rot90(squeeze(img_gp(dropoint:end-dropointen,nslcd-shiftg+stepg*3,dropoint:end-dropointen))),rot90(squeeze(img_gp(dropoint:end-dropointen,nslcd-shiftg+stepg*4,dropoint:end-dropointen))); ...
    %     rot90(squeeze(img_gp(dropoint:end-dropointen,nslcd-shiftg+stepg*5,dropoint:end-dropointen))),rot90(squeeze(img_gp(dropoint:end-dropointen,nslcd-shiftg+stepg*6,dropoint:end-dropointen))),rot90(squeeze(img_gp(dropoint:end-dropointen,nslcd-shiftg+stepg*7,dropoint:end-dropointen))); ...
    %     rot90(squeeze(img_gp(dropoint:end-dropointen,nslcd-shiftg+stepg*8,dropoint:end-dropointen))),rot90(squeeze(img_gp(dropoint:end-dropointen,nslcd-shiftg+stepg*9,dropoint:end-dropointen))),rot90(squeeze(img_gp(dropoint:end-dropointen,nslcd-shiftg+stepg*10,dropoint:end-dropointen)))]);

    axis image off; colormap gray; %caxis([-.1 .8]);

    drawnow()

    gas(n,:,:,:) = img_gp;
    save wlet_ph1 gas


    %     set(gcf,'units','normalized','outerposition',[0 0 1 1])
end
save wlet_ph1 gas


toc

% img_gp = img_gp/max(img_gp(:));
% dropoint = 5;
% 
% figure(500+n+1), %imshow(image,[]), drawnow
% 
% nslcd = round(MS_recon_desired/2);
% 
% stepg = round(lengul/resul);
% 
% imagesc([fliplr(flipud(squeeze(img_gp(nslcd,dropoint:end-dropoint,dropoint:end-dropoint)))),rot90(squeeze(img_gp(dropoint:end-dropoint,nslcd-shiftg,dropoint:end-dropoint))),rot90(squeeze(img_gp(dropoint:end-dropoint,nslcd-shiftg+stepg*1,dropoint:end-dropoint))); ...
%     rot90(squeeze(img_gp(dropoint:end-dropoint,nslcd-shiftg+stepg*2,dropoint:end-dropoint))), rot90(squeeze(img_gp(dropoint:end-dropoint,nslcd-shiftg+stepg*3,dropoint:end-dropoint))),rot90(squeeze(img_gp(dropoint:end-dropoint,nslcd-shiftg+stepg*4,dropoint:end-dropoint))); ...
%     rot90(squeeze(img_gp(dropoint:end-dropoint,nslcd-shiftg+stepg*5,dropoint:end-dropoint))),rot90(squeeze(img_gp(dropoint:end-dropoint,nslcd-shiftg+stepg*6,dropoint:end-dropoint))),rot90(squeeze(img_gp(dropoint:end-dropoint,nslcd-shiftg+stepg*7,dropoint:end-dropoint))); ...
%     rot90(squeeze(img_gp(dropoint:end-dropoint,nslcd-shiftg+stepg*8,dropoint:end-dropoint))),rot90(squeeze(img_gp(dropoint:end-dropoint,nslcd-shiftg+stepg*9,dropoint:end-dropoint))),rot90(squeeze(img_gp(dropoint:end-dropoint,nslcd-shiftg+stepg*10,dropoint:end-dropoint)))]);
% 
% axis image off; colormap gray; caxis([-0.025 .325]);
% 
% drawnow()