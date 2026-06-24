% note : to run this code, first run setup.m in irt folder to setup paths
% for nufft toolbox; then select mapVBVD-main folder and right click to add
% path with subfolders
% voronoidens : add sparseMRI_v0.2/utils path

% VD_yjs_1H.m and spiralVD_G_scanner.m : generate the Archemidia spiral
% with optimal time

% close all
clc
clear
dropoint = 3;
MS_recon_desired = 100;            %%% CHANGED (was 90). Wavelet needs power-of-two size.

msnow = MS_recon_desired;
%
data_frames = load('ACR_test.mat');

ktrajs = data_frames.ktrajs;
kdatas = data_frames.kdatas;
kcomps = data_frames.kcomps;

%
frame = 1;
dkx = ktrajs(frame, 1, dropoint:end);
dky = ktrajs(frame, 2, dropoint:end);
dkz = ktrajs(frame, 3, dropoint:end);

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

w = kcomps(frame,dropoint:end)';
w = w/max(w(:));

data = kdatas(frame,dropoint:end)';
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
% nslcd = round(MS_recon_desired/2);
% fasel = 3;
%
% imagesc([fliplr(flipud(squeeze(img_gp(nslcd,:,:)))),rot90(squeeze(img_gp(:,nslcd+fasel-20,:))),rot90(squeeze(img_gp(:,nslcd+fasel*2-20+2,:))); ...
%     rot90(squeeze(img_gp(:,nslcd+fasel*3-20+1,:))), rot90(squeeze(img_gp(:,nslcd+fasel*4-20+2,:))),rot90(squeeze(img_gp(:,nslcd+fasel*5-20+2,:))); ...
%     rot90(squeeze(img_gp(:,nslcd+fasel*6-20+2,:))),rot90(squeeze(img_gp(:,nslcd+fasel*7-20+3,:))),rot90(squeeze(img_gp(:,nslcd+fasel*8-20+4,:))); ...
%     rot90(squeeze(img_gp(:,nslcd+fasel*9-20+4,:))),rot90(squeeze(img_gp(:,nslcd+fasel*10-20+5,:))),rot90(squeeze(img_gp(:,nslcd+fasel*11-20+6,:)))]);
%
% axis image off; colormap gray; caxis([-0.05 .7]);

% Version 240
nslcd = round(MS_recon_desired/2);
fasel = 4;
res = 2;

stepg = round(16.5/res);
shiftg = 32;

imagesc([fliplr(flipud(squeeze(img_gp(nslcd,:,:)))),rot90(squeeze(img_gp(:,nslcd-shiftg,:))),rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*1,:))); ...
    rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*2,:))), rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*3,:))),rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*4,:))); ...
    rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*5,:))),rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*6,:))),rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*7,:))); ...
    rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*8,:))),rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*9,:))),rot90(squeeze(img_gp(:,nslcd-shiftg+stepg*10,:)))]);

axis image off; colormap gray; caxis([-.1 .8]);

drawnow()

% Wavelet operator: Daubechies-4, 4 levels
% XFM = Wavelet('Daubechies',4,4);   % db4, 4 levels (per-slice wavelets in many toolboxes)

%%%% CS
XFM = 1;				% Identity transform
% TVWeight = 0.01; 	% Weight for TV penalty
TVWeight0 = 0.00; 	% Weight for TV penalty
xfmWeight = 0.00;	    % Weight for Transform L1 penalty
Itnlim = 5;		    % Number of iterations

% initialize Parameters for reconstruction
param = init;
param.FT = FT; % image to k-space transform
param.XFM = XFM;
param.TV = TVOP3D;
param.data = data;
param.TVWeight = TVWeight0;   % TV penalty (will be updated per-iteration)
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

% scale w?
% tmp=zeros(N); tmp(end/2+1,end/2+1)=1; tmp=FT'*(w.*(FT*tmp)); w = w/max(abs(tmp(:)));

% im_dc = FT'*(data.*w);	% init with zf-w/dc (zero-fill with density compensation)
% figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;               % This now works because size is 128x128 per slice.

%%% ADDED: build a smooth circular 3D taper (raised-cosine) to kill FOV-edge ringing
%%%        Centered at image center; radius 0.95 keeps most FOV, 0.05 cosine roll-off.
[Nx, Ny, Nz] = size(im_dc);
[yy, xx, zz] = ndgrid( ( (1:Ny) - (Ny+1)/2 )/((Ny-1)/2), ...
                       ( (1:Nx) - (Nx+1)/2 )/((Nx-1)/2), ...
                       ( (1:Nz) - (Nz+1)/2 )/((Nz-1)/2) );
rr = sqrt(xx.^2 + yy.^2 + zz.^2);
r0 = 0.85;                 % flat pass region (fraction of radius)
bw = 0.5;                 % taper width (cosine roll-off)
edgeWin = ones(size(rr), 'like', rr);
maskT = (rr >= r0) & (rr <= r0 + bw);
edgeWin(rr > r0 + bw) = 0;
edgeWin(maskT) = 0.5*(1 + cos(pi*(rr(maskT) - r0)/bw));   % raised-cosine to zero
% Note: multiplies image each outer loop; keeps data model intact but discourages edge energy.

%%% ADDED: continuation schedule for TV — decays every few outer loops
decayEvery = 5;      % every 5 outer loops
decayFactor = 0.25;   % multiply TV by 0.7
maxOuter = 15;

% do iterations
tic
for n=1:maxOuter
    % continuation on TV
    if n > 1 && mod(n-1, decayEvery) == 0
        param.TVWeight = param.TVWeight * decayFactor;   %%% ADDED
    end

    res = fnlCg(res,param);
    im_res = XFM'*res;

    % --- apply smooth FOV taper to suppress boundary brightening ---
    if n==1
    im_res = im_res .* edgeWin;         %%% ADDED: soft-edge projection
    end
    res    = XFM * im_res;              %%% ADDED: keep iterate consistent (XFM=1 so res=im_res)

    % image = squeeze(abs(im_res(:,:,45)));
    figure(100+n), %imshow(image,[]), drawnow
    % figure(800), %imshow(image,[]), drawnow
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
    % nslcd = round(MS_recon_desired/2);
    % fasel = 3;
    %
    % imagesc([fliplr(flipud(squeeze(img_gp(nslcd,:,:)))),rot90(squeeze(img_gp(:,nslcd+fasel-20,:))),rot90(squeeze(img_gp(:,nslcd+fasel*2-20+2,:))); ...
    %     rot90(squeeze(img_gp(:,nslcd+fasel*3-20+1,:))), rot90(squeeze(img_gp(:,nslcd+fasel*4-20+2,:))),rot90(squeeze(img_gp(:,nslcd+fasel*5-20+2,:))); ...
    %     rot90(squeeze(img_gp(:,nslcd+fasel*6-20+2,:))),rot90(squeeze(img_gp(:,nslcd+fasel*7-20+3,:))),rot90(squeeze(img_gp(:,nslcd+fasel*8-20+4,:))); ...
    %     rot90(squeeze(img_gp(:,nslcd+fasel*9-20+4,:))),rot90(squeeze(img_gp(:,nslcd+fasel*10-20+5,:))),rot90(squeeze(img_gp(:,nslcd+fasel*11-20+6,:)))]);
    %
    % axis image off; colormap gray; caxis([-0.05 .7]);

    % Version 240
    % Version 240
    nslcd = round(MS_recon_desired/2);
    fasel = 4;
    resul = 2;

    stepg = round(16.5/resul);
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
save wlet_ph1_adopt gas

toc