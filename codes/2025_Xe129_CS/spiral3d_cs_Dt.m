% note : to run this code, first run setup.m in irt folder to setup paths
% for nufft toolbox; then select mapVBVD-main folder and right click to add
% path with subfolders
% voronoidens : add sparseMRI_v0.2/utils path

% VD_yjs_1H.m and spiralVD_G_scanner.m : generate the Archemidia spiral
% with optimal time

close all
clc
clear

%
data_frames = load('xe129_frames.mat');

ktrajs = data_frames.ktrajs;
kdatas = data_frames.kdatas;
kcomps = data_frames.kcomps;

% generate circular mask (spirals hav a circular FOV support
% [xx,yy] = meshgrid(linspace(-1,1,N));
% ph = double(sqrt(xx.^2 + yy.^2)<1);
ph = 1;

slices = 128;
frames = 16;
imSize = [256, slices, 256, frames];
FT = NUFFT4D(ktrajs, 1, ph, 0, imSize, 2); % NUFFT1Dt = NUFFT4D

% scale and density compensation
for fr=1:frames
    kdatas(fr,:) = kdatas(fr,:) / max(abs(kdatas(fr,:))); % scale kdata, NOT scaleing initial image
end

im_dc = FT'*(kdatas.*kcomps);

% 
% for sl=1:slices
%     for fr=1:frames
%         figure;
%         imshow(flipud(squeeze(abs(im_dc(:,sl,:,fr)))'), []);
%         title(strcat('slice', int2str(sl), '-frame', int2str(fr)));
%         path = strcat('C:\Users\P53-LOCAL\Desktop\AI\compressed_sensing\xe129\D4\regular\slice', int2str(sl));
%         mkdir(path);
%         if ~exist(path, 'dir')
%             mkdir(path)
%         end
%         saveas(gcf,strcat(path,'\slice',int2str(sl), '-frame', int2str(fr), '.png'));
%     end
%     close all;
% end


%%%% CS
XFM = 1;				% Identity transform
TVWeight = 0.01; 	% Weight for TV penalty
xfmWeight = 0.0;	% Weight for Transform L1 penalty
Itnlim = 25;		% Number of iterations

% initialize Parameters for reconstruction
param = init;
param.FT = FT; % image to k-space transform
param.XFM = XFM;
param.TV = TVOPDt;
param.data = kdatas;
param.TVWeight =TVWeight;     % TV penalty
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

% 
res = XFM*im_dc;

clear im_dc;

% do iterations
tic
for n=1:3
    res = fnlCg(res,param);
    im_res = XFM'*res;
    image = squeeze(abs(im_res(:,55,:,7)));
    figure(100+n), imshow(image,[]), drawnow
end
toc


for sl=1:slices
    for fr=1:frames
        figure;
        imshow(flipud(squeeze(abs(im_res(:,sl,:,fr)))'), []);
        title(strcat('CS-slice', int2str(sl), '-frame', int2str(fr)));
        path = strcat('C:\Users\P53-LOCAL\Desktop\AI\compressed_sensing\xe129\Dt\CS\slice', int2str(sl));
        mkdir(path);
        if ~exist(path, 'dir')
            mkdir(path)
        end
        saveas(gcf,strcat(path,'\CS-slice',int2str(sl), '-frame', int2str(fr), '.png'));
    end
    close all;
end

save('C:\Users\P53-LOCAL\Desktop\AI\compressed_sensing\xe129\spiral3d_CSDt.mat', 'im_res', '-v7.3');