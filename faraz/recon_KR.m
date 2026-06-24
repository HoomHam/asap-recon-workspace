% spiral_human_template.m
% Subject:  
% Methods:
%   Free-breathing
%   UVA coil
%   
% Notes: 
%
%
%
% clc;

clear all;
close all;
% addpath("C:\Users\kevin\Box\Human_MRI_LTX\Lib\Faraz Recon");
% addpath("C:\Users\kevin\Box\Human_MRI_LTX\Lib\mapVBVD");
% addpath("C:\Users\kevin\Box\Human_MRI_LTX\Lib\Faraz Recon\fsroot");
% addpath("C:\Users\kevin\Box\Human_MRI_LTX\Lib\extra_code");

%% Control parameters
GradRasterTime_us = 10;    % Gradient raster time

%% Load Images
tic;
files = [
    % Free-breathing
"meas_MID00305_FID03483_fa_spiral_dyn_fancy_v2_20230131.dat"
];

[data,twixs] = recon_20210622(files,'nStudies',1:length(files),'recon',0);


%% Add to log
logfile = "C:\Users\kevin\Box\Human MRI LTX\patient_log.xlsx"; % File path for log
log = cell(length(data),10);
for i = 1:length(data)
    log{i,1} = 20230426; % Current Date (as a number: YYYYMMDD)
    log{i,2} = '023DB'; % Subject ID
    log{i,3} = 'LTX'; % Subject Cohort
    log{i,4} = data(i).protocol;
    if contains(data(i).protocol,'hybrid')
        log{i,5} = 'MB';
        log{i,6} = 'EE';
    else
        if data(i).numImg < 400 % assume SB if less than 400 images
            log{i,5} = 'SB';
            log{i,6} = 'EI';
        else
            log{i,5} = 'FB';
            log{i,6} = 'NA';
        end
    end
    log{i,7} = data(i).frequency;
    log{i,8} = data(i).voltage;
    log{i,9} = ''; % misc notes
    log{i,10} = char(files(i));
    
    data(i).date = datestr(datenum(num2str(log{i,1}),'yyyymmdd'),'yyyy-mm-dd'); 
    data(i).cohort = char(log{i,3}); 
    data(i).notes = char(log{i,9});
    data(i).patient = char(log{i,2}); 
    
end
% writecell(log,logfile,'WriteMode','append'); clearvars log logfile;

%% Get fancy spiral data 
i = 1;
nsamples = data(i).nCol;
nleaves = data(i).nLin;
fov = data(i).fovPE;
matsize = data(i).MatSize;
nreps = data(i).GAperiod;
nrepetitions = data(i).nRep;
nchannels = data(i).nCha;
nres = data(i).numRes;
if data(i).bSpectra
    nspec = data(i).numSpec;
end
if strcmpi(data(i).acqOrder,'GP-DP-DP') && data(i).numRes > 1
    nres = 3;
end

KSpaceCoor = loadtrajectory3D(data(i).frequency/1.498/1e6,fov,nsamples,nleaves,nreps,matsize,'axial',[],[],'calibrations_3D_20220308.mat');
if nres > 1
    KSpaceCoor_DP = loadtrajectory3D(data(i).frequency/1.498/1e6,fov,nsamples,nleaves,nreps,matsize,'axial','dissolved',[],'calibrations_3D_20220308.mat');
end

fnames = fieldnames(twixs);
rawdata = squeeze(twixs.(fnames{i}).image(:,:,:,:,:,:,:,:,:,:,:));
rawdata = reshape(rawdata,nsamples,nchannels,[]);

if data(i).bSpectra
    rawspec = rawdata(:,:,1:nspec);
    rawdata = rawdata(:,:,nspec+1:end);
end

if isfield(twixs,'file2')
    rawdata_sb = squeeze(twixs.file2.image(:,:,:,:,:,:,:,:,:,:,:));
    rawdata_sb = reshape(rawdata_sb,nsamples,nchannels,[]);
end

intind = [];
dpind = [];
for j = 1:ceil(size(rawdata,3)/nres)
   intind = [intind; repmat(j,nres,1)]; 
   dpind = [dpind; (0:nres-1)']; 
end
intind(size(rawdata,3)+1:end) = [];
dpind(size(rawdata,3)+1:end) = [];

absval = reshape(squeeze(abs(sum(rawdata(1,:,:),2))),[],1);
figure; 
plot(absval);

if exist('rawdata_sb','var')
    intind_sb = [];
    dpind_sb = [];
    for j = 1:ceil(size(rawdata_sb,3)/nres)
        intind_sb = [intind_sb; repmat(j,nres,1)];
        dpind_sb = [dpind_sb; (0:nres-1)'];
    end
    intind_sb(size(rawdata_sb,3)+1:end) = [];
    dpind_sb(size(rawdata_sb,3)+1:end) = [];

    absval_sb = reshape(squeeze(abs(sum(rawdata_sb(1,:,:),2))),[],1);
    figure;
    plot(absval_sb);
end

%% Filter Noise Spikes

med_noise = median(mean(abs(rawdata(nsamples-10:nsamples,:,:)),1));
rawnorm = abs(rawdata) ./ med_noise;

if exist('rawdata_sb','var')
    med_noise_sb = median(mean(abs(rawdata_sb(nsamples-10:nsamples,:,:)),1));
    rawnorm_sb = abs(rawdata_sb) ./ med_noise_sb;
end

for i = 0:nres-1
    movavg_norm = rawnorm(:,:,dpind==i) ./ movmean(rawnorm(:,:,dpind==i),3,3);
    thresh = 1 + 3 * std(movavg_norm,[],3,'omitnan');
    bad_inds = find(movavg_norm > thresh);
    disp(['Removed ',num2str(numel(bad_inds)),' (',num2str(numel(bad_inds)/numel(movavg_norm)*100,2),'%) noisy data points']);
    raw_tmp = rawdata(:,:,dpind==i);
    raw_tmp(bad_inds) = nan;
    for j = 1:length(bad_inds)
        [x,y,z] = ind2sub(size(raw_tmp),bad_inds(j));
        raw_tmp(x,y,z) = mean([raw_tmp(max(1,x-1),y,z); 
                               raw_tmp(min(nsamples,x+1),y,z);
                               raw_tmp(x,y,max(1,z)); 
                               raw_tmp(x,y,min(size(raw_tmp,3),z))],'omitnan');
    end
    rawdata(:,:,dpind==i) = raw_tmp;
end

absval = reshape(squeeze(abs(sum(rawdata(1,:,:),2))),[],1);
figure; 
plot(absval);

if exist('rawdata_sb','var')
    for i = 0:nres-1
        movavg_norm = rawnorm_sb(:,:,dpind_sb==i) ./ movmean(rawnorm_sb(:,:,dpind_sb==i),3,3);
        thresh = 1 + 3 * std(movavg_norm,[],3,'omitnan');
        bad_inds = find(movavg_norm > thresh);
        disp(['Removed ',num2str(numel(bad_inds)),' (',num2str(numel(bad_inds)/numel(movavg_norm)*100,2),'%) noisy data points']);
        raw_tmp = rawdata_sb(:,:,dpind_sb==i);
        raw_tmp(bad_inds) = nan;
        for j = 1:length(bad_inds)
            [x,y,z] = ind2sub(size(raw_tmp),bad_inds(j));
            raw_tmp(x,y,z) = mean([raw_tmp(max(1,x-1),y,z);
                raw_tmp(min(nsamples,x+1),y,z);
                raw_tmp(x,y,max(1,z));
                raw_tmp(x,y,min(size(raw_tmp,3),z))],'omitnan');
        end
        rawdata_sb(:,:,dpind_sb==i) = raw_tmp;
    end

    absval_sb = reshape(squeeze(abs(sum(rawdata_sb(1,:,:),2))),[],1);
    figure;
    plot(absval_sb);
end

%% Find peaks (NEEDS TO BE MANUALLY/INDIVIDUALLY EVALUATED)
range = 859:21349;
%range = 1:length(absval);
threshold = 2.5e-3;
findpeaks(absval(range(1):1:range(end)),'MinPeakDistance',15*nleaves,'MinPeakHeight',threshold);
[~,loc] = findpeaks(absval(range(1):1:range(end)),'MinPeakDistance',20*nleaves,'MinPeakHeight',threshold); % Change the last as need
loc = loc + range(1) - 1;
scantime = (loc(end)-loc(1))*(data(1).TR/nres)/1000;

disp([num2str(length(loc)),' breaths over ',num2str(scantime),' seconds (~',num2str(length(loc)/scantime*60),' br/min)']);

%% Bin based on interleaves
nphases = 16;
binned = cell(nphases,1);
for i = 1:length(loc)-1
    pts = linspace(loc(i),loc(i+1)-1,nphases+1);
    pts = pts(2:end);
    for j = loc(i):loc(i+1)-1
        ind = find(j <= pts, 1);
        binned{ind} = [binned{ind} j];
    end
end

cmap = hsv(nphases);
figure; 
plot(absval(1:1:end),'k','LineWidth',2); hold on;
for i = 1:nphases
    scatter(binned{i}(:),absval(binned{i}(:)),15,cmap(i,:),'filled');
end
hold off; 
axis off;

%% Combine binned kspace data
rawdata2 = zeros(nsamples,nleaves,nreps,nchannels,nphases,'single');
weights = zeros(nsamples,nleaves,nreps,nchannels,nphases,'single');
if nres > 1
    rawdata2_dp = zeros(nsamples,nleaves,nreps,nchannels,nphases,'single');
    weights_dp = zeros(nsamples,nleaves,nreps,nchannels,nphases,'single');
end

t = ticker('binning',0,nphases);
for i = 1:nphases
    d = abs(i-(1:nphases));
    d(d>(nphases/2)) = abs(d(d>(nphases/2)) - nphases);
    e = reshape(1./exp(1).^d,1,1,nphases);
    
    for j = 1:length(binned{i})
        rep = 0;
        while intind(binned{i}(j)) - (rep+1)*nleaves > 0
            rep = rep + 1;
        end
        int = intind(binned{i}(j)) - (rep)*nleaves;
        nrep = mod(rep,nreps)+1;
        
        if dpind(binned{i}(j)) == 0
            weights(:,int,nrep,:,:) = weights(:,int,nrep,:,:) + repmat(reshape(e,1,1,1,1,nphases),[nsamples 1 1 nchannels 1]);       
            rawdata2(:,int,nrep,:,:) = rawdata2(:,int,nrep,:,:) + reshape(rawdata(1:nsamples,:,binned{i}(j)) .* e,nsamples,1,1,nchannels,nphases);
        elseif dpind(binned{i}(j)) >= 1 
            weights_dp(:,int,nrep,:,:) = weights_dp(:,int,nrep,:,:) + repmat(reshape(e,1,1,1,1,nphases),[nsamples 1 1 nchannels 1]);
            rawdata2_dp(:,int,nrep,:,:) = rawdata2_dp(:,int,nrep,:,:) + reshape(rawdata(1:nsamples,:,binned{i}(j)) .* e,nsamples,1,1,nchannels,nphases);
        end
    end
    t = ticker(t);
end
rawdata2 = reshape(rawdata2 ./ weights ,nsamples*nleaves*nreps,nchannels,nphases);
if nres > 1
    rawdata2_dp = reshape(rawdata2_dp ./ weights_dp ,nsamples*nleaves*nreps,nchannels,nphases);
end

clearvars weights weights_dp

%% Shift Raw data
shiftk = [0 0 0];
% shiftk = [0 0 80]; % Rabbit 4 March 
res = fov/matsize;

for i = 1:3
    rawdata2 = rawdata2 .* exp(complex(0,(shiftk(i)/res)*2*pi*KSpaceCoor(:,i).*res));
    if nres > 1
        rawdata2_dp = rawdata2_dp .* exp(complex(0,(shiftk(i)/res)*2*pi*KSpaceCoor_DP(:,i).*res));
    end
end

%% Recon

% flags
res_flag = 0; % 1 = high, 0 = low
save_flag = 1; % 1 = save, 0 = don't save
img_flag = 2; % 0 = gas, 1 = dp, 2 = both

% recon grid res
if res_flag == 1
    % hires grid
    imgsize = matsize;
    os = 3;
    k = 5;
    beta = [];
    zfill = 1;
    fname = sprintf('img_dyn_%iph',nphases);
elseif res_flag == 0
    % lowres grid
    imgsize = matsize; %[48 24 48]; % sag cor ax
    os = 3;
    k = 2;
    beta = []; 
    zfill = 1;
    fname = sprintf('img_dyn_%iph_lowres',nphases);
end

% image size
if length(imgsize) == 1
    imgsize = repmat(imgsize,[1 3]);
end

if img_flag == 0 || img_flag == 2
    [Ind,Dist,wi] = grid_lookup_20230113(KSpaceCoor,imgsize,fov,'os',os,'kernelsize',k,'beta',beta,'zfill',zfill);
end
if img_flag == 1 || img_flag == 2
    [Ind2,Dist2,wi2] = grid_lookup_20230113(KSpaceCoor_DP,imgsize,fov,'os',os,'kernelsize',k,'beta',beta,'zfill',zfill);
end

if img_flag == 0 || img_flag == 2
    img_gp_ch  = zeros([imgsize.*zfill,nchannels,nphases]);
end
if img_flag == 1 || img_flag == 2
    img_dp_ch = zeros([imgsize.*zfill,nchannels,nphases]);
end
toc;
tic;
parfor i = 1:nphases
% for i = 1:nphases
    fprintf('RECON PHASE %i/%i',i,nphases);

    if img_flag == 0 || img_flag == 2
        img_gp_ch(:,:,:,:,i)  = gridrecon_fa_20230113(KSpaceCoor   ,rawdata2(:,:,i)    ,imgsize,fov,'wi',wi ,'Ind',Ind ,'Dist',Dist ,'os',os,'k',k,'beta',beta,'zfill',zfill);
    end
    if img_flag == 1 || img_flag == 2
        img_dp_ch(:,:,:,:,i) = gridrecon_fa_20230113(KSpaceCoor_DP,rawdata2_dp(:,:,i),imgsize,fov,'wi',wi2,'Ind',Ind2,'Dist',Dist2,'os',os,'k',k,'beta',beta,'zfill',zfill);
    end

    fprintf(' (%.6g s)\n',toc);
end
toc;

% combine channels
fprintf('Combining channels...\n');
if img_flag == 0 || img_flag == 2
    tic;
    [img_gp,~,b] = combinecoils_fa(img_gp_ch);
    toc;
end
if img_flag == 1 || img_flag == 2
    img_dp = combinecoils_fa(img_dp_ch);
end
fprintf('Finished (%.6g s)\n',toc);

% rearrange
if img_flag == 0 || img_flag == 2
    img_gp = fliplr(permute(img_gp,[1 3 2 4]));
    img_gp_ch = fliplr(permute(img_gp_ch,[1 3 2 4 5]));
    b = fliplr(permute(b,[2 4 3 1]));
end
if img_flag == 1 || img_flag == 2
    img_dp = fliplr(permute(img_dp,[1 3 2 4]));
    img_dp_ch = fliplr(permute(img_dp_ch,[1 3 2 4 5]));
end

% shift phases
if img_flag == 0 || img_flag == 2
    [~,shift] = min(squeeze(sum(img_gp,[1 2 3])));
    shift = shift - 1;
    fprintf('Peak Shift = %i\n',shift);
    img_gp = circshift(img_gp,-shift,4);
    img_gp_ch = circshift(img_gp_ch,-shift,5);

    if img_flag == 2
        img_dp = circshift(img_dp,-shift,4);
        img_dp_ch = circshift(img_dp_ch,-shift,5);
    end
end

% save
if save_flag == 1
    fprintf('Saving...\n');
    path = [fileparts(files{1}),'/'];
    if img_flag == 0
%         save([path,fname,'_ch'],'img_gp_ch');
        save([path,fname],'img_gp');
    elseif img_flag == 1
%         save([path,fname,'_ch'],'img_dp_ch');
        save([path,fname],'img_dp');
    else
%         save([path,fname,'_ch'],'img_gp_ch','img_dp_ch');
        save([path,fname],'img_gp','img_dp');
    end
    fprintf('Saved to %s\n',path);
end

% display
% KR: Requires ImageViewer installation
% if img_flag == 1 || img_flag == 2
%     if img_flag == 1
%         imageViewer(img_dp);
%     else
%         imageViewer(img_gp,img_dp);
%     end
% elseif img_flag == 0 
%     imageViewer(img_gp);
% end

 %% Select Indices to Display in Movies
inds_cor = floor(linspace(59,34,10));
inds_ax = floor(linspace(31,58,10));
inds_sag = floor([linspace(13,31,5) linspace(50,65,5)]);

%% Make Movie (gas)

fig = figure('color','black','units','normalized','outerposition',[0 0 1 1]);
if res_flag == 0
    v = VideoWriter('dynGP_lowres','MPEG-4');
else
    v = VideoWriter('dynGP','MPEG-4');
end
v.FrameRate = 10;
v.Quality = 100;
open(v);
for i = 1:size(img_gp,4)
    tmp1 = rot90(img_gp(:,:,inds_cor,i),1);
    tmp2 = rot90(permute(img_gp(:,inds_ax,:,i),[1 3 2 4]),1);
    tmp3 = rot90(permute(img_gp(inds_sag,:,:,i),[2 3 1 4]),2);
    
    imagesc(wrapImage_fa(cat(3,tmp1,tmp2,tmp3),3,length(inds_cor),'hgap',0,'vgap',0,'scale',1),[1 100]);
    axis image; axis off; colormap('gray'); colorbar('box','off','color','w','FontSize',18);
    
    annotation('textbox',[.05 0.52 .2 .2],'String','Anterior','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    annotation('textbox',[.74 0.52 .2 .2],'String','Posterior','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    
    annotation('textbox',[.05 0.39 .2 .2],'String','Base','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    annotation('textbox',[.74 0.39 .2 .2],'String','Apex','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    
    annotation('textbox',[.05 0.26 .2 .2],'String','Right','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    annotation('textbox',[.74 0.26 .2 .2],'String','Left','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);

    drawnow;
    frame=getframe(fig);
    writeVideo(v,frame);
    delete(findall(gcf,'type','annotation'));
end
close(v);
close(fig);

%% Make Movie (dissolved)

fig = figure('color','black','units','normalized','outerposition',[0 0 1 1]);
if res_flag == 0
    v = VideoWriter('dynDP_lowres','MPEG-4');
else
    v = VideoWriter('dynDP','MPEG-4');
end
v.FrameRate = 10;
v.Quality = 100;
open(v);
for i = 1:size(img_gp,4)
    tmp1 = rot90(img_dp(:,:,inds_cor,i),1);
    tmp2 = rot90(permute(img_dp(:,inds_ax,:,i),[1 3 2 4]),1);
    tmp3 = rot90(permute(img_dp(inds_sag,:,:,i),[2 3 1 4]),2);
    
    imagesc(wrapImage_fa(cat(3,tmp1,tmp2,tmp3),3,length(inds_cor),'hgap',0,'vgap',0,'scale',1),[0.5 30]);
    axis image; axis off; colormap('gray'); colorbar('box','off','color','w','FontSize',18);
    
    annotation('textbox',[.05 0.52 .2 .2],'String','Anterior','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    annotation('textbox',[.74 0.52 .2 .2],'String','Posterior','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    
    annotation('textbox',[.05 0.39 .2 .2],'String','Base','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    annotation('textbox',[.74 0.39 .2 .2],'String','Apex','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    
    annotation('textbox',[.05 0.26 .2 .2],'String','Right','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    annotation('textbox',[.74 0.26 .2 .2],'String','Left','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);

    drawnow;
    frame=getframe(fig);
    writeVideo(v,frame);
    delete(findall(gcf,'type','annotation'));
end
close(v);
close(fig);

%% Make Movie (DP:GP)

% cmap = jet;
cmap = crameri('-roma');
cmap(1,:) = 0;

mask = zeros(size(img_gp));
z = zeros(size(img_gp,4),1);
for j = 1:size(img_gp,4)
     zz = reshape(img_gp(:,:,1:15,j),[],1);
     for i = 100:-0.1:1
        if sum(zz(zz>i))<=1e2
            z(j) = i;
        end
     end
     mask(:,:,:,j) = img_gp(:,:,:,j) > z(j);
end

dpgp = img_dp ./ img_gp .* mask .* (sind(data(1).FA)/sind(data(1).FA2));

fig = figure('color','black','units','normalized','outerposition',[0 0 1 1]);
if res_flag == 0
    v = VideoWriter('dynDPGP_lowres','MPEG-4');
else
    v = VideoWriter('dynDPGP','MPEG-4');
end
v.FrameRate = 10;
v.Quality = 100;
open(v);
for i = 1:size(dpgp,4)
    tmp1 = rot90(dpgp(:,:,inds_cor,i),1);
    tmp2 = rot90(permute(dpgp(:,inds_ax,:,i),[1 3 2 4]),1);
    tmp3 = rot90(permute(dpgp(inds_sag,:,:,i),[2 3 1 4]),2);
    
    imagesc(wrapImage_fa(cat(3,tmp1,tmp2,tmp3),3,length(inds_cor),'hgap',0,'vgap',10,'scale',1),[0 0.08]);
    axis image; axis off; colormap(cmap); colorbar('box','off','color','w','FontSize',18);
    
    annotation('textbox',[.05 0.54 .2 .2],'String','Anterior','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    annotation('textbox',[.74 0.54 .2 .2],'String','Posterior','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    
    annotation('textbox',[.05 0.39 .2 .2],'String','Base','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    annotation('textbox',[.74 0.39 .2 .2],'String','Apex','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    
    annotation('textbox',[.05 0.24 .2 .2],'String','Right','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    annotation('textbox',[.74 0.24 .2 .2],'String','Left','FitBoxToText','on',...
        'Color',[1 1 1],'LineStyle','none','HorizontalAlignment','center','FontSize',16);
    
    drawnow;
    frame=getframe(fig);
    writeVideo(v,frame);
    delete(findall(gcf,'type','annotation'));
end
close(v);
close(fig);
%% Get Spectra

iDataSet = 1;

delays = 1:1:9; 
TEs = (delays-0.5) .* data(iDataSet).dtSpec; % us. Define TE at the center of the sampling point (i.e. 50% of the dwell time)
TEs = TEs + (data(iDataSet).RFdur2/2)*1000;
TEs = ceil(TEs / GradRasterTime_us) * GradRasterTime_us;    % Round up to the nearest gradient raster time
MeasTE_us = data(iDataSet).TE * 1000;

spectra = complex(0,zeros([size(rawspec) length(TEs)],'single'));
tic;
t = ticker('iter',0,length(delays)*size(spectra,3));
for i = 1:length(delays)
    for j = 1:size(spectra,3)        
        for k = 1:size(spectra,2)
            spectra(:,k,j,i) = fftshift(fft(rawspec(1+delays(i):end,k,j),512));
        end

        t = ticker(t);
    end
end
toc;

spectra_avg = squeeze(mean(spectra(:,:,390:end,:),[2 3]));

%% Fit Spectra
% calculate frequency
fs = 1/(data(iDataSet).dtSpec*1e-6); % Hz
freq = (-ceil(data(iDataSet).nCol/2):ceil(data(iDataSet).nCol/2)-1) * fs / data(iDataSet).nCol / (data(iDataSet).frequency*(1+data(iDataSet).freq2*1e-6)*1e-6); % ppm
IdxNeg100ppm = find(freq >=  -100, 1, 'first');
[~, GpIdx] = max(sum(abs(spectra_avg(1:IdxNeg100ppm,:)),2));      % Find current GP frequency and shift it to 0 ppm
freq = freq - freq(GpIdx);

window = 170:340;
pklocs = [68 76 88];

startPoint = [20, freq(170+pklocs(1)-1), 3, 1, -pi, 0, 0, 0;
              20, freq(170+pklocs(2)-1), 3, 1, -pi, 0, 0, 0;
              10, freq(170+pklocs(3)-1), 5, 1, -pi, 0, 0, 0];
startPoint = [];


fits_rbc = zeros(size(spectra_avg,2),9);
fits_mem = zeros(size(spectra_avg,2),9);
fits_pk3 = zeros(size(spectra_avg,2),9);
spec_plots = zeros(length(window),size(spectra_avg,2));
spec_areas = zeros(length(pklocs),size(spectra_avg,2));

lsq_opts = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
    'MaxFunctionEvaluations',5e5,'MaxIterations',5e5,'Display','off',...
    'FunctionTolerance',1e-8,'StepTolerance',1e-9);

tic;
parfor i = 1:size(spectra_avg,2)
    fit = fitVoigt(freq(window),spectra_avg(window,i),peaks_locs=pklocs,opts=lsq_opts,...
        startPoint=startPoint,...
        ub=[ Inf, Inf, 20, 1,  pi, 0,  Inf,  Inf], ...
        lb=[   0,-Inf,  0, 1, -pi, 0, -Inf, -Inf]);
    fits_pk3(i,:) = fit(1,:);
    fits_mem(i,:) = fit(2,:);
    fits_rbc(i,:) = fit(3,:);
    [spec_plots(:,i), spec_areas(:,i)] = evalVoigt(fit,freq(window));
end
NumDpRandomize    = 5;
NumDpFitParameters= 18;                           % Number of stored fitting parameters for DP fitting
rng('shuffle');                                   % Seed random number generator based on time
Idx160ppm         = find(freq >=  160, 1, 'first');
Idx200ppm         = find(freq >=  200, 1, 'first');
Idx260ppm         = find(freq >=  260, 1, 'first');
FreqVec           = [freq(Idx160ppm:Idx260ppm)];
FitDP_DataVec     = zeros(length(FreqVec), size(spectra_avg,2));    % Spectroscopic data to be fitted
FitDP_FittedCurve = zeros(length(FreqVec), size(spectra_avg,2));    % Best fit of the spectroscopic data
TE_Best_us        = zeros(size(spectra_avg,2), 1);
Theta_deg         = zeros(size(spectra_avg,2), 1);
RBC_Mem_Ratio     = zeros(size(spectra_avg,2), 1);
% LowerBounds = [   0, 0.1, 212, -2*pi,   0, 0.1, 197.1, -2*pi,   0, 0.1, 191, -2*pi,   0, 0.1, 186.5, -2*pi,-inf, -inf];
% UpperBounds = [ inf,  10, 225,  2*pi, inf,  10, 206.0,  2*pi, inf,  10, 197,  2*pi, inf, 10, 190.5,  2*pi, inf,  inf];
LowerBounds = [   0, 0.1, 213, -2*pi,   0, 0.1, 212.5, -2*pi,   0, 0.1, 197.1, -2*pi,   0, 0.1, 191, -2*pi,-inf, -inf];
UpperBounds = [ inf,  10, 225,  2*pi, inf,  10, 205.0,  2*pi, inf,  10, 206.0,  2*pi, inf,  10, 197,  2*pi, inf,  inf];
LowerBounds = [   0, 0.1, 213, -2*pi,   0, 0.1, 212.5, -2*pi,   0, 0.1, 195.1, -2*pi,   0, 0.1, 191, -2*pi,-inf, -inf];
UpperBounds = [ inf,  10, 225,  2*pi, 0,  10, 205.0,  2*pi, inf,  10, 206.0,  2*pi, 0,  10, 197,  2*pi, inf,  inf];
warning('off', 'all')
parfor i = 1:size(spectra_avg,2)
    MinError = inf;
    DataVec = [spectra_avg(Idx160ppm:Idx260ppm, i)]';
    FitDP_DataVec(:,i) = DataVec;
    fprintf('Fitting spectrum %d of %d.\n', i, size(spectra_avg,2));
    for RandLoop = 1:NumDpRandomize
        % Create randomized starting parameters within plausible limits
        StartPointMultiLor = zeros (1,NumDpFitParameters);
        StartPointMultiLor(1) = max(abs(spectra_avg(Idx200ppm:Idx260ppm, i)));      % DP Peak 1 amplitude
        StartPointMultiLor(2) = 3.0;                                            % DP Peak 1 half-width
        StartPointMultiLor(2) = 5.0 + 8.0 * (rand()-0.5);                       % DP Peak 1 half-width
        StartPointMultiLor(3) = 218;                                            % DP Peak 1 chemical shift
        StartPointMultiLor(4) = 0.0;                                            % DP Peak 1 phase
        StartPointMultiLor(5) = max(abs(spectra_avg(Idx160ppm:Idx200ppm, i)));      % DP Peak 2 amplitude
        StartPointMultiLor(5) = 0;      % DP Peak 2 amplitude
        StartPointMultiLor(6) = 3.0;                                            % DP Peak 2 half-width
        StartPointMultiLor(7) = 203;                                            % DP Peak 2 chemical shift
        StartPointMultiLor(8) = 0.0;                                            % DP Peak 2 phase
        StartPointMultiLor(9) = max(abs(spectra_avg(Idx160ppm:Idx200ppm, i)));      % DP Peak 3 amplitude
        StartPointMultiLor(10)= 3.0;                                            % DP Peak 3 half-width
        StartPointMultiLor(11)= 196;                                            % DP Peak 3 chemical shift
        StartPointMultiLor(12)= 0.0;                                            % DP Peak 3 phase
        StartPointMultiLor(13) = max(abs(spectra_avg(Idx160ppm:Idx200ppm, i)));      % DP Peak 4 amplitude
        StartPointMultiLor(14)= 1.0;                                            % DP Peak 4 half-width
        StartPointMultiLor(15)= 189;                                            % DP Peak 4 chemical shift
        StartPointMultiLor(16)= 0.0;                                            % DP Peak 4 phase
        StartPointMultiLor(17)= 0.0;                                            % Offset (Real)
        StartPointMultiLor(18)= 0.0;                                            % Offset (Imag)

        % Fit spectral data to multiple complex Lorentzian functions
        DpPeakEstimates = Fit4LorentzianCplxPh_Dixon_Con_20230126(FreqVec, DataVec, StartPointMultiLor, LowerBounds, UpperBounds);

        Amp1        = DpPeakEstimates(1);
        HalfWidthL1 = DpPeakEstimates(2);
        Shift1      = DpPeakEstimates(3);
        Phase1      = DpPeakEstimates(4);
        Amp2        = DpPeakEstimates(5);
        HalfWidthL2 = DpPeakEstimates(6);
        Shift2      = DpPeakEstimates(7);
        Phase2      = DpPeakEstimates(8);
        Amp3        = DpPeakEstimates(9);
        HalfWidthL3 = DpPeakEstimates(10);
        Shift3      = DpPeakEstimates(11);
        Phase3      = DpPeakEstimates(12);
        Amp4        = DpPeakEstimates(13);
        HalfWidthL4 = DpPeakEstimates(14);
        Shift4      = DpPeakEstimates(15);
        Phase4      = DpPeakEstimates(16);
        Offset_R    = DpPeakEstimates(17);
        Offset_I    = DpPeakEstimates(18);
        Offset      = complex(Offset_R, Offset_I);

        FitLorentzian1 = LorentzianFun(FreqVec, DpPeakEstimates(1:4));
        FitLorentzian2 = LorentzianFun(FreqVec, DpPeakEstimates(5:8));
        FitLorentzian3 = LorentzianFun(FreqVec, DpPeakEstimates(9:12));
        FitLorentzian4 = LorentzianFun(FreqVec, DpPeakEstimates(13:16));
        FittedCurve = squeeze(FitLorentzian1(3,:) + FitLorentzian2(3,:) + FitLorentzian3(3,:) + FitLorentzian4(3,:)) + Offset;
        ErrorVector = FittedCurve - DataVec;
        sse = sum(abs(ErrorVector) .^ 2);

        if (MinError/sse > 1.001)
            MinError                = sse;
            FitDP_FittedCurve(:, i) = FittedCurve;
            PhaseDiff_Best          = Phase3 - Phase1;
            FreqDiff_Best_ppm       = Shift1 - Shift3;
            FreqDiff_Best_Hz        = FreqDiff_Best_ppm * data(1).frequency / 1e6;

            % Calculate TE to achieve a phase difference of pi/2
            TE_Best_us(i)           = TEs(i) + 1e6 * (pi/2 - PhaseDiff_Best) / (2 * pi * FreqDiff_Best_Hz);

            % Calculate phase angle Theta at the imaging TE
            Theta_deg(i)            = 360 * (MeasTE_us - TEs(i)) * FreqDiff_Best_Hz / 1e6;
            Theta_deg(i)            = Theta_deg(i) + PhaseDiff_Best * 360 / (2 * pi);
            RBC_Signal              = sum(abs(FitLorentzian1(3,:)));
            Mem_Signal              = sum(abs(FitLorentzian3(3,:)));
            RBC_Mem_Ratio(i)        = RBC_Signal / Mem_Signal;
        end
    end
end
warning('on', 'all')
toc;

figure;
for i = 1:length(delays)
    subplot(3,3,i); 
    hold on;
    plot(freq(Idx160ppm:Idx260ppm),abs(spectra_avg(Idx160ppm:Idx260ppm,i)),'k','LineWidth',2);
    plot(freq(Idx160ppm:Idx260ppm),abs(FitDP_FittedCurve(:,i)),'r','LineWidth',2);
    title(['TE = ',num2str(TEs(i))]);
%     xlabel('Frequency (ppm)','FontSize',18,'FontWeight','bold');
%     ylabel('Signal Intensity','FontSize',18,'FontWeight','bold');
    set(gca,'XLim',[160 260],'YLim',[0 inf]);
    hold off;
end


%% Plot Phase Differences
% ginds = 3:6;
% 
% pd = nan(length(delays),1);
% for i = 1:length(delays)
%     pd(i) = (fits_rbc(i,5)-fits_mem(i,5))*180/pi;
% end
% TE_est = linspace(TEs(1),TEs(end),100);
% pd_est = polyval(polyfit(TEs(ginds),pd(ginds),1),TE_est);
% [~, m] = min(abs(abs(pd_est)-90));
% theta = polyval(polyfit(TEs(ginds),pd(ginds),1),MeasTE_us);
% R = polyval(polyfit(TEs(ginds),spec_areas(3,ginds)./spec_areas(2,ginds),1),MeasTE_us);
% disp(['Optimal TE = ',num2str(TE_est(m),4),' us']);
% disp(['Angle at 620 us = ',num2str(theta,4)]);
% disp(['RBC:Mem at 620 us = ',num2str(R)]);
% figure;
% hold on;
% scatter(TEs,abs(pd),64,'filled');
% plot(TE_est,abs(pd_est),'--r','LineWidth',2);
% title('LSQNONLIN');
% ylabel('Phase Difference (°)');
% xlabel('TE (us)');
% set(gca,'FontWeight','bold','FontSize',16);
% hold off;

% Update theta and R based on Lorentzian fitting
theta = Theta_deg(1);
R = RBC_Mem_Ratio(1);

%% Phasemap and Mask

mask = zeros(size(img_gp));
z = zeros(size(img_gp,4),1);
for j = 1:size(img_gp,4)
    zz = reshape(img_gp(:,:,1:15,j),[],1);
    for i = 40:-0.1:5
        if sum(zz(zz>i))<=1e1
            z(j) = i;
        end
    end
    mask(:,:,:,j) = img_gp(:,:,:,j) > z(j);
end
mask(mask==0) = nan;

% pmap = angle(squeeze(img_gp_ch));
% pmap = pmap - mean(pmap.*mask,[1 2 3],'omitnan');
% 
% % img_dp_ch2 = squeeze(img_dp_ch) .* exp(-1i.*pmap) .* mask;
% img_dp_ch2 = squeeze(img_dp_ch) .* conj(b) .* mask;

maskdp = single(img_dp >= 6);
maskdp(maskdp==0) = nan;

img_dp_ch2 = squeeze(img_dp_ch) .* mask .* maskdp;

%% 1 Point Dixon

max_ratio = zeros(1,2);
max_ratio(:,1) = 1e6;

irange = -3.14:0.01:3.14;
vals = zeros(length(irange),6);
rdiff = zeros(length(irange),1);

t = ticker('iter',0,length(irange));
tic;
parfor i = 1:length(irange)
% for i = 1:length(irange)
    vals_temp = zeros(1,6);

    vals_temp(3) = 1;
    vals_temp(4) = 1;

    tmp = img_dp_ch2.*exp(1i*irange(i));

    tmprbc = real(tmp);
    tmpmem = imag(tmp);

    vals_temp(3) = mean(tmprbc(:),'omitnan');
    vals_temp(4) = mean(tmpmem(:),'omitnan');

    tmpmem = tmpmem ./ sind(theta);
    tmprbc = tmprbc - (tmpmem .* cosd(theta));

    tmpmem = abs(tmpmem);
    tmprbc = abs(tmprbc);

    vals_temp(5) = mean(tmprbc(:),'omitnan');
    vals_temp(6) = mean(tmpmem(:),'omitnan');

    tmpratio = tmprbc./tmpmem;
    tmpratio(abs(tmpratio) > 2) = nan;

    vals_temp(1) = mean(tmpratio(:),'omitnan');
    vals_temp(2) = sum(tmprbc(:),'omitnan') ./ sum(tmpmem(:),'omitnan');
    vals(i, :) = vals_temp;
    rdiff(i) = sqrt((R - mean(tmpratio(:),'omitnan')).^2);
end
for i = 1:length(irange)
    if rdiff(i) <= max_ratio(1)
        max_ratio(1) = rdiff(i);
        max_ratio(2) = irange(i);
    end
    t = ticker(t);
end
toc;

figure; 
subplot(1,2,1); hold on;
plot(irange,vals(:,1),'b');
plot(irange,vals(:,2),'r');
hold off; legend({'Mean','Sums'});
subplot(1,2,2); hold on;
plot(irange,vals(:,3),'r');
plot(irange,vals(:,4),'b');
plot(irange,vals(:,5),'m');
plot(irange,vals(:,6),'c');
hold off; legend({'Real','Imag','RBC','Mem'});

%% Phase Correct and Display
% img_dp_ph = zeros(80,80,80,16,20);
% phs = linspace(0,pi,20);
% for i = 1:20
%     img_dp_ph(:,:,:,:,i) = img_dp_ch2 .* exp(1i.*phs(i));
% end

img_dp_ph = squeeze(img_dp_ch) .* exp(1i.*max_ratio(2)) .* mask;

img_rbc = real(img_dp_ph);
img_mem = imag(img_dp_ph);

img_mem = img_mem ./ sind(theta);
img_rbc = img_rbc - img_mem .* cosd(theta);

img_rbc = abs(img_rbc);
img_mem = abs(img_mem);

img_rbcmem = img_rbc./img_mem;
img_rbcmem(img_rbcmem<0 | img_rbcmem>2) = nan;

img_rbcgp = img_rbc./img_gp.*100.*sind(data(1).FA2)/sind(data(1).FA);
img_memgp = img_mem./img_gp.*100.*sind(data(1).FA2)/sind(data(1).FA);

% imageViewer(img_rbcgp,img_memgp);

figure('Color','white');
plot(squeeze(mean(img_rbcgp.*maskdp,[1 2 3],'omitnan')),'r','LineWidth',1.5); hold on;
plot(squeeze(mean(img_gp,[1 2 3],'omitnan'))./max(squeeze(mean(img_gp,[1 2 3],'omitnan'))).*1.3,':','Color',[0.6 0.6 0.6],'LineWidth',1.5');
plot(squeeze(mean(img_memgp.*maskdp,[1 2 3],'omitnan')),'g','LineWidth',1.5); hold off;
set(gca,'XLim',[1 16]);
ylabel('RBC-/Mem-Gas Ratio (%)');
xlabel('Respiratory Phase');
yyaxis right;
plot(squeeze(mean(img_rbcmem.*maskdp,[1 2 3],'omitnan')),'Color',[0.93,0.69,0.13],'LineWidth',1.5);
set(gca,'XLim',[1 16],'FontWeight','bold','FontSize',16,'YColor',[0.93,0.69,0.13]);
ylabel('RBC:Mem Ratio');
legend({'RBC:GP','','Mem:GP','RBC:Mem'},'Location','north','Orientation','horizontal','FontSize',10,'box','off');

% save([fileparts(files{1}),'\','img_dyn_16ph_ratios'],'img_rbcgp','img_memgp','img_rbcmem');

%% Make Movie (Ratios)

cmap_mem = crameri('bamako');
cmap_mem(1,:) = 0;

cmap_rbc = crameri('-lajolla');
cmap_rbc(1,:) = 0;

cmap_rbcmem = crameri('-roma');
cmap_rbcmem(1,:) = 0;

fig = figure('color','black','units','normalized','outerposition',[0 0 1 1]);
if res_flag == 0
    v = VideoWriter('dynRatios_lowres','MPEG-4');
else
    v = VideoWriter('dynRatios','MPEG-4');
end
v.FrameRate = 10;
v.Quality = 100;
open(v);
for i = 1:size(img_gp,4)
    tmp1 = rot90(img_gp(:,:,inds_cor,i),1);
    tmp2 = rot90(img_rbcgp(:,:,inds_cor,i),1);
    tmp3 = rot90(img_memgp(:,:,inds_cor,i),1);
    tmp4 = rot90(img_rbcmem(:,:,inds_cor,i),1);

    ax1 = subplot(4,1,1);
    imagesc(wrapImage_fa(tmp1,1,length(inds_cor),'hgap',0,'vgap',10,'scale',1),[1 100]);
    axis image; axis off; colormap(ax1,gray); colorbar('Box','off','color','w','FontSize',18);
    ylabel('GP','Color','w','FontSize',18,'Rotation',0,'Visible','on');
        
    ax2 = subplot(4,1,2);
    imagesc(wrapImage_fa(tmp2,1,length(inds_cor),'hgap',0,'vgap',10,'scale',1),[0 1]);
    axis image; axis off; colormap(ax2,cmap_rbc); colorbar('Box','off','color','w','FontSize',18);
    ylabel('RBC:GP','Color','w','FontSize',18,'Rotation',0,'Visible','on');

    ax3 = subplot(4,1,3);
    imagesc(wrapImage_fa(tmp3,1,length(inds_cor),'hgap',0,'vgap',10,'scale',1),[0 2]);
    axis image; axis off; colormap(ax3,cmap_mem); colorbar('Box','off','color','w','FontSize',18);
    ylabel('Mem:GP','Color','w','FontSize',18,'Rotation',0,'Visible','on');

    ax4 = subplot(4,1,4);
    imagesc(wrapImage_fa(tmp4,1,length(inds_cor),'hgap',0,'vgap',10,'scale',1),[0 0.6]);
    axis image; axis off; colormap(ax4,cmap_rbcmem); colorbar('Box','off','color','w','FontSize',18);
    ylabel('RBC:Mem','Color','w','FontSize',18,'Rotation',0,'Visible','on');
    
    drawnow;
    frame=getframe(fig);
    writeVideo(v,frame);
end
close(v);
close(fig);

%% Plot 1 slice
figure('Color','black','Position',[64,350,1780,420]);
subplot(3,1,1);
imagesc(wrapImage_fa(squeeze(rot90(img_rbcgp(:,:,47,:),1)),1,16,'hgap',0,'vgap',0),[0 1]);
axis image; axis off; colormap(cmap_rbcmem); colorbar('Box','off','color','w','FontSize',18);
ylabel('RBC:GP','Color','White','FontWeight','Bold','FontSize',16,'Visible','on');
subplot(3,1,2);
imagesc(wrapImage_fa(squeeze(rot90(img_memgp(:,:,47,:),1)),1,16,'hgap',0,'vgap',0),[0 1.5]);
axis image; axis off; colormap(cmap_rbcmem); colorbar('Box','off','color','w','FontSize',18);
ylabel('Mem:GP','Color','White','FontWeight','Bold','FontSize',16,'Visible','on');
subplot(3,1,3);
imagesc(wrapImage_fa(squeeze(rot90(img_rbcmem(:,:,47,:),1)),1,16,'hgap',0,'vgap',0),[0 0.8]);
axis image; axis off; colormap(cmap_rbcmem); colorbar('Box','off','color','w','FontSize',18);
ylabel('RBC:Mem','Color','White','FontWeight','Bold','FontSize',16,'Visible','on');

function Lorentzian = LorentzianFun(Freq, params)
    Lorentzian  = zeros(3, length(Freq));
    Amp         = params(1);
    HalfWidth   = params(2);
    Shift       = params(3);
    Phase       = params(4);

    Lorentzian(1,:) = HalfWidth ./ ((Freq - Shift).^2 + HalfWidth^2);                 % Real part
    % Lorentzian(2,:) = -(Freq - Shift) ./ ((Freq - Shift).^2 + HalfWidth^2);           % Imaginary part
    Lorentzian(2,:) = (Freq - Shift) ./ ((Freq - Shift).^2 + HalfWidth^2);           % Imaginary part
    Lorentzian(3,:) = Amp * complex(Lorentzian(1,:), Lorentzian(2,:)) * exp(complex(0,Phase));    % Complex Lorentzian
end

