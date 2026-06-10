% spiral_human_20230824.m
% Subject: 030DN (LTX)
% Methods:
%   Free-breathing
%   UVA coil
%   
% Notes: 9 mo FU; no spectra and ventilator not synced because trigger wasn't working
%
%
%

clear;
spokes = 26;

%% Load Images

files = [
    % Free-breathing
    
"/Users/hoomham/Hooman/Images/2026-05-13_SyrT1/meas_MID00159_FID16681_fa_spiral_dyn_fancy_v3_20240130.dat"
    % EI breath-hold
% "/Users/hoomham/Hooman/Work/Analysis/2023-11-03_000LL/rec/meas_MID02793_FID05038_fa_spiral_dyn_fancy_v3_20230821.dat"
];

[data,twixs] = recon_20210622(files,'nStudies',1:1,'recon',0);


%% Add to log
logfile = "C:\Users\faraz\Documents\General Lab Stuff\Random Stuff\patient_log.xlsx";
log = cell(length(data),10);
for i = 1:length(data)
    log{i,1} = 20230824;
    log{i,2} = '030DN';
    log{i,3} = 'LTX';
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
    log{i,9} = '9 mo FU';
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
    nspec = data(i).numSpec * spokes/20;
end
if strcmpi(data(i).acqOrder,'GP-DP-DP') && data(i).numRes > 1
    nres = 3;
end

KSpaceCoor = loadtrajectory3D(data(i).frequency/1.498/1e6,fov,nsamples,nleaves,nreps,matsize,'axial',[],[],'calibrations_3D_20220308.mat');
if nres > 1
    KSpaceCoor_DP = loadtrajectory3D(data(i).frequency/1.498/1e6,fov,nsamples,nleaves,nreps,matsize,'axial','dissolved',[],'calibrations_3D_20220308.mat');
    KSpaceCoor_DP = KSpaceCoor;
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

%% Find peaks
range = 1:length(absval);

% findpeaks(absval(range(1):1:range(end)),'MinPeakDistance',20*nleaves,'MinPeakHeight',2e-3);
[~,loc] = findpeaks(absval(range(1):1:range(end)),'MinPeakDistance',8*nleaves,'MinPeakHeight',2e-4);
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
res = fov/matsize;

for i = 1:3
    rawdata2 = rawdata2 .* exp(complex(0,(shiftk(i)/res)*2*pi*KSpaceCoor(:,i).*res));
    if nres > 1
        rawdata2_dp = rawdata2_dp .* exp(complex(0,(shiftk(i)/res)*2*pi*KSpaceCoor_DP(:,i).*res));
    end
end

%% Recon

% flags
res_flag = 1; % 1 = high, 0 = low
save_flag = 1; % 1 = save, 0 = don't save
img_flag = 0; % 0 = gas, 1 = dp, 2 = both

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

tic;
for i = 1:nphases
    fprintf('RECON PHASE %i/%i',i,nphases);

    if img_flag == 0 || img_flag == 2
        img_gp_ch(:,:,:,:,i)  = gridrecon_fa_20230113(KSpaceCoor   ,rawdata2(:,:,i)    ,imgsize,fov,'wi',wi ,'Ind',Ind ,'Dist',Dist ,'os',os,'k',k,'beta',beta,'zfill',zfill);
    end
    if img_flag == 1 || img_flag == 2
        img_dp_ch(:,:,:,:,i) = gridrecon_fa_20230113(KSpaceCoor_DP,rawdata2_dp(:,:,i),imgsize,fov,'wi',wi2,'Ind',Ind2,'Dist',Dist2,'os',os,'k',k,'beta',beta,'zfill',zfill);
    end

    fprintf(' (%.6g s)\n',toc);
end

% combine channels
fprintf('Combining channels...\n');
if img_flag == 0 || img_flag == 2
    [img_gp,~,b] = combinecoils_fa(img_gp_ch);
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
    path = [fileparts(files{1}),'\'];
    if img_flag == 0
        save([path,fname,'_ch'],'img_gp_ch');
        save([path,fname],'img_gp');
    elseif img_flag == 1
        save([path,fname,'_ch'],'img_dp_ch');
        save([path,fname],'img_dp');
    else
        save([path,fname,'_ch'],'img_gp_ch','img_dp_ch');
        save([path,fname],'img_gp','img_dp');
    end
    fprintf('Saved to %s\n',path);
end

% display
% if img_flag == 1 || img_flag == 2
%     if img_flag == 1
%         imageViewer(img_dp);
%     else
%         imageViewer(img_gp,img_dp);
%     end
% elseif img_flag == 0 
%     imageViewer(img_gp);
% end

%% Make Movie (gas)

inds_cor = floor(linspace(65,39,10));
inds_ax = floor(linspace(24,54,10));
inds_sag = floor([linspace(18,33,5) linspace(51,65,5)]);

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
     for i = 100:-0.1:10
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
    
    imagesc(wrapImage_fa(cat(3,tmp1,tmp2,tmp3),3,length(inds_cor),'hgap',0,'vgap',10,'scale',1),[0 0.05]);
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

% calculate frequency
i = 1;
fs = 1/(data(i).dtSpec*1e-6); % Hz
freq = (-ceil(data(i).nCol/2):ceil(data(i).nCol/2)-1) * fs / data(i).nCol / (data(i).frequency*(1+data(i).freq2*1e-6)*1e-6); % ppm
freq = freq - freq(139); % shift 0 to GP frequency

delays = 1:1:9; 
TEs = delays .* data(i).dtSpec; % us
TEs = TEs + (data(i).RFdur2/2)*1000;

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
for i = 1:size(spectra_avg,2)
    fit = fitVoigt(freq(window),spectra_avg(window,i),peaks_locs=pklocs,opts=lsq_opts,...
        startPoint=startPoint,...
        ub=[ Inf, Inf, 20, 1,  pi, 0,  Inf,  Inf], ...
        lb=[   0,-Inf,  0, 1, -pi, 0, -Inf, -Inf]);
    fits_pk3(i,:) = fit(1,:);
    fits_mem(i,:) = fit(2,:);
    fits_rbc(i,:) = fit(3,:);
    [spec_plots(:,i), spec_areas(:,i)] = evalVoigt(fit,freq(window));
end
toc;

figure;
for i = 1:length(delays)
    subplot(3,3,i); 
    hold on;
    plot(freq(window),abs(spectra_avg(window,i)),'k','LineWidth',2);
    plot(freq(window),abs(spec_plots(:,i)),'r','LineWidth',2);
    title(['TE = ',num2str(TEs(i))]);
%     xlabel('Frequency (ppm)','FontSize',18,'FontWeight','bold');
%     ylabel('Signal Intensity','FontSize',18,'FontWeight','bold');
    set(gca,'XLim',[160 240],'YLim',[0 inf]);
    hold off;
end


%% Plot Phase Differences
ginds = 1:7;

pd = nan(length(delays),1);
for i = 1:length(delays)
    pd(i) = (fits_rbc(i,5)-fits_mem(i,5))*180/pi;
end
TE_est = linspace(TEs(1),TEs(end),100);
pd_est = polyval(polyfit(TEs(ginds),pd(ginds),1),TE_est);
[~, m] = min(abs(abs(pd_est)-90));
theta = polyval(polyfit(TEs(ginds),pd(ginds),1),620);
R = polyval(polyfit(TEs(ginds),spec_areas(3,ginds)./spec_areas(2,ginds),1),620);
disp(['Optimal TE = ',num2str(TE_est(m),4),' us']);
disp(['Angle at 620 us = ',num2str(theta,4)]);
disp(['RBC:Mem at 620 us = ',num2str(R)]);
figure;
hold on;
scatter(TEs,abs(pd),64,'filled');
plot(TE_est,abs(pd_est),'--r','LineWidth',2);
title('LSQNONLIN');
ylabel('Phase Difference (°)');
xlabel('TE (us)');
set(gca,'FontWeight','bold','FontSize',16);
hold off;

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

t = ticker('iter',0,length(irange));
for i = 1:length(irange)
    tmp = img_dp_ch2.*exp(1i*irange(i));
        
    tmprbc = real(tmp);
    tmpmem = imag(tmp);

    vals(i,3) = mean(tmprbc(:),'omitnan');
    vals(i,4) = mean(tmpmem(:),'omitnan');

    tmpmem = tmpmem ./ sind(theta);
    tmprbc = tmprbc - (tmpmem .* cosd(theta));

    tmpmem = abs(tmpmem);
    tmprbc = abs(tmprbc);

    vals(i,5) = mean(tmprbc(:),'omitnan');
    vals(i,6) = mean(tmpmem(:),'omitnan');

    tmpratio = tmprbc./tmpmem;
    tmpratio(abs(tmpratio) > 2) = nan;
    
    vals(i,1) = mean(tmpratio(:),'omitnan');
    vals(i,2) = sum(tmprbc(:),'omitnan') ./ sum(tmpmem(:),'omitnan');

    rdiff = sqrt((R - mean(tmpratio(:),'omitnan')).^2);
    if rdiff <= max_ratio(1)
        max_ratio(1) = rdiff;
        max_ratio(2) = irange(i);
    end
    t = ticker(t);
end

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

%% Recon the Single-breath
win = nleaves * nreps * nres;
winnum = floor(size(rawdata_sb,3)/win);
rawdata_sb2 = zeros(nsamples,nleaves,nreps,nchannels,winnum);
weights_sb = zeros(nsamples,nleaves,nreps,nchannels,winnum);

for i = 1:win:size(rawdata_sb,3)-win
    for j = 0:win-1
        rep = 0;
        while intind_sb(i+j) - (rep+1)*nleaves > 0
            rep = rep + 1;
        end
        int = intind_sb(i+j) - (rep)*nleaves;
        nrep = mod(rep,nreps)+1;

        if dpind_sb(i+j) == 0
            weights_sb(:,int,nrep,:,ceil(i/win)) = weights_sb(:,int,nrep,:,ceil(i/win)) + 1;
            rawdata_sb2(:,int,nrep,:,ceil(i/win)) = rawdata_sb2(:,int,nrep,:,ceil(i/win)) + reshape(rawdata_sb(1:nsamples,:,i+j),nsamples,1,1,nchannels,1);
        end
    end
end
rawdata_sb2 = reshape(rawdata_sb2 ./ weights_sb ,nsamples*nleaves*nreps,nchannels,winnum);

clearvars weights_sb 

img_sb_ch  = zeros([imgsize,nchannels,winnum]);
tic;
for i = 1:winnum
    fprintf('RECON PHASE %i/%i',i,winnum);
    img_sb_ch(:,:,:,:,i)  = gridrecon_fa_20230113(KSpaceCoor,rawdata_sb2(:,:,i),imgsize,fov,'wi',wi,'Ind',Ind,'Dist',Dist,'os',os,'k',k,'beta',beta,'zfill',zfill);
    fprintf(' (%.6g s)\n',toc);
end

% combine channels
fprintf('Combining channels...\n');
[img_sb,~,b_sb] = combinecoils_fa(img_sb_ch);
fprintf('Finished (%.6g s)\n',toc);

% rearrange
img_sb = fliplr(permute(img_sb,[1 3 2 4]));
img_sb_ch = fliplr(permute(img_sb_ch,[1 3 2 4 5]));
b_sb = fliplr(permute(b_sb,[2 4 3 1]));

% save
if save_flag == 1
    fprintf('Saving...\n');
    path = [fileparts(files{1}),'\'];
    if res_flag == 1
        save([path,'img_sb_ch'],'img_sb_ch');
        save([path,'img_sb'],'img_sb');
    elseif res_flag == 0
        save([path,'img_sb_lowres_ch'],'img_sb_ch');
        save([path,'img_sb_lowres'],'img_sb');
    end
    fprintf('Saved to %s\n',path);
end

% display
imageViewer(img_sb);
