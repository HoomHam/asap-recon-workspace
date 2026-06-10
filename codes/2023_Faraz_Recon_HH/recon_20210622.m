function [allimages,twixs] = recon_20210622(aFileNames,Options)
% aFileNames            String or String array of input filenames (include full path)
% Options               name-value pairs with flags and parameters

%% Initialize

% parse arguments
arguments
    aFileNames string
    Options.nStudies double = 1:length(aFileNames)                          % selects which studies to reconstruct/analyze
    Options.shifts double = 0                                               % circshift for wrap-around (array to specify for each file, otherwise defaults to first value)
    Options.projection double = 0                                           % 1 to make images into projections
    Options.threshold = 'none'                                              % percent of max signal or 'hist' for thresholding
    Options.normalize double = 0                                            % 1 to convert images to grayscale
    Options.dcf double = 1                                                  % 0 for gridding DCF, 1 for Meyer DCF, 2 or 3 for Voronoi, 4 for no DCF
    Options.loadtraj double = 1                                             % 1 to check for trajectory calibration
    Options.spcalibfile string = 'calibrations_20201109.mat'                % mat file containing calibrated spiral trajectories
    Options.sp3Dcalibfile string = 'calibrations_3D_20220308.mat'           % mat file containing calibrated fancy 3D spiral trajectories
    Options.gridfile string = 'grid_lookup_20230113.mat'                    % mat file containing saved grids
    Options.zfill double = 1                                                % spiral gridding zero filling factor
    Options.GridOSFactor double = 3                                         % spiral gridding oversample factor
    Options.KernelSize double = 5                                           % spiral gridding kernel size
    Options.cmap = gray                                                     % colormap 
    Options.cscale double = []                                              % colorscale for images
    Options.combine double = 1                                              % 1 to combine slices/partitions into single image
    Options.scale double = 1                                                % interpolation factor for displaying images
    Options.dispimages double = 0                                           % 1 to display raw images
    Options.progressbar double = 1                                          % 1 to display progressbar
    Options.savegif double = 0                                              % 1 to save images as gif
    Options.saveraw double = 0                                              % 1 to save raw images   
    Options.savekspace double = 0                                           % 1 to save raw kspace
    Options.recon double = 1                                                % 1 to reconstruct images
    Options.verbose double = 1                                              % 1 to printout parameters
    Options.parseSpecial double = 1                                         % 1 to get parameters from the "special" tab
end     

% flags
projection = Options.projection;
normalize = Options.normalize;
dcf = Options.dcf;
loadtraj = Options.loadtraj;
combine = Options.combine;
dispimages = Options.dispimages; 
progressbar = Options.progressbar;
savegif = Options.savegif;
saveraw = Options.saveraw;
savekspace = Options.savekspace;
recon = Options.recon;
verbose = Options.verbose;
parseSpecial = Options.parseSpecial;

% recon parameters
spcalibfile = Options.spcalibfile;
sp3Dcalibfile = Options.sp3Dcalibfile;
gridfile = Options.gridfile;
zfill = Options.zfill;
GridOSFactor = Options.GridOSFactor;
KernelSize = Options.KernelSize;

% display parameters
threshold = Options.threshold; 
shifts = Options.shifts; 
scale = Options.scale; 

% colormaps
cscale = Options.cscale; 
cmap = Options.cmap;
cmap(1,:) = 0; % background color

% studies to analyze
nStudies = Options.nStudies; 

% initialize output
allimages = struct; 
twixs = struct;

%% Recon
tic;
enum = 1;
N = length(nStudies);
for jj = nStudies
    % initialize progress bar
    if progressbar
        if ~exist('h', 'var')
            h = waitbar(0,sprintf('Loading Data...'),'name',sprintf('Study #%i of %i',jj,length(aFileNames)));
        else
            waitbar((enum-1)/N,h,sprintf('Loading Data...'));
            h.Name = sprintf('Study #%i of %i',jj,length(aFileNames));
        end
    end

    % extract header data
    twix = mapVBVD(char(aFileNames(jj)));
    if progressbar, waitbar((0.05+enum-1)/N,h,sprintf('Loading Data...')); end
    ncol = twix.image.NCol;
    nlin = twix.image.NLin;
    nslices = twix.image.NSli;
    npartitions = twix.image.NPar;
    naverages = twix.image.NAve;
    nechoes = twix.image.NEco;
    nrepetitions = twix.image.NRep;
    nchannels = twix.image.NCha;    
    fovPE = twix.hdr.Config.PhaseFoV; % mm
    fovRO = twix.hdr.Config.ReadFoV; % mm
    thickness = twix.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness / npartitions; % mm
    os = twix.hdr.Dicom.flReadoutOSFactor; % readout oversample factor
    dwelltime = twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9; % s
    orient = cell2mat(fieldnames(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal));
    frequency = twix.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1}.lFrequency; % Hz
    larmor = 1e-6 * frequency / twix.hdr.Dicom.flMagneticFieldStrength; % MHz/T 
    protocol = twix.hdr.Config.ProtocolName; % name of imaging protocol
    TR = twix.hdr.MeasYaps.alTR{1} / 1000; % ms
    TE = twix.hdr.MeasYaps.alTE{1} / 1000; % ms
    FA = twix.hdr.MeasYaps.adFlipAngleDegree{1}; % degrees
    if parseSpecial
        special = parseheader(twix, protocol);
    else
        special = struct;
    end
    
    % get reference voltage
    if isfield(twix.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1},'flReferenceAmplitude')
        voltage = twix.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1}.flReferenceAmplitude;
    else
        voltage = 0;
    end
    
    % scale voltage
    if isfield(special,'txfactor')
        voltage = voltage * special.txfactor;
        special = rmfield(special,'txfactor');
    end
    
    % determine 2D or 3D
    if npartitions > 1
        mode = '3D';
    else
        mode = '2D';
    end
    
    % determine spiral or cartesian
    if contains(char(aFileNames(jj)),'spiral')
        traj = 'Spiral';
        nleaves = nlin;
        nsamples = ncol;
        if contains(char(aFileNames(jj)),'fancy')
            sptraj = 'fancy';
            mode = '3D';

            if isfield(special,'MatSize') && ~isempty(special.MatSize)
                imgsize = special.MatSize;
            else
                imgsize = 80; % default imgsize of 80
            end

            thickness = fovPE / imgsize; 
            
            if contains(spcalibfile,'3D')
                calibfile = spcalibfile;
            else
                calibfile = sp3Dcalibfile;
            end
        else
            sptraj = 'var';
            calibfile = spcalibfile;
        end
    elseif contains(char(aFileNames(jj)),'radial')
        traj = 'Radial';
        sptraj = 'NA';
    else
        traj = 'Cartesian';
        sptraj = 'NA';
    end
    
    % determine orientation for correctly rotating images
    flip = 0;
    if strcmp(orient,'dTra')
        if strcmp(traj,'Cartesian')
            rot = 1;
        else
            rot = 3;
        end
        orientation = 'Axial';
    elseif strcmp(orient,'dCor')
        if strcmp(traj,'Cartesian')
            rot = 2;
        else
            rot = 3;
            flip = 1;
        end
        orientation = 'Coronal';
    elseif strcmp(orient,'dSag')
        if strcmp(traj,'Cartesian')
            rot = 2;
        else
            rot = 1;
        end
        orientation = 'Sagittal';
    else
        disp('Could not find orientation... Defaulted to axial');
        if strcmp(traj,'Cartesian')
            rot = 1;
        else
            rot = 3;
        end
        orientation = 'Axial';
    end
    
    if strcmpi(sptraj,'fancy') && strcmpi(mode,'3D')
        flip = 0;
        rot = 0;
    end
    
    % save parameters to output structure
    twixs.(['file',num2str(jj)]) = twix;
    allimages(jj).frequency = frequency;
    allimages(jj).FA = FA;
    allimages(jj).nCol = ncol;
    allimages(jj).nLin = nlin;
    allimages(jj).nSli = nslices;
    allimages(jj).nPar = npartitions;
    allimages(jj).nAve = naverages;
    allimages(jj).nRep = nrepetitions;
    allimages(jj).nCha = nchannels;
    allimages(jj).TR = TR;
    allimages(jj).TE = TE;
    allimages(jj).fovRO = fovRO;
    allimages(jj).fovPE = fovPE;
    allimages(jj).thickness = thickness;
    allimages(jj).dwelltime = dwelltime;
    allimages(jj).protocol = protocol;
    allimages(jj).mode = mode;
    allimages(jj).orientation = orientation;
    allimages(jj).trajectory = traj;
    allimages(jj).voltage = voltage;
    fnames = fieldnames(special);
    for i = 1:length(fnames)
       allimages(jj).(fnames{i}) = special.(fnames{i}); 
    end
    
    % reconstruct images   
    if recon 
        if strcmp(traj,'Cartesian')
            % update progressbar
            if progressbar, waitbar((0.05+enum-1)/N,h,'Reconstructing...0%'); end
            
            % allocate memory
            rawimg = zeros(ncol,nlin,max([nslices,npartitions]),nchannels,nrepetitions);
            rawkspace = zeros(ncol,nlin,max([nslices,npartitions]),nchannels,nrepetitions);
            imgcount = 1;
            
            % loop through raw data
            for m = 1:nrepetitions
                for k = 1:nchannels
                    if npartitions > 1
                        %                    Col Cha Lin Par Sli Ave Phs Eco Rep Set Seg
                        rawdata = twix.image(  :, k,  :,  :,  1,  1,  1,  1,  m,  1,  1);
                        rawkspace(:,:,:,k,m) = squeeze(rawdata);
                        rawimg(:,:,:,k,m) = fftshift(fftn(fftshift(rawkspace(:,:,:,k,m))));
                        if progressbar
                            percent = round(imgcount/(nrepetitions*nchannels),3);
                            totpercent = ((percent*0.9)+0.05+enum-1)/N;
                            waitbar(totpercent,h,sprintf('Reconstructing...%3.0f%%',percent*100));
                            imgcount = imgcount + 1;
                        end
                    else
                        for j = 1:nslices
                            data = 0;
                            for n = 1:naverages
                                %                    Col Cha Lin Par Sli Ave Phs Eco Rep Set Seg
                                rawdata = twix.image(  :, k,  :,  :,  j,  n,  1,  1,  m,  1,  1);
                                data = data + rawdata;
                            end
                            rawimg(:,:,j,k,m) = fftshift(fftn(fftshift(squeeze(data))));
                            rawkspace(:,:,j,k,m) = squeeze(data);
                            if progressbar
                                percent = round(imgcount/(nrepetitions*nchannels*nslices),3);
                                totpercent = ((percent*0.9)+0.05+enum-1)/N;
                                waitbar(totpercent,h,sprintf('Reconstructing...%3.0f%%',percent*100));
                                imgcount = imgcount + 1;
                            end
                        end
                    end
                end
            end
            
            
            % update progressbar
            if progressbar, waitbar((0.95+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',95)); end
            
            % combine channels
            [img,kspace] = combinecoils_fa(rawimg,rawkspace);
            
            % update progressbar
            if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end
            
            
        elseif strcmp(traj,'Spiral')
            % calculate or load trajectory
            if strcmp(sptraj,'var')
                if loadtraj
                    [KSpaceCoor,KSpaceDiameter,imgsize] = loadtrajectory(larmor,fovPE,nsamples,nleaves,orientation,[],calibfile);
                    if isempty(imgsize)
                        imgsize = KSpaceDiameter;
                    end
                end
                if ~loadtraj || isempty(KSpaceCoor)
                    [KSpaceCoor,KSpaceDiameter] = calcspiralgrad_20190918(nsamples,nleaves,fovPE,orientation,larmor,sptraj);
                    imgsize = KSpaceDiameter;
                end
                imgsize = floor(imgsize);
                
                % initialize/calculate density correction
                if isempty(dcf) || dcf == 0
                    wi = [];
                elseif dcf == 1
                    gx = zeros(nsamples*nleaves,1);
                    gy = zeros(nsamples*nleaves,1);
                    for i = 1:nleaves
                        gx(1+(i-1)*nsamples:(i*nsamples),1) = diff([0; KSpaceCoor(1+(i-1)*nsamples:(i*nsamples),1)]);
                        gy(1+(i-1)*nsamples:(i*nsamples),1) = diff([0; KSpaceCoor(1+(i-1)*nsamples:(i*nsamples),2)]);
                    end
                    gx = gx ./ (larmor * 100 .* 1e-6);
                    gy = gy ./ (larmor * 100 .* 1e-6);
                    wi = sqrt(gx.^2 + gy.^2) .* sin(atan2(gy,gx) - atan2(KSpaceCoor(:,2),KSpaceCoor(:,1)));
                elseif dcf == 2
                    [~,ia,ic] = unique(KSpaceCoor,'rows','stable');
                    wi = zeros(nsamples*nleaves,1);
                    wi(ia,1) = density_comp_voronoi(KSpaceCoor(ia,:));
                    wi(~ismember(ic,ia),1) = wi(find(~ismember(ic,ia))-1,1);
                elseif dcf == 3
                    wi = zeros(nsamples*nleaves,1);
                    for i = 1:nleaves
                        wi(1+(i-1)*nsamples:(i*nsamples),1) = density_comp_voronoi(KSpaceCoor(1+(i-1)*nsamples:(i*nsamples),:));
                    end
                elseif dcf == 4
                    wi = ones(size(KSpaceCoor,1));
                else
                    wi = [];
                end
                
            elseif strcmp(sptraj,'fancy')
                if isfield(special,'MatSize') && ~isempty(special.MatSize)
                    imgsize = special.MatSize;
                else
                    imgsize = 80;
                end
                if isfield(special,'GAperiod') && ~isempty(special.GAperiod)
                    nreps = special.GAperiod;
                else
                    nreps = 32;
                end
                wi = [];
                if loadtraj
                    KSpaceCoor = loadtrajectory3D(larmor,fovPE,nsamples,nleaves,nreps,imgsize,orientation,[],[],calibfile);
                end
                allimages(jj).thickness = fovPE / imgsize;
            end
            
            % get raw data and reorder                       Col Lin Par Sli Cha Ave Phs Eco Rep Set Seg ...
            rawdata = permute(twix.image(:,:,:,:,:,:,:,:,:,:,:),[1 3 4 5 2 6 7 8 9 10 11 12 13 14 15 16]);
            
            % sum averages
            rawdata = sum(rawdata,6);
            
            % reconstruct based on protocol
            if contains(protocol, 'fa_spiral_dyn')
                % reshape raw data
                rawdata = reshape(rawdata,nlin*ncol,max([nslices,npartitions]),nchannels,nechoes,nrepetitions);
                
                % store recon params
                reco_params = struct;
                reco_params.imgsize = imgsize;
                reco_params.fovPE = fovPE;
                reco_params.zfill = zfill;
                reco_params.GridOSFactor = GridOSFactor;
                reco_params.KernelSize = KernelSize;
                reco_params.savekspace = savekspace;
                reco_params.saveraw = saveraw;
                
                % add any additional parameters
                if contains(protocol, 'fancy')
                    reco_params.nreps = nreps;
                end
                
                % get echo data
                if nechoes == 1
                    rawdata = reshape(rawdata,nlin*ncol,max([nslices,npartitions]),nchannels,nrepetitions);
                elseif nechoes > 1
                    if strcmpi(special.acqOrder,'Radial-Interleaved') 
                        [~,reco_params.VD13_Data] = mapVBVDExtractRawData_fa(char(aFileNames(jj)));
                    end
                end
                
                % store/update progressbar params
                if progressbar
                    waitbar((0.05+enum-1)/N,h,'Reconstructing...0%');
                    progbar = struct;
                    progbar.h = h;
                    progbar.enum = enum;
                    progbar.N = N;
                else
                    progbar = [];
                end
                
                % reconstruct images and kspace
                [img, kspace, rawimg, rawkspace] = fa_spiral_dyn_recon(rawdata,KSpaceCoor,wi,allimages(jj),reco_params,progbar);
                
            else
                % update progressbar
                if progressbar, waitbar((0.05+enum-1)/N,h,'Reconstructing...0%'); end
                
                % reshape raw data
                rawdata = reshape(rawdata,nlin*ncol,max([nslices,npartitions]),nchannels,nrepetitions);
                
                % reconstruct images and kspace
                [rawimg, rawkspace] = gridrecon_fa_20220520(KSpaceCoor,rawdata,imgsize,fovPE, ...
                    'zfill',zfill,'os',GridOSFactor,'k',KernelSize,'wi',wi, ...
                    'verbose',verbose,'filename',gridfile); 
                
                % 3D FFT
                if npartitions > 1
                    rawimg = fftshift(fft(rawimg,[],3),3);
                end
                
                % update progressbar
                if progressbar, waitbar((0.95+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',95)); end
                
                % combine channels
                [img,kspace] = combinecoils_fa(rawimg,rawkspace);
                
                % update progressbar
                if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end
                
            end
        end
        toc;
        
        % save memory if possible
        if ~savekspace
            kspace = [];
            rawkspace = [];
        end
        if ~saveraw
            rawimg = [];
            rawkspace = [];
        end
        
        % shift images
        try
            img = circshift(img,shifts(jj),1);
        catch
            img = circshift(img,shifts(1),1);
        end
        
        % threshold
        if ~ischar(threshold)
            try
                img = threshold_faraz(img,threshold(jj));
                
            catch
                img = threshold_faraz(img,threshold(1));
            end
        else
            if ~strcmpi(threshold,'none')
                img = threshold_faraz(img,threshold);
            end
        end
        
        % make projection
        if projection
            img = sum(img,3);
        end
        
        % crop oversampled regions
        if os > 1
            img_cropped = zeros(ncol/os,nlin,size(img,3),size(img,4));
            crop = floor(1+size(img,1)/(2*os)):ceil(size(img,1)-(size(img,1)/(2*os)));
            for i = 1:size(img,4)
                for j = 1:size(img,3)
                    img_cropped(:,:,j,i) = img(crop,:,j,i);
                end
            end
            img = img_cropped;
        end
        
        % rotate images and kspace
        img = rot90(img,rot);
        kspace = rot90(kspace,rot);
        
        % flip images and kspace
        if flip
            img = fliplr(img);
            kspace = fliplr(kspace);
        end
        
        % convert img to grayscale
        if normalize
            img = mat2gray(img);
        end
        
        % save images
        allimages(jj).images = img;  
        
        % save kspace
        if savekspace
            allimages(jj).kspace = kspace;
        end
        
        % save raw images
        if saveraw
            allimages(jj).rawimg = rawimg;
            if savekspace
                allimages(jj).rawkspace = rawkspace;
            end
        end

        % display images
        if dispimages
            dispImg(allimages(jj).images,'scale',scale','cmap',cmap,'combine',combine,...
                'cscale',cscale,'rows',nrepetitions,'cols',max([nslices,npartitions]));
        end
        
        % save gif
        if savegif
            gfile = ['temp_',datestr(now,'yyyymmdd'),'_set',num2str(jj)];
            img2gif(img,'filename',gfile,'scale',scale,'cscale',cscale,'cmap',cmap);
        end
    end
    
    % print acquisition parameters
    if verbose
        printParam(allimages(jj));
    end
     
    % increment enumerator
    enum = enum + 1;
end

% close progress bar and files
if progressbar, close(h); end
fclose all;

end
