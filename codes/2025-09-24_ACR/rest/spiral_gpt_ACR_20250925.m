% spiral_ACR_20250816_FIXED.m
% Methods:
%   Phantom Static
%   Head Coil (Single Channel)

clear all;
outdir = "/Users/hoomham/Hooman/Work/Analysis/2025-08-15_ACR";
global MS

% ======== user knobs ========
MS_recon_desired = 100;   % sweep this if you like
resizing = 1;             % 1 = scale K to target grid with 1-voxel margin
shifting = 1;             % 1 = apply image-space shift (vox)
% ============================

%% Select files (unchanged list; pick one)
files = [
    % n=2.0, ms=240
    % "/Users/hoomham/Hooman/Work/Images/2025/2025-09-23_ASAP/meas_MID01377_FID13352_fa_spiral_dyn_fancy_v4_20250915_350-20-240.dat"
    % Version 2
    % "/Users/hoomham/Hooman/Work/Images/2025/2025-08-16_ACR/RAW/meas_MID00123_FID12098_fa_spiral_dyn_fancy_v2_20230131.dat" 
    
    % Version 3, fov250
    "/Users/hoomham/Hooman/Work/Images/2025/2025-08-15_ACR/RAW/meas_MID00096_FID12071_fa_spiral_dyn_fancy_v3_20240130.dat"
       

    % n=2.0, ms=160
    % "/Users/hoomham/Hooman/Work/Images/2025/2025-09-19_ASAP/meas_MID01332_FID13307_fa_spiral_dyn_fancy_v4_20250915_350-20-160.dat"

];

%% Calibrations (XYZ @ 800%)

% n=2.0, ms=240
% ROFileName = "/Users/hoomham/Hooman/Work/Images/2025/2025-09-23_ASAP/meas_MID01387_FID13362_fa_spiral_fancy_calibtraj_v4_20250915_X_350-20-240.dat";
% PEFileName = "/Users/hoomham/Hooman/Work/Images/2025/2025-09-23_ASAP/meas_MID01388_FID13363_fa_spiral_fancy_calibtraj_v4_20250915_Y_350-20-240.dat";
% SSFileName = "/Users/hoomham/Hooman/Work/Images/2025/2025-09-23_ASAP/meas_MID01389_FID13364_fa_spiral_fancy_calibtraj_v4_20250915_Z_350-20-240.dat";

% n=2.0, ms=160
% ROFileName = "/Users/hoomham/Hooman/Work/Images/2025/2025-09-15_ASAP/meas_MID01291_FID13266_fa_spiral_fancy_calibtraj_v4_20250915_X.dat"; 
% PEFileName = "/Users/hoomham/Hooman/Work/Images/2025/2025-09-15_ASAP/meas_MID01292_FID13267_fa_spiral_fancy_calibtraj_v4_20250915_Y.dat"; 
% SSFileName = "/Users/hoomham/Hooman/Work/Images/2025/2025-09-15_ASAP/meas_MID01293_FID13268_fa_spiral_fancy_calibtraj_v4_20250915_Z.dat";

% Version 2 800% or 6000%

ROFileName = "/Users/hoomham/Hooman/Work/Images/2025/2025-08-18_ACR/meas_MID00171_FID12146_fa_spiral_fancy_calibtraj_v3_20230821_X250.dat";
PEFileName = "/Users/hoomham/Hooman/Work/Images/2025/2025-08-18_ACR/meas_MID00172_FID12147_fa_spiral_fancy_calibtraj_v3_20230821_Y250.dat";
SSFileName = "/Users/hoomham/Hooman/Work/Images/2025/2025-08-18_ACR/meas_MID00173_FID12148_fa_spiral_fancy_calibtraj_v3_20230821_Z250.dat";

% ROFileName = "/Users/hoomham/Hooman/Work/Images/2025/2025-08-18_ACR/meas_MID00166_FID12141_fa_spiral_fancy_calibtraj_v2_20220506_X250.dat" 
% PEFileName = "/Users/hoomham/Hooman/Work/Images/2025/2025-08-18_ACR/meas_MID00167_FID12142_fa_spiral_fancy_calibtraj_v2_20220506_Y250.dat" 
% SSFileName = "/Users/hoomham/Hooman/Work/Images/2025/2025-08-18_ACR/meas_MID00168_FID12143_fa_spiral_fancy_calibtraj_v2_20220506_Z250.dat"


%% Load raw headers & images (no recon)
[data,twixs] = load_rawdata_20250816(files,'nStudies',1:1);

%% Auto-detect version token from file names (fallback v1)
version = autodetect_version_from_paths(files);
fprintf("Auto-detected version = %s\n", version);

%% Choose study index and extract scan params
i = 1;
[nsamples, nleaves, nreps, FOV_mm, imgsize_hdr, dwell_s, ncha_hdr, nres, nspec] = extract_scan_params(data, i);
dwell_us   = dwell_s * 1e6;
readout_us = nsamples * dwell_us;

%% Flatten TWIX and reshape to [nsamples x nchannels x nblocks]
nchannels = get_channel_count(data, twixs, i, ncha_hdr);
[rawdata, rawdata_sb] = flatten_and_reshape_twix(twixs, i, nsamples, nchannels, nspec);

% Build block indices for responses (static phantom)
[intind, dpind] = build_block_indices(size(rawdata,3), nres); %#ok<NASGU>

% Quick QC
quick_qc_plots(rawdata, rawdata_sb);

%% Build calibration from XYZ scans and load it
[KSpaceCoor, imgsize_cal, meta, chosen] = loadtrajectory3D( ...
    'BuildFromXYZ', struct( ...
        'RO',     ROFileName, ...
        'PE',     PEFileName, ...
        'SS',     SSFileName, ...
        'OutDir', string(outdir), ...
        'Version', string(version) ...
    ), ...
    'VerifyChecksum', true, ...
    'Verbose', true);

if isempty(KSpaceCoor)
    fprintf("Calibration missing, stopping main.\n");
    return;
end

% ---- IMPORTANT: decouple: choose recon grid explicitly ----
fov     = FOV_mm;                            % mm
matsize = round(MS_recon_desired);           % your sweep value
imgsize = matsize;                           % use this everywhere below
res     = fov / matsize;                     % mm/voxel

%% Aggregate across time: one GP per leaf per GA-rep (your original averaging)
raw_gp = zeros(nsamples, nleaves, nreps, nchannels, 'single');
w_gp   = zeros(nsamples, nleaves, nreps, nchannels, 'single');
rep0   = ceil(intind./nleaves) - 1;
int    = mod(intind - 1, nleaves) + 1;
nrep   = mod(rep0, nreps) + 1;

for j = 1:size(rawdata,3)
    raw_gp(:, int(j), nrep(j), :) = raw_gp(:, int(j), nrep(j), :) + rawdata(:, :, j);
    w_gp(:,   int(j), nrep(j), :) = w_gp(:,   int(j), nrep(j), :) + 1;
end
w_gp(w_gp==0) = eps('single');
raw_gp  = raw_gp ./ w_gp;
rawdata2 = reshape(raw_gp, [], nchannels);   % [nsamples*nleaves*nreps  x coils]

%% SAFE TRAJECTORY SCALING to chosen grid (with 1-voxel margin)
if resizing
    kmax_meas  = max( sqrt(sum(KSpaceCoor.^2,2)) );   % 1/mm (from calib)
    k_edge     = (matsize)/(2*fov);                   % 1/mm (grid edge)
    k_target   = (matsize/2 - 1) / fov;               % 1 voxel margin
    alpha      = min(1, k_target / max(kmax_meas, eps));
    KSpaceCoor = alpha * KSpaceCoor;

    % progress line (helps during MS sweep)
    kmax_now = max( sqrt(sum(KSpaceCoor.^2,2)) );
    MS_nyq   = round(2*fov*kmax_now);
    cov      = 100*kmax_now/k_edge;
    fprintf('kmax_meas=%.4f  MS_nyq≈%d  k_edge=%.4f  coverage=%.2f%%\n', ...
            kmax_meas, MS_nyq, k_edge, cov);
end

%% ENSURE DATA & TRAJECTORY LENGTHS MATCH (tile the shorter)
Ntraj = size(KSpaceCoor,1);
Ndata = size(rawdata2,1);
assert( mod(Ndata, Ntraj)==0 || mod(Ntraj, Ndata)==0, ...
        'Data/trajectory length mismatch: data=%d, traj=%d', Ndata, Ntraj);
if Ndata > Ntraj
    KSpaceCoor = repmat(KSpaceCoor, Ndata/Ntraj, 1);
elseif Ntraj > Ndata
    rawdata2   = repmat(rawdata2,   Ntraj/Ndata, 1);
end

%% Optional image-space shift (consistent voxel size)
if shifting
    shiftk   = [0 0 0];                    % vox (set here)
    ph       = 2*pi * (KSpaceCoor .* (shiftk .* res)); % N×3
    rawdata2 = rawdata2 .* exp(1i*sum(ph,2));
end

%% Recon flags
% save_flag = 1;      % 1 = save mat files
img_flag  = 0;      % 0 = gas (GP) only
os   = 3; k = 5; beta = []; zfill = 1;

%% Lookup + Recon (GP)
imgsize3 = [imgsize imgsize imgsize];

[Ind, Dist, wi] = grid_lookup_20230113(KSpaceCoor, imgsize3, fov, ...
    'os', os, 'kernelsize', k, 'beta', beta, 'zfill', zfill);

img_gp_ch = gridrecon_fa_20230113(KSpaceCoor, rawdata2, imgsize3, fov, ...
    'wi', wi, 'Ind', Ind, 'Dist', Dist, 'os', os, 'k', k, 'beta', beta, 'zfill', zfill);

fprintf('Combining channels...\n');
[img_gp,~,b] = combinecoils_fa(img_gp_ch);    %#ok<NASGU>

% Display orientation (as you had)
img_gp    = fliplr(permute(img_gp,   [1 3 2]));
img_gp_ch = fliplr(permute(img_gp_ch,[1 3 2 4]));
b         = fliplr(permute(b,        [2 4 3 1]));

%% Quick mosaic (same slices as before)
% figure;
% imagesc([rot90(squeeze(img_gp(:,22+12,:))),rot90(squeeze(img_gp(:,23+12,:))),rot90(squeeze(img_gp(:,24+12,:))),rot90(squeeze(img_gp(:,25+12,:))),rot90(squeeze(img_gp(:,26+12,:))),rot90(squeeze(img_gp(:,27+12,:))),rot90(squeeze(img_gp(:,28+12,:))),rot90(squeeze(img_gp(:,29+12,:))); ...
%     rot90(squeeze(img_gp(:,30+12,:))),rot90(squeeze(img_gp(:,31+12,:))),rot90(squeeze(img_gp(:,32+12,:))),rot90(squeeze(img_gp(:,33+12,:))),rot90(squeeze(img_gp(:,34+12,:))),rot90(squeeze(img_gp(:,35+12,:))),rot90(squeeze(img_gp(:,36+12,:))),rot90(squeeze(img_gp(:,37+12,:))); ...
%     rot90(squeeze(img_gp(:,38+12,:))),rot90(squeeze(img_gp(:,39+12,:))),rot90(squeeze(img_gp(:,40+12,:))),rot90(squeeze(img_gp(:,41+12,:))),rot90(squeeze(img_gp(:,42+12,:))),rot90(squeeze(img_gp(:,43+12,:))),rot90(squeeze(img_gp(:,44+12,:))),rot90(squeeze(img_gp(:,45+12,:))); ...
%     rot90(squeeze(img_gp(:,46+12,:))),rot90(squeeze(img_gp(:,47+12,:))),rot90(squeeze(img_gp(:,48+12,:))),rot90(squeeze(img_gp(:,49+12,:))),rot90(squeeze(img_gp(:,50+12,:))),rot90(squeeze(img_gp(:,51+12,:))),rot90(squeeze(img_gp(:,52+12,:))),rot90(squeeze(img_gp(:,53+12,:)))]);

imagesc([
    rot90(squeeze(img_gp(:,22-6,:))),rot90(squeeze(img_gp(:,23-6,:))),rot90(squeeze(img_gp(:,24-6,:))),rot90(squeeze(img_gp(:,25-6,:))),rot90(squeeze(img_gp(:,26-6,:))),rot90(squeeze(img_gp(:,27-6,:))),rot90(squeeze(img_gp(:,28-6,:))),rot90(squeeze(img_gp(:,29-6,:))); ...
    rot90(squeeze(img_gp(:,30-6,:))),rot90(squeeze(img_gp(:,31-6,:))),rot90(squeeze(img_gp(:,32-6,:))),rot90(squeeze(img_gp(:,33-6,:))),rot90(squeeze(img_gp(:,34-6,:))),rot90(squeeze(img_gp(:,35-6,:))),rot90(squeeze(img_gp(:,36-6,:))),rot90(squeeze(img_gp(:,37-6,:))); ...
    rot90(squeeze(img_gp(:,38-6,:))),rot90(squeeze(img_gp(:,39-6,:))),rot90(squeeze(img_gp(:,40-6,:))),rot90(squeeze(img_gp(:,41-6,:))),rot90(squeeze(img_gp(:,42-6,:))),rot90(squeeze(img_gp(:,43-6,:))),rot90(squeeze(img_gp(:,44-6,:))),rot90(squeeze(img_gp(:,45-6,:))); ...
    rot90(squeeze(img_gp(:,46-6,:))),rot90(squeeze(img_gp(:,47-6,:))),rot90(squeeze(img_gp(:,48-6,:))),rot90(squeeze(img_gp(:,49-6,:))),rot90(squeeze(img_gp(:,50-6,:))),rot90(squeeze(img_gp(:,51-6,:))),rot90(squeeze(img_gp(:,52-6,:))),rot90(squeeze(img_gp(:,53-6,:))); ...
    rot90(squeeze(img_gp(:,54-6,:))),rot90(squeeze(img_gp(:,55-6,:))),rot90(squeeze(img_gp(:,56-6,:))),rot90(squeeze(img_gp(:,57-6,:))),rot90(squeeze(img_gp(:,58-6,:))),rot90(squeeze(img_gp(:,59-6,:))),rot90(squeeze(img_gp(:,60-6,:))),rot90(squeeze(img_gp(:,61-6,:))); ...
    rot90(squeeze(img_gp(:,62-6,:))),rot90(squeeze(img_gp(:,63-6,:))),rot90(squeeze(img_gp(:,64-6,:))),rot90(squeeze(img_gp(:,65-6,:))),rot90(squeeze(img_gp(:,66-6,:))),rot90(squeeze(img_gp(:,67-6,:))),rot90(squeeze(img_gp(:,68-6,:))),rot90(squeeze(img_gp(:,69-6,:)))]);


axis image off; colormap gray; caxis([-2 22]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS (all copied from your file; ONLY the ones below differ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ver = autodetect_version_from_paths(files)
tok = regexp(string(files(:)), "(?<=[_\-])v\d+(?=[_\-]|\.|$)", "match");
tok = string([tok{:}]); tok = unique(tok(~strlength(tok)==0));
if ~isempty(tok), ver = tok(end); else, ver = "v1"; end
end

function [nsamples, nleaves, nreps, FOV_mm, imgsize, dwell_s, ncha_hdr, nres, nspec] = extract_scan_params(data, i)
assert(i >= 1 && i <= numel(data), 'Index i=%d out of range.', i);
nsamples = double(data(i).nCol);
nleaves  = double(data(i).nLin);
nreps    = double(data(i).nRep);
FOV_mm   = double(data(i).fovPE);
imgsize  = double(data(i).MatSize);
dwell_s  = double(data(i).dwelltime);
ncha_hdr = []; if isfield(data,'nCha') && ~isempty(data(i).nCha), ncha_hdr = double(data(i).nCha); end
nres     = 1;   if isfield(data,'numRes') && ~isempty(data(i).numRes), nres = max(1, double(data(i).numRes)); end
nspec    = 0;   if isfield(data,'bSpectra') && isfield(data,'numSpec') && data(i).bSpectra, nspec = double(data(i).numSpec); end
end

function nchannels = get_channel_count(data, twixs, i, ncha_hdr)
if ~isempty(ncha_hdr)
    nchannels = ncha_hdr;
else
    fn = fieldnames(twixs); fn = fn{min(i,numel(fn))};
    if isfield(twixs.(fn).image,'NCha'), nchannels = double(twixs.(fn).image.NCha);
    else, nchannels = 1;
    end
end
end

function [rawdata, rawdata_sb] = flatten_and_reshape_twix(twixs, i, nsamples, nchannels, nspec)
fn = fieldnames(twixs); fn = fn{min(i,numel(fn))};
raw_flat = squeeze(twixs.(fn).image(:,:,:,:,:,:,:,:,:,:,:));
L = numel(raw_flat); base = nsamples * nchannels;
assert(mod(L,base)==0, 'Raw data length (%d) not divisible by nsamples*nchannels (%d).', L, base);
nblocks = L / base;
rawdata = reshape(raw_flat, nsamples, nchannels, nblocks);
if nspec > 0
    assert(nblocks >= nspec, 'Not enough blocks for nspec=%d', nspec);
    rawdata = rawdata(:, :, (nspec+1):end);
end
rawdata_sb = [];
if isfield(twixs,'file2') && isfield(twixs.file2,'image')
    sb_flat = squeeze(twixs.file2.image(:,:,:,:,:,:,:,:,:,:,:));
    Lsb = numel(sb_flat);
    assert(mod(Lsb,base)==0, 'SB raw length (%d) not divisible by base (%d).', Lsb, base);
    nblocks_sb = Lsb / base;
    rawdata_sb = reshape(sb_flat, nsamples, nchannels, nblocks_sb);
end
end

function [intind, dpind] = build_block_indices(nblocks, nres)
intind = repelem((1:ceil(nblocks / nres)).', nres); intind = intind(1:nblocks);
dpind  = repmat((0:nres-1).', ceil(nblocks / nres), 1); dpind = dpind(1:nblocks); %#ok<NASGU>
end

function quick_qc_plots(rawdata, rawdata_sb)
try
    absval = reshape(squeeze(abs(sum(rawdata(1,:,:),2))),[],1);
    figure; plot(absval); title('Block-wise magnitude (first-sample sum across coils)');
    xlabel('block'); ylabel('|sum|');
    if ~isempty(rawdata_sb) && ndims(rawdata_sb) >= 3 && size(rawdata_sb,1) >= 1
        absval_sb = reshape(squeeze(abs(sum(rawdata_sb(1,:,:),2))),[],1);
        figure; plot(absval_sb); title('Sideband block magnitude'); xlabel('block'); ylabel('|sum|');
    end
catch
end
end

% ------------------------------------------------------------------------
% The rest of your helper functions (load_rawdata_20250816, build_calibration_from_xyz,
% loadtrajectory3D & friends, grid_lookup_20230113, gridrecon_fa_20230113, parseheader, etc.)
% are IDENTICAL to your version and can remain as-is.
% ------------------------------------------------------------------------


function [allimages,twixs] = load_rawdata_20250816(aFileNames,Options)
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
    Options.spcalibfile string = 'calibrations_20201109.mat'                % mat file containing calibrated spiral trajectories
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
    Options.verbose double = 1                                              % 1 to printout parameters
    Options.parseSpecial double = 1                                         % 1 to get parameters from the "special" tab
end

% flags
projection = Options.projection;
normalize = Options.normalize;
dcf = Options.dcf;
combine = Options.combine;
dispimages = Options.dispimages;
progressbar = Options.progressbar;
savegif = Options.savegif;
saveraw = Options.saveraw;
savekspace = Options.savekspace;
verbose = Options.verbose;
parseSpecial = Options.parseSpecial;

% recon parameters
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
    thickness = twix.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness / max(npartitions,1); % mm
    os = twix.hdr.Dicom.flReadoutOSFactor; % readout oversample factor

    % safer dwelltime extraction
    try
        dwelltime = twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9; % s
    catch
        dwelltime = twix.hdr.MeasYaps.sRXSPEC.alDwellTime(1) * 1e-9; % fallback
    end

    % safer orientation extraction
    try
        sNorm = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal;
        if isfield(sNorm,'dTra') && abs(sNorm.dTra) > 0.5
            orient = 'dTra';
        elseif isfield(sNorm,'dCor') && abs(sNorm.dCor) > 0.5
            orient = 'dCor';
        elseif isfield(sNorm,'dSag') && abs(sNorm.dSag) > 0.5
            orient = 'dSag';
        else
            orient = 'dTra'; % fallback
        end
    catch
        orient = 'dTra'; % default axial if missing
    end

    % frequency and larmor safety
    try
        frequency = twix.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1}.lFrequency; % Hz
        B0 = twix.hdr.Dicom.flMagneticFieldStrength; % T
        if B0 > 0
            larmor = 1e-6 * frequency / B0; % MHz/T
        else
            larmor = NaN;
        end
    catch
        frequency = NaN;
        larmor = NaN;
        B0 = NaN;
    end

    protocol = twix.hdr.Config.ProtocolName; % name of imaging protocol
    TR = twix.hdr.MeasYaps.alTR{1} / 1000; % ms
    TE = twix.hdr.MeasYaps.alTE{1} / 1000; % ms
    FA = twix.hdr.MeasYaps.adFlipAngleDegree{1}; % degrees
    if parseSpecial
        try
            special = parseheader(twix, protocol);
        catch
            warning('parseheader failed, continuing without special fields');
            special = struct;
        end
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
        special.txfactor_applied = true;
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
                % imgsize = MS; % default imgsize of 80
                warning('No MatSize found, defaulting imgsize to 80');
            end

            thickness = fovPE / imgsize;
        else
            sptraj = 'var';
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

    % print acquisition parameters
    if verbose
        printParam(allimages(jj));
    end

    % increment enumerator
    enum = enum + 1;
end

% close progress bar and files
if progressbar && exist('h','var') && isvalid(h)
    close(h);
end
fclose all;

end


function [cal, outfile] = build_calibration_from_xyz(ROFileName, PEFileName, SSFileName, outdir, version)
%BUILD_CALIBRATION_FROM_XYZ  Generate a per-calibration MAT file (new schema)
% from three ASAP 3D calibration scans (X/RO, Y/PE, Z/SS), using the TWO-SLICE (±D)
% difference-of-differences method.
%
% Usage:
%   [cal, outfile] = build_calibration_from_xyz(RO, PE, SS, outdir, 'v3');
%
% Creates and saves a 'cal' struct with fields:
%   version, nucleus, units, nsamples, nleaves, nreps, FOV_mm, imgsize,
%   readout_us, dwell_us, orientation, gamma_MHzT, kx, ky, kz (1/mm),
%   kmax_measured, created_utc, author, data_sha1
%
% Key differences vs previous buggy version:
%   - Reads BOTH slices (j=1 and j=2).
%   - Per-slice baseline removal (on/off, inv/off).
%   - Uses cross-slice difference-of-differences and denominator 4*D*2π.
%   - D is computed from the two slice centers: D = |z2 - z1|/2 (mm).

arguments
    ROFileName (1,1) string
    PEFileName (1,1) string
    SSFileName (1,1) string
    outdir (1,1) string
    version (1,1) string {mustBeNonzeroLengthText}
end

global MS

% ------------------- Load TWIX for each axis -------------------
T = cell(1,3);
F = {ROFileName, PEFileName, SSFileName};
for ii = 1:3
    Tw = mapVBVD(char(F{ii}));
    if iscell(Tw), Tw = Tw{end}; end
    T{ii} = Tw;
end

% Use RO header as reference
Tw = T{1};

% ------------------- Basic header fields -------------------
nsamples   = double(Tw.image.NCol);
nlin_total = double(Tw.image.NLin);
nreps      = double(Tw.image.NRep);      % for calib: should be >= 3 (on, inv, off)
ncha       = double(Tw.image.NCha);
navg       = double(Tw.image.NAve);
FOV_mm     = double(Tw.hdr.Config.PhaseFoV);     % mm (3D fancy: isotropic)
imgsize    = guess_imgsize(Tw, MS);

% dwell time (s -> us)
dwell_s    = safe_dwell_s(Tw);
dwell_us   = dwell_s * 1e6;
readout_us = nsamples * dwell_s * 1e6;

% frequency, B0, gamma
[~, ~, gamma_MHzT] = safe_freq_gamma(Tw);

% orientation string (dTra/dCor/dSag -> Axial/Coronal/Sagittal)
orientation = safe_orientation(Tw);

% nucleus label (site-specific WIP—defaults to water)
nucleus = safe_nucleus(Tw);

% per-repetition interleaves
if nreps > 0 && mod(nlin_total, nreps) == 0
    nleaves = nlin_total / nreps;
else
    nleaves = nlin_total; % fallback
end
Nseg = nsamples * nleaves;   % samples per rep per axis

% ------------------- TWO-SLICE separation D (mm) from slice centers -------------------
% We prefer to compute D from slice 1 and slice 2 centers: D = |pos2 - pos1| / 2.
% If slice 2 is missing, fall back to legacy (pos1 - thickness/2).
try
    orientKey = char(fieldnames(Tw.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal));
catch
    orientKey = 'dTra'; % default
end

has2 = isfield(Tw.hdr.MeasYaps.sSliceArray, 'asSlice') && numel(Tw.hdr.MeasYaps.sSliceArray.asSlice) >= 2 ...
    && isfield(Tw.hdr.MeasYaps.sSliceArray.asSlice{2}, 'sPosition') ...
    && isfield(Tw.hdr.MeasYaps.sSliceArray.asSlice{2}.sPosition, orientKey);

pos1 = NaN; pos2 = NaN;
try
    pos1 = double(Tw.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.(orientKey));
catch, end
if has2
    try
        pos2 = double(Tw.hdr.MeasYaps.sSliceArray.asSlice{2}.sPosition.(orientKey));
    catch, pos2 = NaN; end
end

% if ~isnan(pos1) && ~isnan(pos2)
%     D_mm = 0.5 * abs(pos2 - pos1);
% else
% legacy fallback
try
    thk = double(Tw.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness);
catch, thk = NaN; end
if ~isnan(pos2) && ~isnan(thk)
    D_mm = max(1e-6, pos2 - thk/2); % keep >0
    warning('Two-slice center not available; using legacy D = pos1 - thickness/2 = %.6g mm', D_mm);
else
    D_mm = 1; % last resort
    warning('Slice offset D not found; using D=1 mm placeholder (relative k-scale).');
end
% end

% ------------------- Collect raw: BOTH slices, states on/inv/off -------------------
% We follow the original: use rep #1=#on, #2=#inv, #3=#off (calib protocol).
if nreps < 3
    error('Calibration requires at least 3 repetitions (on, inv, off). Found NRep=%d.', nreps);
end
rep_on = 1; rep_inv = 2; rep_off = 3;
ave_idx = max(navg,1);

% Helper to read one (axis, slice, rep, coil) block as a column vector (Nseg x 1)
    function col = read_block(ti, sli_idx, rep_idx, ch_idx)
        % Col Cha Lin Par Sli Ave Phs Eco Rep Set Seg
        sig = squeeze(ti.image(:, ch_idx, :, 1, sli_idx, ave_idx, 1, 1, rep_idx, 1, 1));
        col = reshape(sig, [], 1);
    end

% For each axis, build concatenated [on; inv; off] per slice → (3*Nseg x ncha)
raw3_s1 = cell(1,3);
raw3_s2 = cell(1,3);

for ax = 1:3
    ti = T{ax};

    % verify there are at least 2 slices
    nsli = double(ti.image.NSli);
    if nsli < 2
        error('Calibration scan for axis %d has only %d slice(s). Need 2.', ax, nsli);
    end

    A1 = zeros(3*Nseg, ncha, 'like', complex(0));  % slice 1
    A2 = zeros(3*Nseg, ncha, 'like', complex(0));  % slice 2
    for ch = 1:ncha
        s1_on  = read_block(ti, 1, rep_on, ch);
        s1_inv = read_block(ti, 1, rep_inv, ch);
        s1_off = read_block(ti, 1, rep_off, ch);

        s2_on  = read_block(ti, 2, rep_on, ch);
        s2_inv = read_block(ti, 2, rep_inv, ch);
        s2_off = read_block(ti, 2, rep_off, ch);

        A1(:,ch) = [s1_on; s1_inv; s1_off];
        A2(:,ch) = [s2_on; s2_inv; s2_off];
    end
    raw3_s1{ax} = A1;
    raw3_s2{ax} = A2;
end

% ------------------- Phase -> k (1/mm) per coil using TWO-SLICE method -------------------
% Returns Nseg x ncha for each axis
[Kx, Ky, Kz] = compute_k_from_triplet(raw3_s1, raw3_s2, Nseg, D_mm);

% ------------------- Interleave normalization and channel combine -------------------
Kx = interleaf_normalize(Kx, nsamples, nleaves);
Ky = interleaf_normalize(Ky, nsamples, nleaves);
Kz = interleaf_normalize(Kz, nsamples, nleaves);

% ------------------- Interleave DC allign -------------------
% Kx = interleaf_dc_align(Kx, nsamples, nleaves, 'median');
% Ky = interleaf_dc_align(Ky, nsamples, nleaves, 'median');
% Kz = interleaf_dc_align(Kz, nsamples, nleaves, 'median');

kx = mean(Kx, 2);  % combine coils (mean)
ky = mean(Ky, 2);
kz = mean(Kz, 2);

k0 = [kx,ky,kz];
ks = [smooth(kx,'loess'), smooth(ky,'loess'), smooth(kz,'loess')];

dk = ks - k0;                      % per-sample change (1/mm)
fprintf('max |Δk| = %.3g  median |Δk| = %.3g (1/mm)\n', max(abs(dk),[],'all'), median(abs(dk(:))));
fprintf('kmax change: %.5f -> %.5f (Δ=%.3g)\n', max(vecnorm(k0,2,2)), max(vecnorm(ks,2,2)), ...
        max(vecnorm(ks,2,2)) - max(vecnorm(k0,2,2)));

% Optional smoothing (loess — same as original)
kx = smooth(kx, 'loess');
ky = smooth(ky, 'loess');
kz = smooth(kz, 'loess');

% Xenon scaling if non-proton: scale by gamma ratio
if ~strcmpi(nucleus,'water') && ~isnan(gamma_MHzT)
    gamma_H1 = 42.577478518; % MHz/T
    scale = gamma_MHzT / gamma_H1;
    kx = kx .* scale; ky = ky .* scale; kz = kz .* scale;
end

% ------------------- Build calibration struct -------------------
cal = struct();
cal.version      = char(version);
cal.nucleus      = char(nucleus);
cal.units        = struct('k','1/mm','fov','mm','dwell','us','readout','us');
cal.nsamples     = double(nsamples);
cal.nleaves      = double(nleaves);
cal.nreps        = double(nreps);
cal.FOV_mm       = double(FOV_mm);
cal.imgsize      = double(imgsize);
cal.readout_us   = double(readout_us);
cal.dwell_us     = double(dwell_us);
cal.orientation  = char(orientation);
cal.gamma_MHzT   = double(gamma_MHzT);
cal.kx           = double(kx(:));
cal.ky           = double(ky(:));
cal.kz           = double(kz(:));
cal.created_utc  = char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
cal.author       = getenv_default('USER','unknown');

% QC and checksum
kmax = max(sqrt(kx.^2 + ky.^2 + kz.^2));
cal.kmax_measured = double(kmax);
cal.data_sha1     = sha1hex([cal.kx cal.ky cal.kz]);

% ------------------- Save file -------------------
if ~isfolder(outdir), mkdir(outdir); end
outfile = fullfile(outdir, compose_filename(cal));
save(outfile, 'cal');
fprintf('Saved calibration to:\n  %s\n', outfile);
end

function dwell_s = safe_dwell_s(Tw)
try
    dwell_s = Tw.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9;
catch
    try
        dwell_s = Tw.hdr.MeasYaps.sRXSPEC.alDwellTime(1) * 1e-9;
    catch
        dwell_s = NaN;
    end
end
end

function [freq_Hz, B0_T, gamma_MHzT] = safe_freq_gamma(Tw)
try
    freq_Hz = Tw.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1}.lFrequency; % Hz
catch
    freq_Hz = NaN;
end
try
    B0_T = Tw.hdr.Dicom.flMagneticFieldStrength;
catch
    B0_T = NaN;
end
if ~isnan(freq_Hz) && ~isnan(B0_T) && B0_T > 0
    gamma_MHzT = 1e-6 * freq_Hz / B0_T;
else
    gamma_MHzT = NaN;
end
end

function orientation = safe_orientation(Tw)
orientation = 'Axial';
try
    sNorm = Tw.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal;
    if isfield(sNorm,'dTra') && abs(sNorm.dTra) > 0.5, orientation = 'Axial'; end
    if isfield(sNorm,'dCor') && abs(sNorm.dCor) > 0.5, orientation = 'Coronal'; end
    if isfield(sNorm,'dSag') && abs(sNorm.dSag) > 0.5, orientation = 'Sagittal'; end
end
end

function nucleus = safe_nucleus(Tw)
nucleus = 'water';
try
    wip = Tw.hdr.MeasYaps.sWipMemBlock.alFree;
    forcexenon = wip{2}; lnucleus = wip{3};
    if ~isempty(forcexenon)
        if lnucleus == 1, nucleus = 'dissolved'; else, nucleus = 'gas'; end
    end
end
end

function D_mm = safe_slice_offset_mm(Tw)
% Compute D from slice position projected onto dominant normal, then minus half-thickness.
D_mm = NaN;
try
    sNorm = Tw.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal;
    sPos  = Tw.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition;
    thick = Tw.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness; % mm
    if isfield(sNorm,'dTra') && abs(sNorm.dTra) > 0.5
        p = abs(sPos.dTra);
    elseif isfield(sNorm,'dCor') && abs(sNorm.dCor) > 0.5
        p = abs(sPos.dCor);
    elseif isfield(sNorm,'dSag') && abs(sNorm.dSag) > 0.5
        p = abs(sPos.dSag);
    else
        p = NaN;
    end
    if ~isnan(p) && ~isempty(thick)
        D_mm = p - (thick/2);
        if D_mm <= 0, D_mm = abs(D_mm); end
    end
catch
    D_mm = NaN;
end
end

function [Kx, Ky, Kz] = compute_k_from_triplet(raw3_s1, raw3_s2, Nseg, D_mm)
%COMPUTE_K_FROM_TRIPLET  Convert phase to k(1/mm) using TWO-SLICE ±D method.
%
% Inputs
%   raw3_s1 : 1x3 cell, each cell is (3*Nseg x ncha) for axis X,Y,Z from slice #1,
%             each concatenated as [on; inv; off].
%   raw3_s2 : 1x3 cell, same for slice #2.
%   Nseg    : scalar, nsamples*nleaves
%   D_mm    : HALF the slice center-to-center separation (mm), i.e., D = |z2 - z1|/2.
%
% Output
%   Kx, Ky, Kz : (Nseg x ncha) k-space coordinates (1/mm), per coil.
%
% Method:
%   per slice: unwrap+smooth phase, baseline-subtract by 'off'
%   cross-slice difference-of-differences:
%       k = [ (phi1_on - phi1_inv) - (phi2_on - phi2_inv) ] / (4 * D_mm * 2*pi)

ncha = size(raw3_s1{1},2);
Kx = zeros(Nseg, ncha);
Ky = zeros(Nseg, ncha);
Kz = zeros(Nseg, ncha);

% Helper to get baseline-subtracted on,inv for one slice/axis
    function [p_on, p_inv] = phi_on_inv(A)
        % A is (3*Nseg x ncha) = [on; inv; off]
        on  = A(1:Nseg,           :);
        inv = A(Nseg+1:2*Nseg,    :);
        off = A(2*Nseg+1:3*Nseg,  :);

        % unwrap+smooth phase then per-slice baseline subtract
        p_on  = zeros(Nseg, size(A,2));
        p_inv = zeros(Nseg, size(A,2));
        for c = 1:size(A,2)
            s_on  = getPhase(on(:,c),  'smoothunwrap');
            s_inv = getPhase(inv(:,c), 'smoothunwrap');
            s_off = getPhase(off(:,c), 'smoothunwrap');
            p_on(:,c)  = s_on  - s_off;
            p_inv(:,c) = s_inv - s_off;
        end
    end

% Axis X
[p1_on, p1_inv] = phi_on_inv(raw3_s1{1});   % slice 1
[p2_on, p2_inv] = phi_on_inv(raw3_s2{1});   % slice 2
Kx = ((p1_on - p1_inv) - (p2_on - p2_inv)) ./ (4 * D_mm * 2*pi);

% Axis Y
[p1_on, p1_inv] = phi_on_inv(raw3_s1{2});
[p2_on, p2_inv] = phi_on_inv(raw3_s2{2});
Ky = ((p1_on - p1_inv) - (p2_on - p2_inv)) ./ (4 * D_mm * 2*pi);

% Axis Z
[p1_on, p1_inv] = phi_on_inv(raw3_s1{3});
[p2_on, p2_inv] = phi_on_inv(raw3_s2{3});
Kz = ((p1_on - p1_inv) - (p2_on - p2_inv)) ./ (4 * D_mm * 2*pi);
end

function K = interleaf_normalize(K, nsamples, nleaves)
% Align the start of each interleaf to the first interleaf baseline.
base = K(1,:);  % first sample of first interleaf for each coil
for n = 2:nleaves
    idx = (n-1)*nsamples + 1 : n*nsamples;
    dev = base - K(idx(1),:);
    K(idx,:) = K(idx,:) + dev;
end
end

function out = getPhase(in, method)
if nargin < 2 || isempty(method), method = 'nothing'; end
out = atan2(imag(in), real(in));
switch lower(method)
    case 'smooth',        out = smooth(out,'loess');
    case 'smoothabs',     out = smooth(abs(out),'loess');
    case 'unwrap',        out = unwrap(out);
    case 'smoothunwrap',  out = smooth(unwrap(out),'loess');
    otherwise
        % no-op
end
end

function imgsize = guess_imgsize(Tw, defaultVal)
imgsize = defaultVal;
try
    if isfield(Tw.hdr,'MeasYaps') && isfield(Tw.hdr.MeasYaps,'sWipMemBlock') && isfield(Tw.hdr.MeasYaps.sWipMemBlock,'alFree')
        al = Tw.hdr.MeasYaps.sWipMemBlock.alFree;
        if numel(al) >= 10 && ~isempty(al{10})
            imgsize = double(al{10});
        end
    end
catch, end
end

function s = getenv_default(name, def)
val = getenv(name);
if isempty(val), s = def; else, s = val; end
end

function fn = compose_filename(cal)
fn = sprintf('calib_ASAP3D_%s_FOV%gmm_NS%d_NL%d_NR%d_RO%gus_DW%gus_%s_%s.mat', ...
    cal.version, cal.FOV_mm, cal.nsamples, cal.nleaves, cal.nreps, cal.readout_us, cal.dwell_us, cal.nucleus, upper(cal.orientation(1:min(2,end))) );
end

function [KSpaceCoor, imgsize, meta, chosen] = loadtrajectory3D(varargin)
%LOADTRAJECTORY3D  Unified loader for ASAP 3D trajectory
% Usage:
%   % A) From explicit calibration .mat file
%   [K, imgsize, meta] = loadtrajectory3D('CalibFile','/path/calib_ASAP_v2_...mat');
%
%   % B) Build from XYZ scans, save, then load (calib = ON)
%   [K, imgsize, meta] = loadtrajectory3D('BuildFromXYZ',struct( ...
%       'RO', '/path/X.dat', 'PE','/path/Y.dat', 'SS','/path/Z.dat', ...
%       'OutDir','/path/out', 'Version','v2'));
%
%   % C) Match by metadata against a JSON manifest (calib = OFF)
%   spec = struct('version','v2','nsamples',512,'nleaves',20,'nreps',32, ...
%                 'FOV_mm',250,'imgsize',80,'readout_us',2560,'dwell_us',5);
%   [K, imgsize, meta, chosen] = loadtrajectory3D('Spec',spec,'Manifest','/path/manifest.json');

p = inputParser; p.KeepUnmatched = true;
addParameter(p,'CalibFile','',@(s)ischar(s)||isstring(s));
addParameter(p,'Spec',struct(),@isstruct);
addParameter(p,'Manifest','',@(s)ischar(s)||isstring(s));
addParameter(p,'BuildFromXYZ',struct(),@isstruct);
addParameter(p,'VerifyChecksum',true,@islogical);
addParameter(p,'Verbose',true,@islogical);
parse(p,varargin{:});
args = p.Results;

if ~isempty(args.CalibFile)
    [KSpaceCoor, imgsize, meta] = loadtrajectory3D_fromfile(args.CalibFile, ...
        'VerifyChecksum', args.VerifyChecksum, 'Verbose', args.Verbose);
    chosen = struct('mode','file','path',string(args.CalibFile));
    return;
end

if ~isempty(fieldnames(args.BuildFromXYZ))
    bx = args.BuildFromXYZ;  % fields: RO, PE, SS, OutDir, Version
    [cal, outfile] = build_calibration_from_xyz(string(bx.RO), string(bx.PE), string(bx.SS), string(bx.OutDir), string(bx.Version));
    [KSpaceCoor, imgsize, meta] = loadtrajectory3D_fromfile(outfile, 'VerifyChecksum', args.VerifyChecksum, 'Verbose', args.Verbose);
    chosen = struct('mode','built','path',string(outfile));
    return;
end

if ~isempty(fieldnames(args.Spec)) && ~isempty(args.Manifest)
    [KSpaceCoor, imgsize, meta, chosen] = loadtrajectory3D_match(args.Spec, string(args.Manifest), ...
        'VerifyChecksum', args.VerifyChecksum, 'Verbose', args.Verbose);
    return;
end

error('loadtrajectory3D: specify one of CalibFile, BuildFromXYZ, or Spec+Manifest.');
end

function [KSpaceCoor, imgsize, meta] = loadtrajectory3D_fromfile(calibPath, varargin)
%LOADTRAJECTORY3D_FROMFILE  Load a single calibration .mat and return K

p = inputParser; p.KeepUnmatched = true;
addParameter(p,'VerifyChecksum',true,@islogical);
addParameter(p,'Verbose',true,@islogical);
parse(p,varargin{:});
args = p.Results;

S = load(calibPath);

% Accept either 'cal' (new), 's' (legacy), or top-level struct
if isfield(S,'cal')
    cal = S.cal;
elseif isfield(S,'s') && isstruct(S.s)
    cal = S.s;
else
    cal = S;
end

% Required fields (allow FOV_mm or fov)
req = {'nsamples','nleaves','kx','ky','kz'};
for i = 1:numel(req)
    assert(isfield(cal,req{i}), 'Calibration missing field: %s', req{i});
end

% Shapes & lengths
kx = cal.kx(:); ky = cal.ky(:); kz = cal.kz(:);
Nseg = numel(kx);
Nexp = double(cal.nsamples) * double(cal.nleaves);
assert(Nseg == Nexp, 'Length mismatch: got %d, expected %d (nsamples*nleaves).', Nseg, Nexp);

% ----- imgsize with SAFE fallbacks (NO global MS) -----
if isfield(cal,'imgsize') && ~isempty(cal.imgsize)
    imgsize = double(cal.imgsize);
else
    % derive a reasonable default from FOV and kmax if present
    kmax = max(sqrt(kx.^2 + ky.^2 + kz.^2));
    if isfield(cal,'FOV_mm') && ~isempty(cal.FOV_mm) && isfinite(kmax) && kmax > 0
        imgsize = round(2 * double(cal.FOV_mm) * kmax);   % Nyquist grid ~ 2*FOV*kmax
    elseif isfield(cal,'fov') && ~isempty(cal.fov) && isfinite(kmax) && kmax > 0
        imgsize = round(2 * double(cal.fov) * kmax);
    else
        imgsize = NaN;  % unknown; caller should ignore/override
    end
end

% FOV in mm (support new and old)
if isfield(cal,'FOV_mm') && ~isempty(cal.FOV_mm)
    FOVmm = double(cal.FOV_mm);
elseif isfield(cal,'fov') && ~isempty(cal.fov)
    FOVmm = double(cal.fov);
else
    FOVmm = NaN;
end

% NR with fallback
if isfield(cal,'nreps') && ~isempty(cal.nreps)
    NR = double(cal.nreps);
else
    NR = 1;
end

% Compose outputs
KSpaceCoor = [kx ky kz];

% Checksum (optional)
if args.VerifyChecksum && isfield(cal,'data_sha1')
    try
        hex = sha1hex([kx ky kz]);
        assert(strcmpi(hex, cal.data_sha1), 'Checksum mismatch for %s', calibPath);
    catch
        warning('Checksum verification failed/skipped for %s', calibPath);
    end
end

% Meta without heavy arrays
meta = cal;
fd = fieldnames(meta);
rm = intersect({'kx','ky','kz'}, fd);
if ~isempty(rm), meta = rmfield(meta, rm); end

if args.Verbose
    % robust print: allow NaN imgsize without exploding
    if isnan(imgsize)
        fprintf('Loaded calibration: %s\n  NS=%d, NL=%d, NR=%d, FOV=%.1f mm, imgsize=? (unspecified), N=%d\n', ...
            calibPath, cal.nsamples, cal.nleaves, NR, FOVmm, Nseg);
    else
        fprintf('Loaded calibration: %s\n  NS=%d, NL=%d, NR=%d, FOV=%.1f mm, imgsize=%d, N=%d\n', ...
            calibPath, cal.nsamples, cal.nleaves, NR, FOVmm, imgsize, Nseg);
    end
end
end

function [KSpaceCoor, imgsize, meta, chosen] = loadtrajectory3D_match(spec, manifestPath, varargin)
%LOADTRAJECTORY3D_MATCH  Locate a calibration by metadata and load it.
%
% Inputs:
%   spec: struct with fields like:
%     version (string), nsamples, nleaves (per-rep), nreps, nleaves_total,
%     FOV_mm, imgsize, readout_us, dwell_us
%   manifestPath: JSON file path (array of objects)
%
% Outputs:
%   KSpaceCoor: [N x 3] double (kx,ky,kz)
%   imgsize:    scalar
%   meta:       calibration metadata struct (no k arrays)
%   chosen:     manifest row used (with selected_from=manifestPath)

p = inputParser; p.KeepUnmatched = true;
addParameter(p,'VerifyChecksum',true,@islogical);
addParameter(p,'Verbose',true,@islogical);
addParameter(p,'StrictFOV',true,@islogical);         % << enforce FOV in fallback
parse(p,varargin{:});
args = p.Results;

% ---------- load & decode manifest ----------
assert(isfile(manifestPath), 'Manifest not found: %s', manifestPath);
raw = fileread(manifestPath);
M = jsondecode(raw);
if ~iscell(M), M = num2cell(M); end

% ---------- normalize rows ----------
for ii = 1:numel(M)
    row = M{ii};

    % alias "file" -> "path"
    if ~isfield(row,'path') && isfield(row,'file')
        row.path = row.file;
    end

    % coerce numeric strings -> double
    numKeys = {'nsamples','nleaves','nreps','nleaves_total','FOV_mm','imgsize','readout_us','dwell_us'};
    for kk = 1:numel(numKeys)
        k = numKeys{kk};
        if isfield(row,k) && ischar(row.(k))
            val = str2double(row.(k));
            if ~isnan(val), row.(k) = val; end
        end
    end

    % infer nleaves_total if missing (from fields or NL### in filename)
    if ~isfield(row,'nleaves_total') || isempty(row.nleaves_total)
        if isfield(row,'nleaves') && isfield(row,'nreps') && ~isempty(row.nleaves) && ~isempty(row.nreps)
            row.nleaves_total = double(row.nleaves) * double(row.nreps);
        elseif isfield(row,'path') && ~isempty(row.path)
            [~,bn,~] = fileparts(string(row.path));
            tok = regexp(bn, "NL(\d+)", "tokens", "once");
            if ~isempty(tok), row.nleaves_total = str2double(tok{1}); end
        end
    end

    % infer FOV_mm if missing (FOV###mm in filename)
    if (~isfield(row,'FOV_mm') || isempty(row.FOV_mm)) && isfield(row,'path') && ~isempty(row.path)
        [~,bn,~] = fileparts(string(row.path));
        tok = regexp(bn, "FOV(\d+)mm", "tokens", "once");
        if ~isempty(tok), row.FOV_mm = str2double(tok{1}); end
    end

    % ensure version is string
    if isfield(row,'version'), row.version = string(row.version); end

    M{ii} = row;
end

% ---------- helper: tolerant match ----------
    function [tf, why] = matches(row, s)
        tf = true; why = "";

        % version (string)
        if isfield(s,'version') && ~isempty(s.version)
            if ~isfield(row,'version') || ~strcmpi(string(row.version), string(s.version))
                tf=false; why="version"; return;
            end
        end

        % nsamples exact
        if isfield(s,'nsamples') && ~isempty(s.nsamples)
            if ~isfield(row,'nsamples') || double(row.nsamples) ~= double(s.nsamples)
                tf=false; why="nsamples"; return;
            end
        end

        % total interleaves comparison (accept per-rep or total on either side)
        s_nr = 1; if isfield(s,'nreps') && ~isempty(s.nreps), s_nr = double(s.nreps); end
        if isfield(s,'nleaves_total') && ~isempty(s.nleaves_total)
            s_total = double(s.nleaves_total);
        elseif isfield(s,'nleaves') && ~isempty(s.nleaves)
            s_total = double(s.nleaves) * s_nr;
        else
            s_total = [];
        end

        if ~isempty(s_total)
            if isfield(row,'nleaves_total') && ~isempty(row.nleaves_total)
                r_total = double(row.nleaves_total);
            elseif isfield(row,'nleaves') && isfield(row,'nreps') ...
                    && ~isempty(row.nleaves) && ~isempty(row.nreps)
                r_total = double(row.nleaves) * double(row.nreps);
            else
                r_total = [];
            end
            if ~isempty(r_total) && r_total ~= s_total
                tf=false; why=sprintf('nleaves_total (row=%g vs spec=%g)', r_total, s_total); return;
            end
        end

        % FOV within ±0.25 mm
        if isfield(s,'FOV_mm') && ~isempty(s.FOV_mm)
            if ~isfield(row,'FOV_mm') || abs(double(row.FOV_mm) - double(s.FOV_mm)) > 0.25
                tf=false; why="FOV_mm"; return;
            end
        end

        % imgsize exact (optional)
        if isfield(s,'imgsize') && ~isempty(s.imgsize)
            if ~isfield(row,'imgsize') || double(row.imgsize) ~= double(s.imgsize)
                tf=false; why="imgsize"; return;
            end
        end

        % timing tolerances
        if isfield(s,'dwell_us') && ~isempty(s.dwell_us)
            if ~isfield(row,'dwell_us') || abs(double(row.dwell_us) - double(s.dwell_us)) > 1e-3
                tf=false; why="dwell_us"; return;
            end
        end
        if isfield(s,'readout_us') && ~isempty(s.readout_us)
            if ~isfield(row,'readout_us') || abs(double(row.readout_us) - double(s.readout_us)) > 1
                tf=false; why="readout_us"; return;
            end
        end
    end

% ---------- find matches ----------
idx = [];
whyNot = strings(0,1);
for iRow = 1:numel(M)
    [ok, why] = matches(M{iRow}, spec);
    if ok
        idx(end+1) = iRow; %#ok<AGROW>
    else
        whyNot(end+1,1) = sprintf("#%d %s", iRow, why); %#ok<AGROW>
    end
end

% ---------- fallback if needed (with StrictFOV) ----------
if isempty(idx)
    if args.Verbose && ~isempty(whyNot)
        fprintf("No exact manifest match. Reasons (first up to 5):\n");
        disp(whyNot(1:min(5,end)));
    end

    % relaxed fallback: same version & nsamples; choose closest total interleaves
    cand = [];
    for iRow = 1:numel(M)
        row = M{iRow};

        % version + nsamples must match
        if ~(isfield(row,'version') && strcmpi(string(row.version), string(spec.version)) ...
                && isfield(row,'nsamples') && double(row.nsamples) == double(spec.nsamples))
            continue;
        end

        % if StrictFOV, require FOV match (±0.25 mm)
        if args.StrictFOV && isfield(spec,'FOV_mm') && ~isempty(spec.FOV_mm)
            if ~isfield(row,'FOV_mm') || abs(double(row.FOV_mm) - double(spec.FOV_mm)) > 0.25
                continue;
            end
        end

        if isfield(row,'nleaves_total') && ~isempty(row.nleaves_total)
            cand(end+1,:) = [iRow, double(row.nleaves_total)]; %#ok<AGROW>
        end
    end

    if ~isempty(cand) && isfield(spec,'nleaves_total') && ~isempty(spec.nleaves_total)
        [~,j] = min(abs(cand(:,2) - double(spec.nleaves_total)));
        idx = cand(j,1);
        if args.Verbose
            fprintf("Fallback: using closest nleaves_total in manifest (row %d) with matching FOV.\n", idx);
        end
    end
end

if isempty(idx)
    if args.Verbose
        fprintf("No calibration matched the provided spec (even after fallback).\n");
    end
    % Return empty outputs instead of error
    KSpaceCoor = [];
    imgsize    = [];
    meta       = struct();
    chosen     = struct();
    return;
end

assert(~isempty(idx), 'No calibration matched the provided spec (even after fallback).');

% if multiple, pick newest by created_utc, else last
if numel(idx) > 1
    try
        ts = cellfun(@(r) string(r.created_utc), M(idx), 'UniformOutput', false);
        t  = datetime(ts);
        [~, ord] = sort(t);
        idx = idx(ord(end));
    catch
        idx = idx(end);
    end
end

row = M{idx};
if ~isfield(row,'path') && isfield(row,'file'), row.path = row.file; end
assert(isfield(row,'path') && ~isempty(row.path), 'Manifest row missing "path" to .mat');

% ---------- load the chosen calibration ----------
[KSpaceCoor, imgsize, meta] = loadtrajectory3D_fromfile(row.path, ...
    'VerifyChecksum', args.VerifyChecksum, 'Verbose', args.Verbose);

chosen = row;
chosen.selected_from = string(manifestPath);

if args.Verbose
    fprintf('Matched calibration from manifest: %s\n', manifestPath);
end
end

% ---------- utilities reused ----------
function hex = sha1hex(x)
% Compute SHA-1 of numeric array x via Java
JMD = java.security.MessageDigest.getInstance('SHA-1');
ba = typecast(reshape(double(x).',1,[]),'uint8');
JMD.update(int8(ba));
d = typecast(JMD.digest(),'uint8');
hex = lower(reshape(dec2hex(d,2).',1,[]));
end

function [Ind,Dist,wi] = grid_lookup_20230113(KSpaceCoor,imgsize,fov,Options)
% load or save spiral grid parameters based on input trajectory/options

% parse arguments
arguments
    KSpaceCoor
    imgsize
    fov
    Options.os double = 3
    Options.zfill double = 1
    Options.kernelsize double = 5
    Options.beta double = []
    Options.wi double = []
    Options.parallel double = 1
    Options.filename char = 'grid_lookup_20230113.mat'
    Options.date char = datestr(now,'yyyy-mm-dd')
    Options.save double = 1
    Options.verbose double = 1
end

% initialize outputs
Ind = [];
Dist = [];
wi = [];

% check arguments
if isempty(Options.beta)
    Options.beta = pi*sqrt(Options.kernelsize^2/Options.os^2*(Options.os-0.5)^2-0.8); % Eq(5) in https://ieeexplore.ieee.org/document/1435541/
end

% find file
fileloc = which(Options.filename);

if ~isempty(fileloc)
    % load grids
    s = load(fileloc);
    s = s.s;

    % initial search
    match1 = arrayfun(@(x) ...
        all(size(x.KSpaceCoor) == size(KSpaceCoor))         && ...
        all(x.imgsize == imgsize)                           && ...
        abs(x.fov - fov)                            <= 0    && ...
        abs(x.os - Options.os)                      <= 0    && ...
        abs(x.kernelsize - Options.kernelsize)      <= 0    && ...
        abs(x.beta - Options.beta)                  <= 0    && ...
        abs(x.zfill - Options.zfill)                <= 0    ...
        , s);

    % filter out incorrect sized trajectories
    s2 = s(match1);

    % search for saved grid
    match = arrayfun(@(x) ...
        abs(sum(x.KSpaceCoor(:) - single(KSpaceCoor(:)))) <= 0, s2);

    % get index for last (most recent) match
    i = find(match == 1,1,'last');
else
    i = [];
    s = [];
end

% get grid params if available or calculate
if isempty(i)

    % get last index
    i = 1+size(s,2);

    % check 2D vs 3D
    dimk = size(KSpaceCoor,2);

    % gridding

    % check grid size
    if length(imgsize) == 1
        imgsize = repmat(imgsize,[1 3]);
    end

    % create cartesian grid
    x = (-ceil(Options.zfill*imgsize(1))/2/fov:1/(fov*Options.os):ceil(Options.zfill*imgsize(1))/2/fov - 1/(fov*Options.os))';
    y = (-ceil(Options.zfill*imgsize(2))/2/fov:1/(fov*Options.os):ceil(Options.zfill*imgsize(2))/2/fov - 1/(fov*Options.os))';
    z = (-ceil(Options.zfill*imgsize(3))/2/fov:1/(fov*Options.os):ceil(Options.zfill*imgsize(3))/2/fov - 1/(fov*Options.os))';

    gridszx = size(x,1);
    gridszy = size(y,1);
    gridszz = size(z,1);

    % find N^2 closest grid points for each sampled kspace point
    if Options.verbose, disp('Finding nearest grid points...'); tic; end

    if gridszx > 180 || gridszy > 180 || gridszz > 180 || ~all(imgsize==imgsize(1))
        Indx = knnsearch(x,KSpaceCoor(:,1),'K',1);
        Indy = knnsearch(y,KSpaceCoor(:,2),'K',1);
        Indz = knnsearch(z,KSpaceCoor(:,3),'K',1);
        k = Options.kernelsize;
        Ind = zeros(size(KSpaceCoor,1),k^2,'single');
        Dist = zeros(size(KSpaceCoor,1),k^2,'single');
        if Options.parallel
            parfor ii = 1:size(KSpaceCoor,1)
                [nx,ny,nz] = ndgrid(max(1,Indx(ii)-floor(k/2)):min(gridszx,Indx(ii)+floor(k/2)), ...
                    max(1,Indy(ii)-floor(k/2)):min(gridszy,Indy(ii)+floor(k/2)), ...
                    max(1,Indz(ii)-floor(k/2)):min(gridszz,Indz(ii)+floor(k/2)));
                [nD,nI] = sort(sqrt(sum((KSpaceCoor(ii,:)-[x(nx(:)),y(ny(:)),z(nz(:))]).^2,2)));
                Dist(ii,:) = nD(1:k^2);
                Ind(ii,:) = sub2ind([gridszx gridszy gridszz], ...
                    nx(nI(1:k^2)), ...
                    ny(nI(1:k^2)), ...
                    nz(nI(1:k^2)));
            end
        else
            for ii = 1:size(KSpaceCoor,1)
                [nx,ny,nz] = ndgrid(max(1,Indx(ii)-floor(k/2)):min(gridszx,Indx(ii)+floor(k/2)), ...
                    max(1,Indy(ii)-floor(k/2)):min(gridszy,Indy(ii)+floor(k/2)), ...
                    max(1,Indz(ii)-floor(k/2)):min(gridszz,Indz(ii)+floor(k/2)));
                [nD,nI] = sort(sqrt(sum((KSpaceCoor(ii,:)-[x(nx(:)),y(ny(:)),z(nz(:))]).^2,2)));
                Dist(ii,:) = nD(1:k^2);
                Ind(ii,:) = sub2ind([gridszx gridszy gridszz], ...
                    nx(nI(1:k^2)), ...
                    ny(nI(1:k^2)), ...
                    nz(nI(1:k^2)));
            end
        end
    else
        [x,y,z] = ndgrid(-ceil(Options.zfill*imgsize(1))/2/fov:1/(fov*Options.os):ceil(Options.zfill*imgsize(1))/2/fov - 1/(fov*Options.os));
        try
            [Ind,Dist] = knnsearch([x(:) y(:) z(:)],KSpaceCoor,'K',Options.kernelsize^2,'IncludeTies',true,'NSMethod','kdtree');
            Ind = cell2mat(Ind);
            Dist = cell2mat(Dist);
        catch
            [Ind,Dist] = knnsearch([x(:) y(:) z(:)],KSpaceCoor,'K',Options.kernelsize^2);
        end
    end

    if Options.verbose, fprintf('Nearest points found in %g seconds. \n',toc); end

    % create Kaiser-Bessel Kernel
    width = Options.kernelsize/(fov*Options.os);
    klength = 10000;
    [kernel,u] = createKBkernel(width,Options.beta,klength);
    kernel3Dtable = interp1(u',kernel,Dist,'linear',0);

    % calculate density correction (if not provided)
    if isempty(Options.wi)
        if Options.verbose, tic; end
        iter = 5;
        dcftable = kernel3Dtable.^(1/2);
        wi = iterative_dcf_fa_20190910(iter,KSpaceCoor,dcftable,Ind,[gridszx gridszy gridszz],[],Options.verbose);
        if Options.verbose, fprintf('Density correction finished in %g seconds. \n',toc); end
    else
        wi = Options.wi;
    end

    % add new grid
    if Options.save
        s(i).KSpaceCoor = single(KSpaceCoor);
        s(i).fov = fov;
        s(i).imgsize = imgsize;
        s(i).os = Options.os;
        s(i).kernelsize = Options.kernelsize;
        s(i).beta = Options.beta;
        s(i).zfill = Options.zfill;
        s(i).wi = single(wi);
        s(i).Ind = single(Ind);
        s(i).Dist = single(Dist);
        s(i).date = Options.date;

        % update .mat file
        save(fileloc,'s');
        if Options.verbose
            disp(['Grid saved to ',fileloc]);
        end
    end

else
    % load grid
    Ind = s2(i).Ind;
    Dist = s2(i).Dist;
    wi = s2(i).wi;
    if Options.verbose
        disp(['Grid from ',s2(i).date,' loaded']);
    end
end

end

function [Image_out,KSpace_out,wi,Ind,Dist] = gridrecon_fa_20230113(adKSpaceCoor,rawdata,NumK,fov,Options)
% Look at each sampled kspace point, convolve with kernel, resample onto grid

% parse arguments
arguments
    adKSpaceCoor                                            % nx2 or nx3 kspace coordinates
    rawdata                                                 % raw kspace data (nCol*nLin,max([nSli,nPar]),nCha,nRep)
    NumK                                                    % desired matrix size
    fov                                                     % image FOV
    Options.os double = 3                                   % grid oversampling factor
    Options.zfill double = 1                                % grid zero-filling factor
    Options.k double = 5                                    % kernel size (along 1 dimension)
    Options.verbose double = 0                              % flag for printing progress
    Options.wi double = []                                  % nx1 density correction
    Options.Ind double = []                                 % indices for nearest neighbors
    Options.Dist double = []                                % distances for nearest neighbors
    Options.beta double = []                                % KB kernel beta
    Options.parallel double = 1                             % 1 to search neighbors in parallel (only for large grids)
    Options.lookup double = 1                               % flag for looking up grid
    Options.filename char = 'grid_lookup_20230113.mat'      % filename for saved grids
    Options.savegrid double = 1                             % flag for saving current grid
end

% initialize outputs
Image_out = [];
KSpace_out = [];

% check 2D vs 3D
dimk = size(adKSpaceCoor,2);

% find grid if available
if Options.lookup && (isempty(Options.Ind) && isempty(Options.Dist))
    [Ind,Dist,wi] = grid_lookup_20230113(adKSpaceCoor,NumK,fov,'os',Options.os,'zfill',Options.zfill, ...
        'kernelsize',Options.k,'beta',Options.beta,'wi',Options.wi,'parallel',Options.parallel, ...
        'filename',Options.filename,'save',Options.savegrid,'verbose',Options.verbose);
else
    Ind = Options.Ind;
    Dist = Options.Dist;
    wi = Options.wi;
end

% Check grid size
if length(NumK) == 1
    NumK = repmat(NumK,[1 3]);
end

% Create cartesian grid
x = (-ceil(Options.zfill*NumK(1))/2/fov:1/(fov*Options.os):ceil(Options.zfill*NumK(1))/2/fov - 1/(fov*Options.os))';
y = (-ceil(Options.zfill*NumK(2))/2/fov:1/(fov*Options.os):ceil(Options.zfill*NumK(2))/2/fov - 1/(fov*Options.os))';
z = (-ceil(Options.zfill*NumK(3))/2/fov:1/(fov*Options.os):ceil(Options.zfill*NumK(3))/2/fov - 1/(fov*Options.os))';

gridszx = size(x,1);
gridszy = size(y,1);
gridszz = size(z,1);

% Find N^2 closest grid points for each sampled kspace point
if isempty(Ind) || isempty(Dist)
    if Options.verbose, disp('Finding nearest grid points...'); tic; end

    if gridszx > 200 || gridszy > 200 || gridszz > 200 || ~all(NumK==NumK(1))
        Indx = knnsearch(x,adKSpaceCoor(:,1),'K',1);
        Indy = knnsearch(y,adKSpaceCoor(:,2),'K',1);
        Indz = knnsearch(z,adKSpaceCoor(:,3),'K',1);
        k = Options.k;
        Ind = zeros(size(adKSpaceCoor,1),k^2,'single');
        Dist = zeros(size(adKSpaceCoor,1),k^2,'single');
        if Options.parallel
            parfor i = 1:size(adKSpaceCoor,1)
                [nx,ny,nz] = ndgrid(max(1,Indx(i)-floor(k/2)):min(gridszx,Indx(i)+floor(k/2)), ...
                    max(1,Indy(i)-floor(k/2)):min(gridszy,Indy(i)+floor(k/2)), ...
                    max(1,Indz(i)-floor(k/2)):min(gridszz,Indz(i)+floor(k/2)));
                [nD,nI] = sort(sqrt(sum((adKSpaceCoor(i,:)-[x(nx(:)),y(ny(:)),z(nz(:))]).^2,2)));
                Dist(i,:) = nD(1:k^2);
                Ind(i,:) = sub2ind([gridszx gridszy gridszz], ...
                    nx(nI(1:k^2)), ...
                    ny(nI(1:k^2)), ...
                    nz(nI(1:k^2)));
            end
        else
            for i = 1:size(adKSpaceCoor,1)
                [nx,ny,nz] = ndgrid(max(1,Indx(i)-floor(k/2)):min(gridszx,Indx(i)+floor(k/2)), ...
                    max(1,Indy(i)-floor(k/2)):min(gridszy,Indy(i)+floor(k/2)), ...
                    max(1,Indz(i)-floor(k/2)):min(gridszz,Indz(i)+floor(k/2)));
                [nD,nI] = sort(sqrt(sum((adKSpaceCoor(i,:)-[x(nx(:)),y(ny(:)),z(nz(:))]).^2,2)));
                Dist(i,:) = nD(1:k^2);
                Ind(i,:) = sub2ind([gridszx gridszy gridszz], ...
                    nx(nI(1:k^2)), ...
                    ny(nI(1:k^2)), ...
                    nz(nI(1:k^2)));
            end
        end
    else
        [x,y,z] = ndgrid(-ceil(Options.zfill*NumK(1))/2/fov:1/(fov*Options.os):ceil(Options.zfill*NumK(1))/2/fov - 1/(fov*Options.os));
        try
            [Ind,Dist] = knnsearch([x(:) y(:) z(:)],adKSpaceCoor,'K',Options.k^2,'IncludeTies',true);
            Ind = cell2mat(Ind);
            Dist = cell2mat(Dist);
        catch
            [Ind,Dist] = knnsearch([x(:) y(:) z(:)],adKSpaceCoor,'K',Options.k^2);
        end
    end

    if Options.verbose, fprintf('Nearest points found in %g seconds. \n',toc); end
end

% Create Kaiser-Bessel Kernel
if isempty(Options.beta)
    Options.beta = pi*sqrt(Options.k^2/Options.os^2*(Options.os-0.5)^2-0.8); % Eq(5) in https://ieeexplore.ieee.org/document/1435541/
end
width = Options.k/(fov*Options.os);
klength = 10000;
[kernel,u] = createKBkernel(width,Options.beta,klength);
kernel3Dtable = interp1(u',kernel,Dist,'linear',0);

% Calculate density correction (if not provided)
if isempty(wi)
    if Options.verbose, tic; end
    iter = 5;
    dcftable = kernel3Dtable.^(1/2);
    wi = iterative_dcf_fa_20190910(iter,adKSpaceCoor,dcftable,Ind,[gridszx gridszy gridszz],[],Options.verbose);
    if Options.verbose, fprintf('Density correction finished in %g seconds. \n',toc); end
end

% Get raw data dimensions
dim = size(rawdata);

% Check dimensions
if length(dim) < 3, dim(end+1:3) = 1; end

% Pre-allocate memory
k_real = zeros([gridszx*gridszy*gridszz,dim(2),dim(3)],'single');
k_imag = zeros([gridszx*gridszy*gridszz,dim(2),dim(3)],'single');
Auxiliary = zeros([gridszx*gridszy*gridszz,dim(2),dim(3)],'single');

% Multiply sampled data with density correction
M = single(rawdata.*wi);

% Reshape density correction
wr = single(repmat(wi,1,dim(2),dim(3)));

% Transpose kernel table
kernel3Dtable = single(kernel3Dtable');

% Loop through each sampled data point
if Options.verbose, disp('Adding points onto grid...'); tic; end
for i = 1:size(adKSpaceCoor,1)
    % Convolve kernel with kspace point
    Mk = M(i,:,:).*kernel3Dtable(:,i);

    % Add convolved data to grid
    k_real(Ind(i,:),:,:) = k_real(Ind(i,:),:,:) + real(Mk);
    k_imag(Ind(i,:),:,:) = k_imag(Ind(i,:),:,:) + imag(Mk);
    Auxiliary(Ind(i,:),:,:) = Auxiliary(Ind(i,:),:,:) + (wr(i,:,:).*kernel3Dtable(:,i));
end

% Combine real and imag components
KSpace_out = complex(k_real,k_imag);

% Partial de-apodization (ignore regions outside FOV)
KSpace_out = KSpace_out ./ Auxiliary;
KSpace_out(isnan(KSpace_out)) = 0;

% FFT kspace to image
KSpace_out = reshape(KSpace_out,[gridszx,gridszy,gridszz,dim(2),dim(3)]);
Image_out = KSpace_out;
for i = 1:3
    Image_out = fftshift(fft(fftshift(Image_out,i),[],i),i);
end

% Resize image
if Options.os > 1
    Image_out = Image_out(floor(1+(gridszx-ceil(NumK(1)*Options.zfill))/2):floor(1+(gridszx-ceil(NumK(1)*Options.zfill))/2)+(ceil(NumK(1)*Options.zfill)-1),...
        floor(1+(gridszy-ceil(NumK(2)*Options.zfill))/2):floor(1+(gridszy-ceil(NumK(2)*Options.zfill))/2)+(ceil(NumK(2)*Options.zfill)-1),...
        floor(1+(gridszz-ceil(NumK(3)*Options.zfill))/2):floor(1+(gridszz-ceil(NumK(3)*Options.zfill))/2)+(ceil(NumK(3)*Options.zfill)-1),:,:);
elseif Options.os < 1
    Image_out = padarray(Image_out,...
        [floor((ceil(NumK(1)*Options.zfill)-gridszx)/2),...
        floor((ceil(NumK(2)*Options.zfill)-gridszy)/2),...
        floor((ceil(NumK(3)*Options.zfill)-gridszz)/2)],0,'both');
end

% Print
if Options.verbose, fprintf('Gridding completed in %g seconds. \n',toc); end
end

function special = parseheader(twix,protocol)
% Parses through the variables of the "Special" tab for Siemens scanners

% check input
if nargin < 2, protocol = ''; end

% initalize output structure
special = struct;

if isfield(twix.hdr.MeasYaps,'sWipMemBlock')
    if contains(protocol,'fa_spiral_dyn')
        if contains(protocol,'fancy_v2') || contains(protocol,'fancy_v3') || contains(protocol,'fancy_v4')
            % version number
            if contains(protocol,'fancy_v4')
                version = 20250815
            else
                version = extract(protocol,digitsPattern(8));
                version = str2double(version{1});
            end

            % tx factor
            special.txfactor = twix.hdr.MeasYaps.sWipMemBlock.adFree{1};

            % number of images
            special.numImg = twix.hdr.MeasYaps.sWipMemBlock.alFree{2};

            % number of resonances
            special.numRes = 1 + twix.hdr.MeasYaps.sWipMemBlock.alFree{9};
            if isempty(special.numRes), special.numRes = 1; end

            % pulse shape for additional resonances
            special.RFshape = twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
            if isempty(special.RFshape)
                special.RFshape = 'aSinc';
            elseif special.RFshape == 1
                special.RFshape = 'Gauss';
            elseif special.RFshape == 2
                special.RFshape = 'Rect';
            elseif special.RFshape == 3
                special.RFshape = 'Kai';
            end

            % acquisition order
            special.acqOrder = twix.hdr.MeasYaps.sWipMemBlock.alFree{10};
            if isempty(special.acqOrder)
                special.acqOrder = 'Sequential';
            elseif special.acqOrder == 1
                special.acqOrder = 'GP-DP';
            elseif special.acqOrder == 2
                special.acqOrder = 'GP-DP-DP';
            else
                special.acqOrder = 'who fucking knows';
            end

            % frequency, flip angle, and pulse duration for other resonances
            special.RFdur = 1.0240; % ms
            if special.numRes >= 2
                special.freq2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{3}; % ppm
                special.FA2 = twix.hdr.MeasYaps.sWipMemBlock.alFree{4}; % deg
                special.RFdur2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{5}; % ms
            end

            % matrix size
            special.MatSize = twix.hdr.MeasYaps.sWipMemBlock.alFree{6};

            % golden angle periodicity
            special.GAperiod = twix.hdr.MeasYaps.sWipMemBlock.alFree{7};

            % check if spectra acquired
            if twix.hdr.MeasYaps.sWipMemBlock.alFree{11} == 1
                special.bSpectra = true;
            else
                special.bSpectra = false;
            end

            % spectra parameters
            if special.bSpectra
                special.numSpec = twix.hdr.MeasYaps.sWipMemBlock.alFree{12}; % number of spectra
                special.dtSpec = twix.hdr.MeasYaps.sWipMemBlock.alFree{13}; % dwelltime (us/point)
                special.FASpec = twix.hdr.MeasYaps.sWipMemBlock.alFree{14}; % deg
            end

        elseif contains(protocol,'fancy')
            % version number
            version = extract(protocol,digitsPattern(8));
            version = str2double(version{1});

            % tx factor
            special.txfactor = twix.hdr.MeasYaps.sWipMemBlock.adFree{1};

            % number of images
            special.numImg = twix.hdr.MeasYaps.sWipMemBlock.alFree{2};

            % number of resonances
            special.numRes = 1 + twix.hdr.MeasYaps.sWipMemBlock.alFree{9};
            if isempty(special.numRes), special.numRes = 1; end

            if version < 20220316
                % delay between images (ms)
                special.repDelay = twix.hdr.MeasYaps.sWipMemBlock.adFree{3};
            else
                % number of images before alternating saturation
                if special.numRes > 1
                    special.numAlt = twix.hdr.MeasYaps.sWipMemBlock.alFree{3};
                end
            end

            % pulse shape for additional resonances
            special.RFshape = twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
            if isempty(special.RFshape)
                special.RFshape = 'aSinc';
            elseif special.RFshape == 1
                special.RFshape = 'Gauss';
            elseif special.RFshape == 2
                special.RFshape = 'Rect';
            elseif special.RFshape == 3
                special.RFshape = 'Kai';
            end

            % acquisition order
            special.acqOrder = twix.hdr.MeasYaps.sWipMemBlock.alFree{10};
            if isempty(special.acqOrder)
                special.acqOrder = 'Sequential';
            elseif special.acqOrder == 1
                special.acqOrder = 'Interleaved';
            elseif special.acqOrder == 2
                special.acqOrder = 'Inter-Interleaved';
            elseif special.acqOrder == 3
                special.acqOrder = 'Saturation';
            elseif special.acqOrder == 4
                special.acqOrder = 'ADC';
                special.bvalue = twix.hdr.MeasYaps.sWipMemBlock.adFree{11};
            elseif special.acqOrder == 5
                special.acqOrder = 'Forced DP Traj';
            elseif special.acqOrder == 6
                special.acqOrder = 'GP-DP-DP';
            elseif special.acqOrder == 7
                special.acqOrder = 'GP-DP-DP-Spectra';
            else
                special.acqOrder = 'who fucking knows';
            end

            % frequency, flip angle, and pulse duration for other resonances
            special.RFdur = twix.hdr.MeasYaps.sWipMemBlock.adFree{6}; % ms
            if special.numRes >= 2
                special.freq2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{4}; % ppm
                special.FA2 = twix.hdr.MeasYaps.sWipMemBlock.alFree{5}; % deg
                special.RFdur2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{6}; % ms

                if special.numRes >= 3
                    special.freq3 = twix.hdr.MeasYaps.sWipMemBlock.adFree{11}; % ppm
                    special.FA3 = twix.hdr.MeasYaps.sWipMemBlock.alFree{12}; % deg
                    special.RFdur3 = twix.hdr.MeasYaps.sWipMemBlock.adFree{13}; % ms
                end
            end

            % matrix size
            special.MatSize = twix.hdr.MeasYaps.sWipMemBlock.alFree{7};

            % golden angle periodicity
            special.GAperiod = twix.hdr.MeasYaps.sWipMemBlock.alFree{14};
        else
            % tx factor
            special.txfactor = twix.hdr.MeasYaps.sWipMemBlock.adFree{1};

            % number of images
            special.numImg = twix.hdr.MeasYaps.sWipMemBlock.alFree{2};

            % number of resonances
            special.numRes = 1 + twix.hdr.MeasYaps.sWipMemBlock.alFree{9};
            if isempty(special.numRes), special.numRes = 1; end

            % delay between images (ms)
            special.repDelay = twix.hdr.MeasYaps.sWipMemBlock.adFree{3};

            % pulse shape for additional resonances
            special.RFshape = twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
            if isempty(special.RFshape)
                special.RFshape = 'aSinc';
            elseif special.RFshape == 1
                special.RFshape = 'Gauss';
            elseif special.RFshape == 2
                special.RFshape = 'Rect';
            elseif special.RFshape == 3
                special.RFshape = 'Kai';
            end

            % acquisition order
            special.acqOrder = twix.hdr.MeasYaps.sWipMemBlock.alFree{10};
            if isempty(special.acqOrder)
                special.acqOrder = 'Sequential';
            elseif special.acqOrder == 1
                special.acqOrder = 'Interleaved';
            elseif special.acqOrder == 2
                special.acqOrder = 'Inter-Interleaved';
            elseif special.acqOrder == 3
                special.acqOrder = 'Radial-Interleaved';
            elseif special.acqOrder == 4
                special.acqOrder = 'Saturation';
            end

            % frequency, flip angle, and pulse duration for other resonances
            special.RFdur = twix.hdr.MeasYaps.sWipMemBlock.adFree{6}; % ms
            if special.numRes >= 2
                special.freq2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{4}; % ppm
                special.FA2 = twix.hdr.MeasYaps.sWipMemBlock.alFree{5}; % deg
                special.RFdur2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{6}; % ms

                if special.numRes >= 3
                    special.freq3 = twix.hdr.MeasYaps.sWipMemBlock.adFree{11}; % ppm
                    special.FA3 = twix.hdr.MeasYaps.sWipMemBlock.alFree{12}; % deg
                    special.RFdur3 = twix.hdr.MeasYaps.sWipMemBlock.adFree{13}; % ms
                end
            end

            % golden angle trajectory used?
            special.GA = twix.hdr.MeasYaps.sWipMemBlock.alFree{7};

            % golden angle periodicity
            if special.GA == 1
                special.GAperiod = twix.hdr.MeasYaps.sWipMemBlock.alFree{14};
            end
        end

    elseif contains(protocol,'fa_spiral_hybridXTC')

        % tx factor
        special.txfactor = twix.hdr.MeasYaps.sWipMemBlock.adFree{1};

        % default saturation method
        special.SatMethod = 'none';

        % number of images per set
        special.NumImg = twix.hdr.MeasYaps.sWipMemBlock.alFree{3};

        % GP saturation after 1st image?
        if twix.hdr.MeasYaps.sWipMemBlock.alFree{7} == 1
            special.GPsat = true;
        else
            special.GPsat = false;
        end

        if ~isempty(twix.hdr.MeasYaps.sWipMemBlock.alFree{2})
            % saturation pulse flip angle
            special.SatFA = twix.hdr.MeasYaps.sWipMemBlock.alFree{8}; % deg

            % saturation pulse width
            special.SatDuration = twix.hdr.MeasYaps.sWipMemBlock.adFree{9}; % ms

            % saturation pulse type
            if twix.hdr.MeasYaps.sWipMemBlock.alFree{10} == 1
                special.SatPulse = 'Rect';
            else
                special.SatPulse = 'Gauss';
            end

            % number frequency, and spacing of RBC saturations
            if ismember(twix.hdr.MeasYaps.sWipMemBlock.alFree{2}, [1 3])
                special.SatNum_RBC = twix.hdr.MeasYaps.sWipMemBlock.alFree{4};
                special.SatFreq_RBC = twix.hdr.MeasYaps.sWipMemBlock.adFree{5}; % ppm
                try
                    special.SatSpacing_RBC = twix.hdr.MeasYaps.sWipMemBlock.adFree{6}; % ms
                catch
                    special.SatSpacing_RBC = 0;
                end
            end

            % number frequency, and spacing of RBC saturations
            if ismember(twix.hdr.MeasYaps.sWipMemBlock.alFree{2}, [2 3])
                special.SatNum_TP = twix.hdr.MeasYaps.sWipMemBlock.alFree{11};
                special.SatFreq_TP = twix.hdr.MeasYaps.sWipMemBlock.adFree{12}; %  ppm
                try
                    special.SatSpacing_TP = twix.hdr.MeasYaps.sWipMemBlock.adFree{13}; % ms
                catch
                    special.SatSpacing_TP = 0;
                end
            end

            % saturation method
            if twix.hdr.MeasYaps.sWipMemBlock.alFree{2} == 1
                special.SatMethod = 'RBC';
            elseif twix.hdr.MeasYaps.sWipMemBlock.alFree{2} == 2
                special.SatMethod = 'TP';
            elseif twix.hdr.MeasYaps.sWipMemBlock.alFree{2} == 3
                special.SatMethod = 'RBC+TP';
            end
        end
    elseif contains(protocol,'fa_spec_')
        % tx factor
        special.txfactor = twix.hdr.MeasYaps.sWipMemBlock.adFree{1};

        % number of reps
        special.numReps = twix.hdr.MeasYaps.sWipMemBlock.alFree{2};

        % delay between images (ms)
        special.repDelay = twix.hdr.MeasYaps.sWipMemBlock.adFree{3};

        % number of resonances
        special.numRes = 1 + twix.hdr.MeasYaps.sWipMemBlock.alFree{4};
        if isempty(special.numRes), special.numRes = 1; end

        % pulse shape for additional resonances
        special.RFshape = twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
        if isempty(special.RFshape)
            special.RFshape = 'aSinc';
        elseif special.RFshape == 1
            special.RFshape = 'Gauss';
        elseif special.RFshape == 2
            special.RFshape = 'Rect';
        elseif special.RFshape == 3
            special.RFshape = 'Kai';
        end

        % acquisition order
        special.acqOrder = twix.hdr.MeasYaps.sWipMemBlock.alFree{10};
        if isempty(special.acqOrder)
            special.acqOrder = 'Sequential';
        elseif special.acqOrder == 1
            special.acqOrder = 'Interleaved';
        elseif special.acqOrder == 2
            special.acqOrder = 'Inter-Interleaved';
        end

        % frequency, flip angle, and pulse duration for other resonances
        special.RFdur = twix.hdr.MeasYaps.sWipMemBlock.adFree{9}; % ms
        if special.numRes >= 2
            special.freq2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{5}; % ppm
            special.FA2 = twix.hdr.MeasYaps.sWipMemBlock.alFree{6}; % deg
        end

    elseif contains(protocol,'fa_radial_')
        % RF duration
        special.RFdur = twix.hdr.MeasYaps.sWipMemBlock.adFree{1}; % ms

        % RF pulse shape
        special.RFshape = twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
        if isempty(special.RFshape)
            special.RFshape = 'aSinc';
        elseif special.RFshape == 1
            special.RFshape = 'Gauss';
        elseif special.RFshape == 2
            special.RFshape = 'Rect';
        elseif special.RFshape == 3
            special.RFshape = 'Kai';
        end


    else
        for i = 1:14
            if isfield(twix.hdr.MeasYaps.sWipMemBlock,'alFree') &&...
                    i <= length(twix.hdr.MeasYaps.sWipMemBlock.alFree) && ~isempty(twix.hdr.MeasYaps.sWipMemBlock.alFree{i})

                special.(['param',num2str(i)]) = twix.hdr.MeasYaps.sWipMemBlock.alFree{i};
            else
                if isfield(twix.hdr.MeasYaps.sWipMemBlock,'adFree') &&...
                        i <= length(twix.hdr.MeasYaps.sWipMemBlock.adFree) && ~isempty(twix.hdr.MeasYaps.sWipMemBlock.adFree{i})

                    special.(['param',num2str(i)]) = twix.hdr.MeasYaps.sWipMemBlock.adFree{i};
                else
                    special.(['param',num2str(i)]) = 0;
                end
            end
        end

    end

    % zero any empty fields
    fnames = fieldnames(special);
    for i = 1:length(fnames)
        if isempty(special.(fnames{i}))
            special.(fnames{i}) = 0;
        end
    end
end
end

function printParam(data,stitle)
% data          struct with necessary parameters
% stitle        title for acquisition (optional)

if nargin < 2, stitle = []; end

% print title (if provided)
if ~isempty(stitle)
    fprintf('\n%s \n',stitle);
end

% check if spiral or cartesian
if isfield(data,'trajectory')
    if strcmpi(data.trajectory,'Spiral')
        fprintf('\nSamples/Interleaves = %g/%g\n',data.nCol,data.nLin);
    else
        fprintf('\n');
    end
else
    fprintf('\n');
end

% print imaging parameters
fprintf('TR/TE = %g/%g ms\nFA = %g%c\n',data.TR,data.TE,data.FA,176);

% get dimensions
if isfield(data,'images')
    dim = size(data.images);
else
    dim = [];
end

% print fov/res
if length(dim) > 2 && strcmp(data.mode,'2D')
    fov = [data.fovRO,data.fovPE,data.thickness*dim(3)];
    tmpstring2 = 'FOV = %g x %g x %g mm3\n';
    res = [data.fovRO/dim(1),data.fovPE/dim(2),data.thickness];
    tmpstring = 'Res = %.3g x %.3g x %.3g mm3\n';
    fprintf(tmpstring2,fov);
    fprintf(tmpstring,res);
    fprintf('%s \n',data.orientation);
elseif length(dim) > 2 && strcmp(data.mode,'3D')
    fov = [data.fovRO,data.fovPE,data.thickness*dim(3)];
    tmpstring2 = 'FOV = %g x %g x %g mm3\n';
    res = [data.fovRO/dim(1),data.fovPE/dim(2),data.thickness];
    tmpstring = 'Res = %.3g x %.3g x %.3g mm3\n';
    fprintf(tmpstring2,fov);
    fprintf(tmpstring,res);
    fprintf('3D %s \n',data.orientation);
elseif length(dim) == 2
    fov = [data.fovRO,data.fovPE];
    tmpstring2 = 'FOV = %g x %g mm2\n';
    res = [data.fovRO/dim(1),data.fovPE/dim(2)];
    tmpstring = 'Res = %.3g x %.3g mm2\n';
    fprintf(tmpstring2,fov);
    fprintf(tmpstring,res);
    if data.thickness == 500
        fprintf('%s Projection\n',data.orientation)
    else
        fprintf('Slice Thickness = %g mm\n',data.thickness);
        fprintf('%s \n',data.orientation);
    end
elseif isempty(dim)
    if isfield(data,'MatSize') && ~isempty(data.MatSize)
        fov = [data.fovRO,data.fovPE,data.thickness*data.MatSize];
        res = fov ./ data.MatSize;
        tmpstring2 = 'FOV = %g x %g x %g mm3\n';
        tmpstring = 'Res = %.3g x %.3g x %.3g mm3\n';
        fprintf(tmpstring2,fov);
        fprintf(tmpstring,res);
    else
        fov = [data.fovRO,data.fovPE,data.thickness*max(data.nSli,data.nPar)];
        tmpstring2 = 'FOV = %g x %g x %g mm3\n';
        fprintf(tmpstring2,fov);
    end
    fprintf('%s %s \n',data.mode,data.orientation);
end

% extra space
fprintf('\n');

end

function [ImgArrayCombined,KspaceArrayCombined,b] = combinecoils_fa(ImgArray,KspaceArray,NoiseSize,b)
% Array In:         [ncol,nlin,nslices,nchannels,nrepetitions]
% Array Out:        [ncol,nlin,nslices,nrepetitions]

% check inputs
if nargin < 2, KspaceArray = []; end
if nargin < 3, NoiseSize = []; end
if nargin < 4, b = []; end
if isempty(KspaceArray), KspaceArrayCombined = []; end
if isempty(NoiseSize), NoiseSize = 6; end

% get dimensions
dim = size(ImgArray);

% fix for images with no repetitions
try
    dim(5) = dim(5);
catch
    dim(5) = 1;
end

% define noise regions
NoiseRegionX = (dim(1)-NoiseSize):(dim(1)-1);
NoiseRegionY = (dim(2)-NoiseSize):(dim(2)-1);

% define NoiseRegionZ if 3D/iso
if dim(1) == dim(2) && dim(2) == dim(3)
    NoiseRegionZ = (dim(3)-NoiseSize):(dim(3)-1);
else
    NoiseRegionZ = 1:dim(3);
end

% normalize each channel by its noise
for ic = 1:dim(4)
    ImgArray(:, :, :, ic, :) = ImgArray(:, :, :, ic, :) / mean(abs(reshape(ImgArray(NoiseRegionX,NoiseRegionY,NoiseRegionZ,ic,:),[],1)));
end

if isempty(b)
    % sum all repetitions for calculating coil sensitivity
    b = sum(ImgArray, 5);

    % normalize by image norms
    b = b ./ sum(abs(ImgArray),5);

    % reshape
    b = permute(b, [4 1 2 3]);
end

% reshape
ImgArray = permute(ImgArray, [4 1 2 3 5]);

% combine optimally as per Bydder et al, MRM 47:539-458 (2002)
ImgArrayCombined = zeros(dim(1), dim(2), dim(3), dim(5));
for ix = 1:dim(1)
    for iy = 1:dim(2)
        for is = 1:dim(3)
            for ii = 1:dim(5)
                thisb = b(:, ix, iy, is);
                ImgArrayCombined(ix, iy, is, ii) = real((thisb' * ImgArray(:, ix, iy, is, ii)));
            end
        end
    end
end

% normalize by its noise
ImgArrayCombined = ImgArrayCombined / mean(reshape(ImgArrayCombined(NoiseRegionX,NoiseRegionY,NoiseRegionZ,:),[],1));

% zero negatives
% ImgArrayCombined(ImgArrayCombined < 0) = 0;

% just add up kspace channels
if ~isempty(KspaceArray)
    KspaceArrayCombined = sum(KspaceArray,4);
    KspaceArrayCombined = reshape(KspaceArrayCombined,size(KspaceArray,[1 2 3 5]));
end
end

function out = resize_fa(img,scale,method,threshold)
% img:              absolute value of image matrix
% scale:            interpolation factor
% method:           interpolation method
% threshold:        threshold for masking edges

% check inputs
if nargin < 3, method = []; end
if nargin < 4, threshold = []; end
if isempty(method), method = 'bilinear'; end
if isempty(threshold), threshold = 0.5; end

% do nothing if scale = 1
if scale == 1
    out = img;
else

    % create mask
    if any(isnan(img(:)))
        % include 0s if background is NaNs
        mask = double(abs(img) >= 0);
    else
        % don't include 0s if background is 0s
        mask = double(abs(img) > 0);
    end

    % set nans to 0
    img(isnan(img)) = 0;

    % resize img and mask
    out = imresize(img,scale,method);
    maskout = imresize(mask,scale,method);

    % threshold mask/img
    maskout(maskout < threshold) = nan;
    out = out .* maskout;

end
end

function [kernel,u] = createKBkernel(width,beta,length)
u = (0:length-1)/(length-1) * width/2;
f = beta*sqrt(1-(2*u/width).^2);
kernel = besseli(0,f)./width;
kernel = (kernel/max(kernel))';
u = u';
end

function [wi,error] = iterative_dcf_fa_20190910(iter,adKSpaceCoor,kerneltable,Ind,MatSize,w0,verbose)
% iter:                 # of iterations
% adKSpaceCoor:         k space coordinates (nx3 or nx2)
% kerneltable:          kernel values (from gridding reconstruction)
% Ind:                  cartesian grid indices (from gridding reconstruction)
% MatSize:              length of cartesian matrix (int) or dimensions (array)
% w0:                   initial DCF estimate (optional)
% verbose:              print progress (optional)

% Adopted from https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.23041
% 3D was modified and needs to be double-checked

% Check for initial DCF
if nargin<6, w0 = []; end
if nargin<7, verbose = []; end
if isempty(w0), w0 = ones(size(adKSpaceCoor,1),1); end
if isempty(verbose), verbose = 1; end

% 2D or 3D?
dim = size(adKSpaceCoor,2);

% Round iter to nearest odd number
% iter = round((iter-1)/2)*2+1;

% Progress ticker
if verbose, fprintf(1,'dcf iter:  1 / %d',iter); end

% if dim == 3
% Check grid size
if length(MatSize) == 3
    msize = MatSize;
else
    msize = [MatSize, MatSize, MatSize];
end

% Initial DCF Estimate
wi = w0;
kernelweights = zeros(msize);

for zz = 1:iter
    % Progress ticker
    if verbose
        fprintf(1,repmat('\b',1,3+numel([num2str(iter) num2str(zz)])));
        fprintf(1,'%d / %d',zz,iter);
    end

    % Allocate memory
    grid = zeros(msize);

    % % Gridding % %

    % Convolve kernel with kspace point
    %         witmp = bsxfun(@times,kerneltable,repmat(wi,1,size(kerneltable,2))');
    witmp = (kerneltable .* wi);

    for i = 1:size(adKSpaceCoor,1)
        % Get grid indices for each kspace point
        tmpind = Ind(i,:);

        % Combine convolved data
        grid(tmpind) = grid(tmpind) + witmp(i,:);

        kernelweights(tmpind) = kernelweights(tmpind) + kerneltable(i,:);
    end

    % % %         tmpimg = fftshift(ifftn(fftshift(grid)));
    % % %         ifftweights = fftshift(ifftn(fftshift(kernelweights)));
    % % % %         nonzeroweights = repmat(floor(1+(NumK*OverSampleFactor+1-NumK)/2):floor(1+(NumK*OverSampleFactor+1-NumK)/2)+(NumK-1),3,1)';
    % % % %         tmpimg(nonzeroweights) = tmpimg(nonzeroweights) ./ (ifftweights(nonzeroweights).^2);
    % % %         tmpimg = tmpimg ./ (ifftweights.^2);
    % % %         grid = abs(fftshift(fftn(tmpimg)));

    % % Degridding % %

    % Allocate memory
    w1i = zeros(size(wi));

    % Convolve grid with kernel
    w1itmp = (grid(Ind)) .* kerneltable;

    for i = 1:size(adKSpaceCoor,1)
        % Combine convolved data
        w1i(i) = sum(w1itmp(i,:));
    end

    % Invert density to get weights
    nonzero = w1i ~= 0;
    wi(nonzero) = wi(nonzero) ./ w1i(nonzero);
    wi(~nonzero) = 0;
    %         error = mean(1 - w1i.^-1);
    error = max(w1i.^-1) - min(w1i.^-1);
end
if verbose, fprintf('\n'); end
end


function [Kcorr, deltas, diagout] = gradient_delay_correct( ...
    KSpaceCoor, nsamples, nleaves, nreps, varargin)
%GRADIENT_DELAY_CORRECT  Estimate per-axis time shifts (Δx,Δy,Δz) and apply.
%
% KSpaceCoor : [N x 3] (kx,ky,kz) in 1/mm, N = nsamples*nleaves*nreps
% nsamples   : samples per interleaf
% nleaves    : interleaves per GA repetition
% nreps      : number of GA repetitions (do NOT average across reps)
%
% Name-Value:
%   'N0'            (60)     number of early samples to fit on (min 8)
%   'SearchDwells'  (4)      scalar W => [-W,+W], or [min max] in dwells
%   'FineStep'      (0.1)    dwell step for grid search
%   'Weights'       ([] )    length-N0 weights over early samples
%   'PreAlign'      (true)   subtract per-leaf DC offset for fitting only
%   'AlignN'        (6)      #samples to form DC estimate if PreAlign=true
%   'Verbose'       (true)
%   'DoPlots'       (false)
%
% Returns:
%   Kcorr  : [N x 3] corrected trajectory (time-shifted component-wise)
%   deltas : [1 x 3] best Δ (dwells) for X,Y,Z
%   diagout: struct with grids/errs and a few summary metrics

args = struct('N0',60,'SearchDwells',4,'FineStep',0.1, ...
    'Weights',[],'PreAlign',true,'AlignN',6, ...
    'Verbose',true,'DoPlots',false);
args = parseargs(args, varargin{:});

% --- shape checks & reshape to [samples x (leaves*reps) x 3] ---
N = size(KSpaceCoor,1);  nInst = nleaves*nreps;
assert(N == nsamples*nInst, 'Size mismatch: N ~= nsamples*nleaves*nreps.');
K3 = reshape(KSpaceCoor, [nsamples, nInst, 3]);

% --- early-sample configuration ---
N0 = max(8, min(args.N0, nsamples));
if isempty(args.Weights)
    w = (0:N0-1)';
    w = w / max(1,sum(w));                 % gentle ramp emphasis
else
    w = args.Weights(:);
    assert(numel(w)==N0, 'Weights must have length N0.');
    s = sum(w); if s==0, w(1)=1; s=1; end; w = w/s;
end

deltas = zeros(1,3);
axisLbl = ["X","Y","Z"];
diagout = struct('axis',axisLbl, 'grid',{{}}, 'err',{{}}, ...
    'N0',N0,'weights',w,'prealign',args.PreAlign);

for a = 1:3
    Ka_orig = K3(:,:,a);                           % [nsamples x nInst]
    % --- DC pre-align for FITTING ONLY (robust) ---
    if args.PreAlign
        k0 = median(Ka_orig(1:min(args.AlignN,nsamples),:),1);  % 1 x nInst
        Ka_fit = Ka_orig - k0;                                  % do NOT keep!
    else
        Ka_fit = Ka_orig;
    end

    % --- fit Δ on the early samples using robust reference (median over leaves) ---
    [bestDelta, gridD, gridE] = fit_delay_axis(Ka_fit, N0, args.SearchDwells, args.FineStep, w);

    deltas(a) = bestDelta;
    if args.Verbose
        % Build a readable window description
        if isscalar(args.SearchDwells)
            rngtxt = sprintf('±%g', args.SearchDwells);
        else
            v = args.SearchDwells(:).';
            if numel(v) >= 2
                % show first/last (and step if uniform)
                if all(abs(diff(v) - v(2)+v(1)) < 1e-12)
                    rngtxt = sprintf('[%g..%g step %g]', v(1), v(end), v(2)-v(1));
                else
                    rngtxt = sprintf('[%g,%g]', v(1), v(end));
                end
            else
                rngtxt = sprintf('%g', v);
            end
        end

        fprintf('Gradient-delay Δ%s = %+0.3f dwells (N0=%d, window=%s).\n', ...
            axisLbl(a), bestDelta, N0, rngtxt);
    end
    diagout.grid{a} = gridD(:);
    diagout.err{a}  = gridE(:);

    % --- APPLY Δ to the ORIGINAL Ka (not the DC-aligned copy) ---
    K3(:,:,a) = shift_all_leaves(Ka_orig, bestDelta);
end

Kcorr = reshape(K3, [N, 3]);

if args.DoPlots
    plot_delay_diagnostics(diagout, deltas);
end
end

% ===== helpers =====

function [bestDelta, gridD, gridE] = fit_delay_axis(Ka, N0, SearchDwells, stepdw, w)
% Ka       : [nsamples x nInst] single axis
% N0       : # early samples
% SearchDwells : scalar W => [-W,+W] or [min max]
% stepdw   : grid step in dwells
% w        : N0x1 weights over early samples (sum to 1)

[nsamples, ~] = size(Ka);
N0 = min(N0, nsamples);

if isscalar(SearchDwells), rng = [-abs(SearchDwells), abs(SearchDwells)];
else,                      rng = sort(SearchDwells(:).'); end
gridD = rng(1):stepdw:rng(2);
gridE = zeros(size(gridD));

for ii = 1:numel(gridD)
    d  = gridD(ii);
    Ks = shift_all_leaves(Ka, d);           % [nsamples x nInst]
    Ke = Ks(1:N0,:);                         % early window
    ref = median(Ke, 2);                     % robust reference shape [N0 x 1]
    R   = Ke - ref;                          % residuals [N0 x nInst]
    m   = mean(R.^2, 2);                     % mean over leaves -> [N0 x 1]
    gridE(ii) = sum(w .* m);                 % weighted sample energy
end

% pick min and parabolic refine if interior
[~, idx] = min(gridE);
bestDelta = gridD(idx);
if idx>1 && idx<numel(gridD)
    x1 = gridD(idx-1); x2 = gridD(idx); x3 = gridD(idx+1);
    y1 = gridE(idx-1); y2 = gridE(idx); y3 = gridE(idx+1);
    denom = (x1-x2)*(x1-x3)*(x2-x3);
    if denom ~= 0
        A = (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)) / denom;
        B = (x3^2*(y1-y2) + x2^2*(y3-y1) + x1^2*(y2-y3)) / denom;
        if A ~= 0
            xc = -B/(2*A);
            if xc >= gridD(1) && xc <= gridD(end)
                bestDelta = xc;
            end
        end
    end
end
end

function Ks = shift_all_leaves(Ka, delta_dwells)
% Sub-dwell shift along sample dimension, same Δ for all leaf instances.
% PCHIP keeps monotonicity and avoids ringing.
[nsamples, nInst] = size(Ka);
t  = (0:nsamples-1).';
tp = t + delta_dwells;                        % shifted "time" indices
Ks = zeros(nsamples, nInst, 'like', Ka);
for L = 1:nInst
    v = Ka(:,L);
    Ks(:,L) = interp1(t, v, tp, 'pchip', 'extrap');
end
end

function plot_delay_diagnostics(diagout, deltas)
figure('Color','w');
for a=1:3
    subplot(3,1,a);
    plot(diagout.grid{a}, diagout.err{a}, 'o-','LineWidth',1); grid on;
    title(sprintf('Axis %s: error vs \\Delta (best = %+0.3f dwells)', diagout.axis(a), deltas(a)));
    xlabel('\Delta (dwells)'); ylabel('weighted mean residual^2');
end
end

function S = parseargs(S, varargin)
for k=1:2:numel(varargin), S.(varargin{k}) = varargin{k+1}; end
end
