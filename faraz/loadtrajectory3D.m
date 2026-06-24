function [KSpaceCoor,imgsize] = loadtrajectory3D(frequency,fov,nsamples,nleaves,nreps,imgsize,orientation,nucleus,tolerance,calibfile)
% INPUTS 
% frequency:        center frequency (MHz/T)
% fov:              field of view (mm)
% nsamples:         number of points per interleave
% nleaves:          number of interleaves
% nreps:            number of repetitions/golden-angle rotations
% orientation:      'axial', 'coronal',or 'sagittal'
% nucleus           'water','gas,'dissolved'
% tolerance:        tolerance for matching frequency and fov (optional)
% calibfile:        '*mat' file containing calibrations (optional)

% OUTPUTS
% KSpaceCoor:       corrected trajectory [nsamples*nleaves*nreps,3]
% imgsize:          matrix size from matched calibration entry

% check inputs
if nargin < 10, calibfile = [];  end
if nargin < 9, tolerance = [];  end
if nargin < 8, nucleus = [];  end
if nargin < 7, orientation = [];end
if nargin < 6, imgsize = [];    end
if nargin < 5, nreps = [];      end
if nargin < 4, nleaves = [];    end
if nargin < 3, nsamples = [];   end
if nargin < 2, fov = [];        end
if nargin < 1, frequency = [];  end

% assign defaults
if isempty(calibfile), calibfile = 'calibrations_3D_20220308.mat'; end
if isempty(tolerance), tolerance = 0.05; end
if isempty(imgsize), imgsize = 80; end
if isempty(nucleus)
    if frequency > 30
        nucleus = 'water';
    else
        nucleus = 'gas';
    end
end

% find calibration file 
filename = which(calibfile);

% load calibrations
s = load(filename);
s = s.s;

% backward compatibility
if ~isfield(s,'nucleus')
    for i = 1:length(s)
        if s(i).frequency > 30
            s(i).nucleus = 'water';
        else
            s(i).nucleus = 'gas';
        end
    end
end

% search for calibrated trajectory
match = arrayfun(@(x) ...
    abs(x.frequency - frequency) <= tolerance*frequency && ...
    abs(x.fov - fov)             <= tolerance*fov       && ...
    abs(x.nsamples - nsamples)   <= 0                   && ...
    abs(x.nleaves - nleaves)     <= 0                   && ...
    abs(x.nreps - nreps)         <= 0                   && ...
    abs(x.imgsize - imgsize)     <= 0                   && ...
    strcmpi(x.nucleus,nucleus)                          && ...
    strcmpi(x.orientation,orientation), s);

% get index for last (most recent) match
ind = find(match == 1,1,'last');

% get kspace coordinates and diameter if available
if isempty(ind)
    disp('No 3D trajectory calibration available');
    KSpaceCoor = [];
    imgsize = [];
else
    KSpaceCoor = [s(ind).kx, s(ind).ky, s(ind).kz];
    if isempty(s(ind).imgsize)
        imgsize = [];
    else
        imgsize = s(ind).imgsize;
    end
    disp(['3D Trajectory calibration from ',s(ind).date,' loaded']);
end

end