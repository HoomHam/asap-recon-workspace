function savetrajectory3D(frequency,fov,nsamples,nleaves,nreps,orientation,thickness,imgsize,nucleus,kx,ky,kz,date,calibfile)
% INPUTS 
% frequency:        center frequency (MHz/T)
% fov:              field of view (mm)
% nsamples:         number of points per interleave
% nleaves:          number of interleaves
% nreps:            number of repetitions/golden-angle rotations
% orientation:      'axial', 'coronal',or 'sagittal'
% thickness:        slice thickness used for calibration (mm)
% nucleus:          'gas' or 'dissolved' trajectory for xenon frequencies
% imgsize:          imgsize for reconstruction
% kx:               x gradient data [nsamples*nleaves*nreps,1]
% ky:               y gradient data [nsamples*nleaves*nreps,1]
% kz:               z gradient data [nsamples*nleaves*nreps,1]
% data:             date of trajectory calculation (optional)
% calibfile:        '*mat' file containing calibrations (optional)


% check inputs
if nargin < 14, calibfile = [];     end
if nargin < 13, date = [];          end
if nargin < 12, kz = [];            end
if nargin < 11, ky = [];            end
if nargin < 10, kx = [];            end
if nargin < 9,  nucleus = [];       end
if nargin < 8,  imgsize = [];       end
if nargin < 7,  thickness = [];     end
if nargin < 6,  orientation = [];   end
if nargin < 5,  nreps = [];         end
if nargin < 4,  nleaves = [];       end
if nargin < 3,  nsamples = [];      end
if nargin < 2,  fov = [];           end
if nargin < 1,  frequency = [];     end

% default date is current date
if isempty(date), date = datestr(now,'yyyy/mm/dd'); end

% find calibration file 
if isempty(calibfile), calibfile = 'calibrations3D.mat'; end
filename = which(calibfile);

% check if filename exists
if isempty(filename)
   filename = calibfile; 
end

% load calibrations
try
    s = load(filename);
    s = s.s;

    % get last index
    ind = 1+size(s,2);

catch
    s = struct;
    ind = 1;
end

% add new calibration
s(ind).frequency = frequency;
s(ind).fov = fov;
s(ind).nsamples = nsamples;
s(ind).nleaves = nleaves;
s(ind).nreps = nreps;
s(ind).orientation = orientation;
s(ind).thickness = thickness;
s(ind).imgsize = imgsize;
s(ind).nucleus = nucleus;
s(ind).kx = kx;
s(ind).ky = ky;
s(ind).kz = kz;
s(ind).date = date;

% save new .mat file
save(filename,'s');
disp(['Trajectory saved to ',filename]);

end