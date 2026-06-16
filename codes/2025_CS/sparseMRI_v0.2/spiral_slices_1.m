% note : to run this code, first run setup.m in irt folder to setup paths
% for nufft toolbox; then select mapVBVD-main folder and right click to add
% path with subfolders
% voronoidens : add sparseMRI_v0.2/utils path

% the script tests undersampling for spiral;
close all
clc
clear

% 
% measID = 'C:\Users\P53-LOCAL\Desktop\AI\compressed_sensing\python\meas_MID00670_FID05512_flashback_spiral_multishot_2000p_80shots_shimming';
measID = 'C:\Users\P53-LOCAL\Desktop\AI\compressed_sensing\python\meas_MID00544_FID16608_flashback_spiral_fov250';
twix = mapVBVD(measID);

% return all image-data
% rawdata = twix{2}.image();

% return all image-data with all singular dimensions removed/squeezed:
rawdata = twix{2}.image{''}; % '' necessary due to a matlab limitation

[readouts, channels, lines, slices] = size(rawdata);

turns = 4;
N = lines*turns*2;

% trajectory
points = readouts;
shots = fix(lines/2);
cycle = points/turns;

dkx = zeros(points, shots);
dky = zeros(points, shots);

for indx=1:points
    dkx(indx, 1) = indx/points*sin(2*pi*indx/cycle) * pi; %125 : 4turns, 15.625 : 32turns
    dky(indx, 1) = indx/points*cos(2*pi*indx/cycle) * pi;
end

ga = deg2rad(360/shots);
for indx=2:shots
    dkx(:, indx) = cos(ga) * dkx(:, indx-1) - sin(ga) * dky(:, indx - 1);
    dky(:, indx) = sin(ga) * dkx(:, indx-1) + cos(ga) * dky(:, indx - 1);
end

figure(1);
plot(dkx(:,:), dky(:,:));
axis equal
grid on

dk = dkx + 1j*dky;
k = reshape(dk, [points*shots, 1]);
k = k(:)/max(abs(k(:)))/2;

% generate circular mask (spirals hav a circular FOV support
% [xx,yy] = meshgrid(linspace(-1,1,N));
% ph = double(sqrt(xx.^2 + yy.^2)<1);

ph = 1;
data = squeeze(rawdata(:,1,1:2:end,2));
FT = NUFFT(k, 1, ph, 0, [N, N], 2);

% w = 1;
w = voronoidens(k); % calculate voronoi density compensation function
w = w/max(w(:));
data = reshape(data, [points*shots, 1]);
im_dc = FT'*(data.*w);

figure;
imshow(abs(fliplr(im_dc)'), []);
