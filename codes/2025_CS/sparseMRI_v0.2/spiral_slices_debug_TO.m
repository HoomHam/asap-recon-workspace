close all
clc
clear

% addpath('mapVBVD-main/');

measID = 'C:\Users\P53-LOCAL\Desktop\AI\compressed_sensing\python\meas_MID00670_FID05512_flashback_spiral_multishot_2000p_80shots_shimming';
% measID = 'C:\Users\P53-LOCAL\Desktop\AI\compressed_sensing\python\meas_MID00544_FID16608_flashback_spiral_fov250';
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
shots = lines;
% points = 500;
% shots = 4;
cycle = readouts/turns;

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
plot(dkx(:,1:10:end), dky(:,1:10:end));
axis equal
grid on

dk = dkx + 1j*dky;
k = reshape(dk, [points*shots, 1]);
k = k(:)/max(abs(k(:)))/2;

% generate circular mask (spirals hav a circular FOV support
% [xx,yy] = meshgrid(linspace(-1,1,N));
% ph = double(sqrt(xx.^2 + yy.^2)<1);
ph = 1;

data = squeeze(rawdata(:,1,:,2));
FT = NUFFT(k, 1, ph, 0, [N, N], 2);

% w = 1;
w = voronoidens(k); % calculate voronoi density compensation function
w = w/max(w(:));
data = reshape(data, [points*shots, 1]);
im_dc = FT'*(data.*w);

figure;
imshow(abs(fliplr(im_dc)'), []);

% % % time-optimal k-trajectory
file = load('C:\Users\P53-LOCAL\Desktop\PAV10\rawdata\time_optimal\k_traj.mat');
k_traj = file.k_riv;

kx = k_traj(2:end,1);
ky = k_traj(2:end,2);

turns = 4;
shots = 80;

figure;
plot(k_traj(:,1:2));

N = lines*turns*2;

% trajectory
points = length(k_traj)-1;

dkx = zeros(points, shots);
dky = zeros(points, shots);

dkx(:, 1) = kx;
dky(:, 1) = ky;

ga = deg2rad(360/shots);
for indx=2:shots
    dkx(:, indx) = cos(ga) * dkx(:, indx-1) - sin(ga) * dky(:, indx - 1);
    dky(:, indx) = sin(ga) * dkx(:, indx-1) + cos(ga) * dky(:, indx - 1);
end

figure(1);
plot(dkx(:,1:20:end), dky(:,1:20:end));
axis equal
grid on

dk = dkx + 1j*dky;
k = reshape(dk, [points*shots, 1]);
k = k(:)/max(abs(k(:)))/2;

ph = 1;
FT = NUFFT(k, 1, ph, 0, [N, N], 2);

% % 
k_spiral = FT*im_dc;

w = voronoidens(k); % calculate voronoi density compensation function
w = w/max(w(:));

im_dc1 = FT'*(k_spiral.*w);
% im_dc1 = FT'*(k_spiral);

figure('Name','time-optimal trajectory');
imshow(abs(fliplr(im_dc1)'), []);

k_spiral2d = reshape(k_spiral, [points, shots]);
figure;
plot(abs(k_spiral2d(:,1))/max(abs(k_spiral2d(:,1))));
hold on

data2d = reshape(data, [2000,shots]);
plot(abs(data2d(1:points,1))/max(abs(data2d(1:points,1))));
grid on
title('time-optimal');

% % 
measID = 'C:\Users\P53-LOCAL\Desktop\PAV10\rawdata\time_optimal\meas_MID00640_FID18248_spiral_timeOptimal';
twix = mapVBVD(measID);

% return all image-data with all singular dimensions removed/squeezed:
rawdata = twix{2}.image{''}; % '' necessary due to a matlab limitation

[readouts, channels, lines, slices] = size(rawdata);

data_slice2 = rawdata(:,1,:,2);
figure;
% plot(abs(data_slice2(:,1:5))/max(abs(data_slice2(:,1:5))));
plot(abs(data_slice2(:,1:1:end)));
grid on
title('TO scanner data');

figure;
% plot(abs(k_spiral2d(:,1:5))/max(abs(k_spiral2d(:,1:5))));
plot(abs(k_spiral2d(:,1:1:end)));
grid on
title('simulation');
% legend('scanner', 'simulation');

figure;
% plot(abs(k_spiral2d(:,1:5))/max(abs(k_spiral2d(:,1:5))));
plot(abs(data2d(:,1:1:end)));
grid on
title('scanner CAV');
