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

dkx = reshape(dkx, [points*shots, 1]);
dky = reshape(dky, [points*shots, 1]);
dk = dkx+1j*dky;

% [-0.5, 0.5]
dkx = dkx(:)/max(abs(dk(:)))/2;
dky = dky(:)/max(abs(dk(:)))/2;

%%
om = [dkx, dky]*2*pi;
Nd = [N, N];
Jd = [6,6];
Kd = Nd*2;
n_shift = Nd/2;
st = nufft_init(om, Nd, Jd, Kd, n_shift);

data = squeeze(rawdata(:,1,:,1));
data = reshape(data, [points*shots, 1]);

w = voronoidens(dkx+1j*dky); 
w = w/max(w(:));

image = nufft_adj(data.*w, st);

figure;
imshow(abs(image), []);

% 
% dk = dkx + 1j*dky;
% k = reshape(dk, [points*shots, 1]);
% k = k(:)/max(abs(k(:)))/2;
% 
% % generate circular mask (spirals hav a circular FOV support
% % [xx,yy] = meshgrid(linspace(-1,1,N));
% % ph = double(sqrt(xx.^2 + yy.^2)<1);
% ph = 1;
% 
% data = squeeze(rawdata(:,1,:,1));
% FT = NUFFT(k, 1, ph, 0, [N, N], 2);
% 
% % w = 1;
% w = voronoidens(k); % calculate voronoi density compensation function
% w = w/max(w(:));
% data = reshape(data, [points*shots, 1]);
% im_dc = FT'*(data.*w);
% 
% figure;
% imshow(abs(fliplr(im_dc)'), []);
