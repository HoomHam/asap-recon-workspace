close all
clc
clear

% 
img = phantom(256,256);
figure;
imshow(img, [])

% trajectory
points = 500;
shots = 25;
turns = 4;
cycle = points/turns;

dkx = zeros(points, shots);
dky = zeros(points, shots);

for indx=1:points
    dkx(indx, 1) = indx/points*sin(2*pi*indx/cycle) * 0.5; %
    dky(indx, 1) = indx/points*cos(2*pi*indx/cycle) * 0.5;
end

ga = deg2rad(360/shots);
for indx=2:shots
    dkx(:, indx) = cos(ga) * dkx(:, indx-1) - sin(ga) * dky(:, indx - 1);
    dky(:, indx) = sin(ga) * dkx(:, indx-1) + cos(ga) * dky(:, indx - 1);
end

figure;
plot(dkx(:,1:10:end), dky(:,1:10:end));
axis equal
grid on

dk = dkx + 1j*dky;
k = reshape(dk, [points*shots, 1]);
k = k(:)/max(abs(k(:)))/2;

% 
% generate circular mask (spirals hav a circular FOV support
N = 256;
% [xx,yy] = meshgrid(linspace(-1,1,N));
% ph = double(sqrt(xx.^2 + yy.^2)<1);

ph = 1;
FT = NUFFT(k, 1, ph, 0, [N,N], 2);

data = FT*img;

% 
w = voronoidens(k); % calculate voronoi density compensation function
w = w/max(w(:));

im_dc = FT'*(data.*w);

figure;
imshow(abs(im_dc), []);

% 
TVWeight = 0.01; 	% Weight for TV penalty
xfmWeight = 0.00;	% Weight for Transform L1 penalty
Itnlim = 25;		% Number of iterations

XFM = 1;				% Identity transform 	

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

% im_dc = FT'*(data.*w);	% init with zf-w/dc (zero-fill with density compensation)
% figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
tic
for n=1:8
	res = fnlCg(res,param);
	im_res = XFM'*res;
	figure(100+n), imshow(abs(im_res),[]), drawnow
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
end
toc

