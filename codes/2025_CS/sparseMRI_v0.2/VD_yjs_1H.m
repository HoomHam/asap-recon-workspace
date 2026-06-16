% add nufft toolbox path : C:\Users\P53-LOCAL\Desktop\AI\compressed_sensing\codes\irt

close all
clear
clc
addpath(strcat(pwd,'/utils'));

% generate a variable density spiral with gaussian density.
SIGMA = 3;
FOV = gausswin(40,SIGMA);
FOV = FOV(end/2+1:end);
FOV = FOV*20+4; %cm

% FOV = ones(20,1) * 24;

figure;
plot(FOV);
grid on
ylabel("FOV (cm)");

RADIUS = linspace(0,1,20);
RESOLUTION = 0.5; % in mm
NITLV = 16; % number of spiral interleaves: 16

gamma = 4.257; % 1H

Gmax = 1.8 ; % maximum gradient in [G/CM];
Smax = 5; %Maximum slew-rate; [G/(CM*ms)]
T = 10e-3; % time sampling (in mS);

disp('design spiral')
[k,g,s,time] = vdSpiralDesign_old(NITLV, RESOLUTION,FOV,RADIUS,gamma, Gmax,Smax,T,'cubic');
% k = k(2:end).'*exp(2*pi*i*[1:NITLV]/NITLV);
% k = k(:)/max(abs(k(:)))/2; % scale to range [-0.5,0.5]

% plot k-trajectory in (kx. ky);
readout = length(k);
figure('Name', "k-trajectory (kx, ky)");
plot(real(k), imag(k));
axis equal
grid on
title('K-trajectory (kx vs ky)');

% plot k-trajectory in time domain;
t = linspace(0,time,readout);
figure('Name', "K-trajectory: time domain");
plot(t, real(k));
hold on
plot(t, imag(k));
grid
xlabel('time (ms)');
ylabel('k traj (1/cm)');

% plot gradients in (gx, gy);
g = g*10; % G/cm -> mT/m;
figure('Name', "g-trajectory (gx, gy)");
plot(real(g), imag(g));
axis equal
grid on
title('gradient (gx vs gy)');

% plot gradients in time domain;
figure('Name', "g-trajectory : time domain");
plot(t, real(g));
hold on
plot(t, imag(g));
hold on
grid
xlabel('Time (ms)');
ylabel('Gradient (mT/m)');

% plot slewrate in (sx, sy);
s = s*10; % mT/(m*ms)
figure('Name', "s-trajectory (sx, sy)");
plot(real(s), imag(s));
axis equal
grid

% plot slewrate in time domain;
figure('Name', "s-trajectory : time domain");
plot(t, real(s));
hold on
plot(t, imag(s));
grid
xlabel('Time (ms)');
ylabel('Slew rate (mT/(m*ms))');

save('C:\Users\P53-LOCAL\Desktop\PAV10\rawdata\time_optimal\ktraj_vd_fov240_16shots.mat', 'k');
