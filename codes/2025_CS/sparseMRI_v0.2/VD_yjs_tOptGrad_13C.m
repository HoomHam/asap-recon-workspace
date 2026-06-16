% add nufft toolbox path : C:\Users\P53-LOCAL\Desktop\AI\compressed_sensing\codes\irt

close all
clear
clc
% addpath(strcat(pwd,'/utils'));

% generate a variable density spiral with gaussian density.
SIGMA = 3;
FOV = gausswin(40,SIGMA);
FOV = FOV(end/2+1:end);
FOV = FOV*7.4+0.1; %cm

% FOV = ones(20,1) * 7.5;

figure;
plot(FOV);
grid on
ylabel("FOV (cm)");

RADIUS = linspace(0,1,20);
RESOLUTION = 75/64; % in mm
NITLV = 4; % number of spiral interleaves: 16

% gamma = 4.257; % 1H
gamma = 1.071; % MHz/T (13C)

Gmax = 2.0 ; % maximum gradient in [G/CM];
Smax = 5; %Maximum slew-rate; [G/(CM*ms)]
T = 10e-3; % time sampling (in mS);

disp('design spiral')
[k,g,s,time] = vdSpiralDesign(NITLV, 0, RESOLUTION,FOV,RADIUS,gamma, Gmax,Smax,T,[],'cubic'); % TODO
% k = k(2:end).'*exp(2*pi*i*[1:NITLV]/NITLV);
% k = k(:)/max(abs(k(:)))/2; % scale to range [-0.5,0.5]

% plot k-trajectory in (kx. ky);
readout = length(k);
figure('Name', "k-trajectory (kx, ky)");
plot(k(:,1), k(:,2));
axis equal
grid on
title('K-trajectory (kx vs ky)');

% plot k-trajectory in time domain;
t = linspace(0,time,readout);
figure('Name', "K-trajectory: time domain");
plot(t, k(:,1));
hold on
plot(t, k(:,2));
grid
xlabel('time (ms)');
ylabel('k traj (1/cm)');

% plot gradients in (gx, gy);
g = g*10; % G/cm -> mT/m;
figure('Name', "g-trajectory (gx, gy)");
plot(g(:,1), g(:,2));
axis equal
grid on
title('gradient (gx vs gy)');

% plot gradients in time domain;
figure('Name', "g-trajectory : time domain");
plot(t, g(:,1));
hold on
plot(t, g(:,2));
hold on
grid
xlabel('Time (ms)');
ylabel('Gradient (mT/m)');

% plot slewrate in (sx, sy);
s = s*10; % mT/(m*ms)
figure('Name', "s-trajectory (sx, sy)");
plot(s(:,1), s(:,2));
axis equal
grid

% plot slewrate in time domain;
figure('Name', "s-trajectory : time domain");
t = linspace(0,time,length(s));
plot(t, s(:,1));
hold on
plot(t, s(:,2));
grid
xlabel('Time (ms)');
ylabel('Slew rate (mT/(m*ms))');

% save('C:\Users\P53-LOCAL\Desktop\Trio\rawdata\coil_test_04162025\ktraj_vd_fov75_8shots.mat', 'k');
