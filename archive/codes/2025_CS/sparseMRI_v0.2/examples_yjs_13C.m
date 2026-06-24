close all
clc
clear

%% spiral 
points = 500;
indx = 1:1:points;
turns = 4;
shots = 8;
cycle = points/turns;

FOV = 75; % mm
kmax = turns*shots/(FOV/10);  % 1/cm 
% kmax = turns*shots/(FOV/10) *2;  % for CS, *2 

kx = indx/points.*sin(2*pi*indx/cycle) * kmax; %
ky = indx/points.*cos(2*pi*indx/cycle) * kmax;

figure('Name','Orginal KxKy');
plot(kx, ky);
axis equal
grid on

kz = 0*ky';
C = [kx' ky' kz];

% gamma = 4.257; % 1H
gamma = 1.071; % MHz/T (13C)

gmax = 2.0; % G/cm
smax = 5.0; % G/cm/ms

[C_riv, time_riv, g_riv, s_riv, k_riv] = minTimeGradient(C,0, 0, 0, gmax, smax, gamma, 10e-3);          % Rotationally invariant solution

save('C:\Users\P53-LOCAL\Desktop\Trio\rawdata\coil_test_04162025\ktraj_Arch.mat', 'k_riv');

N = length(g_riv);
t = linspace(0,time_riv, N);

% 
figure('Name', 'gradient (time)');
plot(t, g_riv*10);
xlabel('time (ms)');
ylabel('gradient (mT/m)');
legend('gx', 'gy', 'gz');
title('gradient');
grid on

% 
figure('Name', 'gradient (gx, gy)');
plot(g_riv(:,1)*10, g_riv(:,2)*10);
xlabel('gx (mT/m)');
ylabel('gy (mT/m)');
title('gradient');
grid on
axis equal

% 
N = length(s_riv);
t = linspace(0,time_riv, N);
figure('Name', 'slewrate (time)');
plot(t, s_riv*10);
xlabel('time (ms)');
ylabel('slewrate (mT/m/ms)');
legend('sx', 'sy', 'sz');
title('slew rate');
grid on

% 
figure('Name', 'slewrate (sx, sy)');
plot(s_riv(:,1)*10, s_riv(:,2)*10);
xlabel('sx (mT/m/ms)');
ylabel('sy (mT/m/ms)');
title('slew rate');
grid on
axis equal

% 
N = length(C);
t1 = linspace(0,time_riv, N);
figure('Name', 'k-trajectory comparison (time)');
% plot(t1, C(:,1:2)*100);
% hold on

N = length(C_riv);
t2 = linspace(0,time_riv, N);
plot(t2, C_riv(:,1:2)*100);
xlabel('time (ms)');
ylabel('k_trajectory (1/m)');
% legend('Cx', 'Cy', 'kx', 'ky');
legend('kx', 'ky');
title('k-space trajectory');
grid on

% 
figure('Name', 'K-tractory');
plot(C(:,1)*100, C(:,2)*100);
hold on
plot(C_riv(:,1), C_riv(:,2));
legend('kxy', 'Cxy');
title('k-space trajectory');
grid on
axis equal

