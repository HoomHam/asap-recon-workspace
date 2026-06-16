close all
clc
clear

%% spiral 
points = 2000;
indx = 1:1:points;
turns = 4;
shots = 80;
cycle = points/turns;

FOV = 400; % mm
kmax = turns*shots/(FOV/10);  % 1/cm

kx = indx/points.*sin(2*pi*indx/cycle) * kmax; %
ky = indx/points.*cos(2*pi*indx/cycle) * kmax;

figure('Name','Orginal KxKy');
plot(kx, ky);
axis equal
grid on

kz = 0*ky';
C = [kx' ky' kz];

gamma = 4.257; % 1H
% gamma = 1.071; % MHz/T (13C)
gmax = 2.2;
smax = 10;

[C_riv, time_riv, g_riv, s_riv, k_riv] = minTimeGradient(C,0, 0, 0, gmax, smax, gamma, 10e-3);          % Rotationally invariant solution

save('C:\Users\P53-LOCAL\Desktop\PAV10\rawdata\time_optimal\ktraj_Arch.mat', 'k_riv');

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
plot(t2, k_riv(:,1:2)*100);
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
plot(k_riv(:,1), k_riv(:,2));
legend('Cxy', 'kxy');
title('k-space trajectory');
grid on
axis equal

return;

% %% Line (Trapezoid)
% disp('######################################');
% disp('#### Design a circular trajectory ####');
% disp('####                              ####');
% disp('######################################');
% disp(' ');
% kx = linspace(-5,5, 256)';
% ky = linspace(-5,5, 256)';
% kz = 0*ky;
% C = [kx ky kz];
% 
% [C_riv, time_riv, g_riv, s_riv, k_riv] = minTimeGradient(C,0, 0, 0, 4, 15, 4e-3);          % Rotationally invariant solution
% [C_rv, time_rv, g_rv, s_rv, k_rv] = minTimeGradient(C,1, 0, 0, 4, 15, 4e-3);     % Rotationally variant solution
% 
% L = max(length(s_riv), length(s_rv));
% 
% figure, subplot(2,2,1), plot(C_rv(:,1), C_rv(:,2)); title('k-space'); axis([-5 5 -5 5]);
% subplot(2,2,2), plot(g_rv(:,1), 'r'); axis([0,L,-4.5,4.5]); title('gradient waveforms (R. Variant)'); axis([0 L 0 6]);
% hold on, subplot(2,2,2), plot(g_rv(:,2), '-.');
% legend('gx', 'gy', 'Location', 'NorthEast');
% subplot(2,2,3), plot((g_riv(:,1).^2 + g_riv(:,2).^2).^0.5, '--'), 
% hold on, subplot(2,2,3), plot((g_rv(:,1).^2 + g_rv(:,2).^2).^0.5, 'r');   axis([0 L 0 6]);
% legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthEast'); title('gradient magnitude')
% subplot(2,2,4), plot((s_riv(:,1).^2 + s_riv(:,2).^2).^0.5, '--'); title('slew-rate magnitude'); axis([0 L 0 27]);
% hold on, subplot(2,2,4), plot((s_rv(:,1).^2 + s_rv(:,2).^2).^0.5, 'r'); 
% legend('rotationally invariant', 'rotationally variant', 'Location', 'NorthEast');
% 
% 
% %% Circle
% 
% disp('######################################');
% disp('#### Design a circular trajectory ####');
% disp('####                              ####');
% disp('######################################');
% disp(' ');
% 
% C = exp(i*2*pi*linspace(0,1,512)')*10;
% C = [real(C) imag(C) 0*C];
% [C_rv, time_rv, g_rv, s_rv, k_rv] = minTimeGradient(C,0);          % Rotationally variant solution
% [C_riv, time_riv, g_riv, s_riv, k_riv] = minTimeGradient(C,1, 0);  % Rotationally invariant solution
% 
% L = max(length(s_riv), length(s_rv));
% 
% figure, subplot(2,2,1), plot(C_rv(:,1), C_rv(:,2)); title('k-space'); axis([-10 10 -10 10]);
% subplot(2,2,2), plot(g_riv(:,1)); axis([0,L,-4.5,4.5]); title('gradient waveforms (R. Variant)'); axis([0 L -6 6]);
% hold on, subplot(2,2,2), plot(g_riv(:,2), 'r');
% legend('gx', 'gy', 'Location', 'NorthEast');
% subplot(2,2,3), plot((g_rv(:,1).^2 + g_rv(:,2).^2).^0.5, '--'), 
% hold on, subplot(2,2,3), plot((g_riv(:,1).^2 + g_riv(:,2).^2).^0.5, 'r'); axis([0 L 0 6]);
% legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthEast'); title('gradient magnitude')
% subplot(2,2,4), plot((s_rv(:,1).^2 + s_rv(:,2).^2).^0.5, '--'); title('slew-rate magnitude'); axis([0 L 0 20]);
% hold on, subplot(2,2,4), plot((s_riv(:,1).^2 + s_riv(:,2).^2).^0.5, 'r'); 
% legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthEast');
% 
% return;
% 
% %% Spiral
% 
% disp('######################################');
% disp('#### Design a dual density spiral ####');
% disp('####                              ####');
% disp('######################################');
% disp(' ');
% 
% [k_rv,g_rv,s_rv,time_rv,Ck_rv] = vdSpiralDesign(1, 16, 0.83,[55,55,10,10],[0,0.2,0.3,1],4,15,4e-3,[],'cubic');
% [k_riv,g_riv,s_riv,time_riv,Ck_riv] = vdSpiralDesign(0, 16, 0.83,[55,55,10,10],[0,0.2,0.3,1],4,15,4e-3,[],'cubic');
% 
% L = max(length(s_riv), length(s_rv));
% 
% figure, subplot(2,2,1), plot(k_rv(:,1), k_rv(:,2)); title('k-space'); axis([-6 6 -6 6]);
% subplot(2,2,2), plot(g_riv(:,1)); axis([0,L,-4.5,4.5]); title('gradient waveforms (R. Variant)')
% hold on, subplot(2,2,2), plot(g_riv(:,2), 'r');
% legend('gx', 'gy', 'Location', 'NorthEast');
% subplot(2,2,3), plot((g_rv(:,1).^2 + g_rv(:,2).^2).^0.5, '--'), 
% hold on, subplot(2,2,3), plot((g_riv(:,1).^2 + g_riv(:,2).^2).^0.5, 'r');  axis([0 L 0 6]);
% legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthEast'); title('gradient magnitude')
% subplot(2,2,4), plot((s_rv(:,1).^2 + s_rv(:,2).^2).^0.5, '--'); title('slew-rate magnitude');  axis([0 L 0 20]);
% hold on, subplot(2,2,4), plot((s_riv(:,1).^2 + s_riv(:,2).^2).^0.5, 'r'); 
% legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthWest');
% 
% %% Rosette
% 
% disp('############################################');
% disp('#### Design a rosette trajectory        ####');
% disp('####                                    ####');
% disp('############################################');
% disp(' ');
% 
% Gmx = 4;
% Smx = 15;
% T = 17/Gmx;
% Kmx = 6;
% w1 = 0.147*2*pi*Gmx;
% w2 = 0.087/1.02*2*pi*Gmx;
% t = 0e-3:4e-3:T;
% C = Kmx*sin(w1*t').*exp(i*w2*t');
% C = [real(C) imag(C) 0*C];
% 
% [C_riv, time_riv, g_riv, s_riv, k_riv] = minTimeGradient(C,0);          % Rotationally invariant solution
% [C_rv, time_rv, g_rv, s_rv, k_rv]= minTimeGradient(C,1, 0);  % Rotationally variant solution
% L = max(length(s_riv), length(s_rv));
% 
% figure, subplot(2,2,1), plot(C_rv(:,1), C_rv(:,2)); title('k-space'); axis([-6 6 -6 6]);
% subplot(2,2,2), plot(g_riv(:,1)); axis([0,L,-4.5,4.5]); title('gradient waveforms (R. Variant)')
% hold on, subplot(2,2,2), plot(g_riv(:,2), 'r');
% legend('gx', 'gy', 'Location', 'NorthEast');
% subplot(2,2,3), plot((g_riv(:,1).^2 + g_riv(:,2).^2).^0.5, '--'), 
% hold on, subplot(2,2,3), plot((g_rv(:,1).^2 + g_rv(:,2).^2).^0.5, 'r');  axis([0 L 0 6]);
% legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthEast'); title('gradient magnitude')
% subplot(2,2,4), plot((s_riv(:,1).^2 + s_riv(:,2).^2).^0.5, '--'); title('slew-rate magnitude');  axis([0 L 0 20]);
% hold on, subplot(2,2,4), plot((s_rv(:,1).^2 + s_rv(:,2).^2).^0.5, 'r'); 
% legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthWest');
% %% Cone
% 
% disp('############################################');
% disp('####   Design a cone trajectory         ####');
% disp('####                                    ####');
% disp('############################################');
% disp(' ');
% 
% r = linspace(0,5, 512)';
% th = linspace(0,2*pi, 512)';
% C = r.*exp(3*1i*th);
% C = [real(C) imag(C) r];
% figure, plot3(C(:,1), C(:,2), C(:,3))
% title('k-space trajectory')
% xlabel('k_x'); ylabel('k_y'); zlabel('k_z');
% 
% [C_rv, time_rv, g_rv, s_rv, k_rv] = minTimeGradient(C,0);          % Rotationally variant solution
% [C_riv, time_riv, g_riv, s_riv, k_riv] = minTimeGradient(C,1, 0);  % Rotationally invariant solution
% L = max(length(s_riv), length(s_rv));
% 
% figure, subplot(2,2,1), plot3(C_rv(:,1), C_rv(:,2), C_rv(:,3)); title('k-space'); axis([-6 6 -6 6]);
% subplot(2,2,2), plot(g_riv(:,1)); axis([0,L,-4.5,4.5]); title('gradient waveforms (R. Invariant)')
% hold on, subplot(2,2,2), plot(g_riv(:,2), 'r');
% hold on, subplot(2,2,2), plot(g_riv(:,3), 'g');
% legend('gx', 'gy', 'gz', 'Location', 'NorthEast');
% subplot(2,2,3), plot((g_rv(:,1).^2 + g_rv(:,2).^2).^0.5, '--'), 
% hold on, subplot(2,2,3), plot((g_riv(:,1).^2 + g_riv(:,2).^2).^0.5, 'r');  axis([0 L 0 6]);
% legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthEast'); title('gradient magnitude')
% subplot(2,2,4), plot((s_rv(:,1).^2 + s_rv(:,2).^2+ s_rv(:,3).^2).^0.5, '--'); title('slew-rate magnitude');  axis([0 L 0 20]);
% hold on, subplot(2,2,4), plot((s_riv(:,1).^2 + s_riv(:,2).^2+ s_riv(:,3).^2).^0.5, 'r'); 
legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthWest');