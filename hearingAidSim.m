% EC4450        Array Signal Processing
%               Naval Postgraduate School, Monterey CA
%% Title:       hearingAidSim_v3.m
% Students:     Teddy Herrera
% Desciption:   Initializes a simulation that models FIR based wideband
% MPDR beamforming for a behind the ear hearing aid Microphone Array
% %%----Changes from Previous Version----%%
% None

close all; clear; clc;

%% --- Parameters ---
c = 343;              % Speed of sound [m/s]
fs = 48e3;           % Sampling frequency [Hz]
L = 256;               % FIR length (space-time)
num_samples = 24e3;   % Example signal length
M = 3;                % Number of microphones

%% -- Define Array Geometry (Kayser et al. 2009) --
% Define array geometry with the center of the ear canal 
% being at [0,0]. The +y axis is coming out of the ear canal, the positive x axis is
% the horizontal direction from the ear to the front of the head, and the 
% + z-axis is the vertical direction from the ear to the top of the head.
% Distance between element 1 and 2 is 7.6 mm
% Distance between 2 and 3 is 7.3 mm

% element positons [x; y; z]
pos_BTE = [-13.6,  -21.2,  -28.5;
            zeros(1,3);
            37.4,   35.3,   32.7] * 1e-3; % element 3

%% -- Target Direction 
theta_t = pi/2; % 90 degree elevation -> xy plane
psi_t = pi/2; % 90 degree azimuth -> y-axis (broadside)

u_t = [ sin(theta_t) * cos(psi_t);
      sin(theta_t) * sin(psi_t);
      cos(theta_t) ];  

delays_t = (pos_BTE' * u_t) / c;

%% -- Generate Target Space-Time Steering Vector g
g_target = build_space_time_steering(M, L, delays_t, fs);

%% -- Null Direction
num_nulls = 0;  % Sets # of nulls
theta_nulls = [pi/2 * ones(2,1); pi/4];
psi_nulls = [deg2rad(30); deg2rad(140); pi/2];

%% -- Generate All Null Space-Time Steering Vector g
G_nulls = zeros(M * L, num_nulls);  % Each column is a null constraint
for n = 1:num_nulls
    u_null = [sin(theta_nulls(n)) * cos(psi_nulls(n));
              sin(theta_nulls(n)) * sin(psi_nulls(n));
              cos(theta_nulls(n))];
    delays_null = (pos_BTE' * u_null) / c;
    G_nulls(:, n) = build_space_time_steering(M, L, delays_null, fs);
end

%% -- Build Covariance Matrix
% For demo: white noise 
Rxx = eye(M * L);

%% -- Build Constraint Matrix
if num_nulls == 0
    C = g_target;
    f = 1;
else
    C = [g_target, G_nulls];
    f = [1; zeros(num_nulls,1)];
end
%% -- Compute LCMV Weights
w = (Rxx \ C) / (C' * (Rxx \ C)) * f;

%% ---Azimuth & Elevation Scan to Plot Pattern ---
theta = pi/2;  % Elevation fixed at 90° (x-y plane)
psi_scan = linspace(0, pi, 180);
beamAz = zeros(size(psi_scan));

for idx = 1:length(psi_scan)
    dir = [sin(theta)*cos(psi_scan(idx));
           sin(theta)*sin(psi_scan(idx));
           cos(theta)];
    delays_scan = (pos_BTE' * dir) / c;
    g_AzScan = build_space_time_steering(M, L, delays_scan, fs);
    beamAz(idx) = abs(w' * g_AzScan);
end

beamAz_dB = 20*log10(abs(beamAz) / max(abs(beamAz)));

theta_scan = linspace(0, pi, 180);  
psi = pi/2;
beamEl = zeros(size(theta_scan));

for idx = 1:length(theta_scan)
    dir = [sin(theta_scan(idx))*cos(psi);
           sin(theta_scan(idx))*sin(psi);
           cos(theta_scan(idx))];
    delays_scan = (pos_BTE' * dir) / c;
    g_ElScan = build_space_time_steering(M, L, delays_scan, fs);
    beamEl(idx) = abs(w' * g_ElScan);
end

beamEl_dB = 20*log10(abs(beamEl) / max(abs(beamEl)));

%% -- Plot the Beampattern
figure;
plot(rad2deg(psi_scan), beamAz_dB);
title('LCMV Azimuth Beampattern', 'Null at: \theta = [90; 90; 45] deg \psi = [30; 140; 90] deg');
ylim([-70 0]);

figure;
plot(rad2deg(theta_scan), beamEl_dB);
title('LCMV Elevation Beampattern', 'Null at: \theta = [90; 90; 45] deg \psi = [30; 140; 90] deg');
ylim([-70 0]);

function g = build_space_time_steering(M, L, delays, fs)
    g = zeros(M * L, 1);
    for m = 1:M
        delay_samples = delays(m) * fs;
        h_d = sinc((0:L-1) - delay_samples);  % Fractional delay approximation
        g((m-1)*L + (1:L)) = h_d;
    end
end