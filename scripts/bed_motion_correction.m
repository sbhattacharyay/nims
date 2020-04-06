%% Authors: Shubhayu Bhattacharyay, B.S. Candidate, Matthew Wang, B.S. Candidate
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
% March 2019; Last revision: 06-Apr-2020
%% ------------- BEGIN CODE --------------
%Set Directory and Procure all feature data
addpath('functions/')
tic
load('../motion_feature_data/imputed_complete_sensor_data.mat');
load('../motion_feature_data/feature_thresholds.mat');
toc

%% Bed Motion Correction Algorithm:
% 1. In SMA domain, locate periods of time when bed contributes significant
% motion
% 2. Save those indices as markers for bed motion
% 3. Analyze SMA profiles directly preceeding those markers. If any motion
% sensor's SMA profile is significantly active prior to bed motion, then we
% 4. If no motion directly preceeds bed motion in those sensors, subtract
% bed SMA from those points.