%==========================================================================
% run_MU_population_model.m
% Author: Akira Nagamori
% Last update: 3/1/19
% Descriptions:
%   Run simulations of a motor unit population model
%==========================================================================

load('modelParameter')
modelParameter.CV_MU = 0.2;

%% Run simulation
Fs = 30000; % sampling frequency
time = 0:1/Fs:5; % time vector

amp = 0.1; % 10% of maximum synaptic input
synaptic_drive = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];

output = MU_population_model(Fs,time,synaptic_drive,modelParameter,1);
% output = MU_population_model_no_tendon(Fs,time,input,modelParameter,1); %
% simulation for no tendon condition in Fig.6-8


