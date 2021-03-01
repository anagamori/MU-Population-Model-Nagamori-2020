%==========================================================================
% run_MU_population_model.m
% Author: Akira Nagamori
% Last update: 7/10/2020
% Descriptions:
%   Run simulations of a motor unit population model
%==========================================================================

load('modelParameter')

%% Run simulation
Fs = 10000; % sampling frequency
time = 0:1/Fs:5; % time vector

amp = 1; % 10% of maximum synaptic input
synaptic_input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];

tic
output = MU_population_model(Fs,time,synaptic_input,modelParameter,1);
toc
% output = MU_population_model_no_tendon(Fs,time,input,modelParameter,1); %
% simulation for no tendon condition in Fig.6-8


