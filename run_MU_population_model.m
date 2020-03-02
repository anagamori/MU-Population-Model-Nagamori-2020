%==========================================================================
% run_MU_population_model.m
% Author: Akira Nagamori
% Last update: 3/1/19
% Descriptions:
%   Run simulations of a motor unit population model
%==========================================================================
close all
clear all
clc


%%
data_folder = '/Volumes/DATA2/New_Model/withTendon/Model_8_drift';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_8';

%%
cd(model_parameter_folder)
load('modelParameter')
cd(code_folder)
%% MU simulation parameters
modelParameter.CV_MU = 0.2;

Fs = 30000; % sampling frequency
time = 0:1/Fs:15; % time vector

amp = 0.1; % 10% of maximum synaptic input
input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
output = MU_population_model(Fs,time,input,modelParameter,1);


