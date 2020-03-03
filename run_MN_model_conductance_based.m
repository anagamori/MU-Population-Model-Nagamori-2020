%==========================================================================
% run_MN_model_conductance_based.m
% Author: Akira Nagamori
% Last update: 3/3/20
% Descriptions:
%   Generate conductance-based synaptic input based on stochastic
%   excitatory and inhibitory EPSP and IPSPs
%   Run a conductance based MN model (MN_model_conductance_based)
%==========================================================================

Fs = 30000; % sampling frequency
time = 0:1/Fs:5; % time vector

%% generate vectors of excitatory and inhibitory synaptic conductances (Eq.48)
f = 1; % discharge rate of pre-synaptic neurons [Hz]
g_max = 0.0001; % maximum conductance [mS]
n_exc = 100;
n_inh = 100;
r_exc = zeros(n_exc,length(0:1/Fs:1));
r_inh = zeros(n_exc,length(0:1/Fs:1));
% 
for k = 1:4 % create four 1-sec vectors and concatenate them 
    % generate vectors of normalized synaptic conductance of n-neurons, each
    % discharging a specified rate, f, with duration of 1-sec
    [input_exc_temp,input_inh_temp] = input_generator(n_exc,n_inh,[1:1*Fs]./Fs,Fs,1,f);
    r_exc = [r_exc input_exc_temp];
    r_inh = [r_inh input_inh_temp];
end
% turn cumulative spike trains into synaptic conductance 
r_exc = r_exc*g_max;
r_inh = r_inh*g_max;
%% Run the model 
pltOpt = 0;
output = MN_model_conductance_based(time,r_exc,r_inh,Fs,pltOpt,1);

%% Compute statistical properties of synaptic input
I_eff = output.I_eff*10^3; % convert the unit into nA
mean_I = mean(I_eff(2*Fs+1:end)); % mean synaptic current
var_I = var(I_eff(2*Fs+1:end)); % variance of synaptic current
std_I = std(I_eff(2*Fs+1:end)); % standard deviation 

%% Compute discharge variability of a motoneuron
spike_time = find(output.binary(2*Fs:end));
ISI = diff(spike_time)/(Fs/1000); % vector of inter-spike intervals (ISIs)
mean_DR = mean(1./ISI*1000); % mean discharge rate
CoV_ISI = std(1./ISI*1000)/mean_DR*100; % coefficient of variation of ISIs

%% input generator function
% generate a vector of cumulative spike trains of pre-synaptic neurons,
% each of which follow a Poisson process with a specfied dischage rate
function [r_exc,r_inh] = input_generator(n_exc,n_inh,time,Fs,tStim,f)

dt = 1/Fs; % time step
tVec = 0: dt:tStim-dt; % length of signal
r_exc = zeros(n_exc,length(time));
r_inh = zeros(n_inh,length(time));

for i = 1:n_exc    
    nBins = floor(tStim/dt); % number of samples/bins in the specified length of time
    spike = rand(1,nBins) < f*dt; % generate spike index
    
    r_exc_temp  = zeros(1,length(time)); % assign it in a vector
    r_exc_temp(end-length(tVec)+1:end) = spike;
    r_exc(i,:) = kinetic_model(time,r_exc_temp,Fs);   
end

for i = 1:n_inh
    nBins = floor(tStim/dt);
    spike = rand(1,nBins) < f*dt;
    
    r_inh_temp  = zeros(1,length(time));
    r_inh_temp(end-length(tVec)+1:end) = spike;
    r_inh(i,:) = kinetic_model(time,r_inh_temp,Fs);
    
end
end

%% Conductance kinetics described in Destexhe et al. 1994
function r_vec = kinetic_model(time,input,Fs)

% model parameters (Destexhe et al. 1994)
alpha = 2e3;
beta = 1e3;
T_max = 1;
r_inf = alpha*T_max/(alpha*T_max+beta);
tau_r = 1/(alpha*T_max+beta);

% initialization
index_s = zeros(1,length(time));
r_vec = zeros(1,length(time));

t_0_end = 0;
r = 0;
r_0 = 0;

for t = 1:length(time)
    if t > 1
        if input(t) == 1
            index_s(t:t+round(0.0006*Fs)) = 1; %pulse duration of .6 ms
        end
    end
    %% State variables for ionic conductances
    if t > 1
        if index_s(t) == 1 && index_s(t-1) == 0
            t_0_start = t/Fs;
            t_0_end = (t+round(0.0006*Fs))/Fs;
            r_0 = r;
            
        elseif index_s(t) == 0 && index_s(t-1) == 1
            r_0 = r;
            
        end
        if index_s(t) == 1
            T = t/Fs;
            r = r_inf + (r_0-r_inf)*exp(-(T-t_0_start)/tau_r);
            
        else
            T = t/Fs;
            r = r_0*exp(-beta*(T-t_0_end));
        end
        
    end
    r_vec(t) = r;
end

end
