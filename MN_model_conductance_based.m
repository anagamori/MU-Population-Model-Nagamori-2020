%==========================================================================
% MN_model_conductance_based.m
% Author: Akira Nagamori
% Last update: 3/3/20
% Model descriptions:
%   input: Fs = sampling frequency
%             time = time vector
%             input_exc = excitatory synaptic conductance
%             input_inh = inhibitory synaptic conductance
%             Fs = sampling frequency 
%             pltOpt = plotting option (0: no plot generated, 1: plot generated)
%             PIC = PIC amplitude [nA]
%   output: a data structure ('output') that contains binary spike trains, somatic voltage, and effective synaptic current (see l.177)
%==========================================================================
function output = MN_model_conductance_based(time,input_exc,input_inh,Fs,pltOpt,PIC)

step = 1/Fs;

%% Synaptic reversal potentials
E_exc = 70; 
E_inh = -16;
%% Geometric parameters
R_i = 0.07; %[kOhm*cm]
param_s.area_s = 5.943e-5; % the surface area of soma [cm^2]
r_s = sqrt(param_s.area_s/(2*pi)); % radius [cm]
l_s = sqrt(param_s.area_s/(2*pi)); % length [cm]

param_d.area_d = param_s.area_s*38; % the surface area of dendrite [cm^2]
r_d = r_s*0.5; % radius [cm]
l_d = param_d.area_d./(2*pi*r_d); % length [cm]

g_c = 2./(R_i.*l_d./(pi.*r_d.^2)+R_i.*l_s./(pi.*r_s.^2)); % coupling conductance [mS]

%% Dendtritic parameters
param_d.V_l = 0; % leak reversal potential [mV]
param_d.g_c = g_c; % coupling conductance [mS]
param_d.g_l = param_d.area_d./15; % leak conductance [mS]
param_d.C_m = 1*param_d.area_d; % membrane capacitance [microF]

%% Somatic parameters
param_s.V_Na = 120; % sodium reversal potential [mV]
param_s.V_K = -10; % potassium reversal potential [mV]
param_s.V_l = 0; % leak reversal potential [mV]

param_s.g_Na = 30.*param_s.area_s ; % sodium conductance [mS]
param_s.g_Kf = 4.*param_s.area_s ; % fast potassium conductance [mS]
param_s.g_Ks = 167.5861.*param_s.area_s ; % slow potassium conductance [mS]
param_s.g_c = g_c; % coupling conductance [mS]
param_s.g_l = param_s.area_s./0.6; % leak conductance [mS]
param_s.C_m = 1.*param_s.area_s; % membrane capacitance [microF]

param_s.V_th = 4.3834e6.*9.642e-10*10^3; % voltage threshold [mV]

%% State variables
alpha_m = 22000;
beta_m = 13000;
alpha_h = 500;
beta_h = 4000;
alpha_n = 1500;
beta_n = 100;
alpha_q = 1.6884e+03;
beta_q = 24.2778;

%% Initialize variables
V_s = 0;
V_d = 0;

m = 0;
h = 1;
n = 0;
q = 0;

m_0 = 0;
h_0 = 1;
n_0 = 0;
q_0 = 0;
%% Vectors to store data
V_s_vec =  zeros(1,length(time));
V_d_vec =  zeros(1,length(time));
m_vec =  zeros(1,length(time));
h_vec =  zeros(1,length(time));
n_vec =  zeros(1,length(time));
q_vec =  zeros(1,length(time));

binary = zeros(1,length(time));

I_eff_vec = zeros(1,length(time));

index_th = zeros(1,length(time));
index_s = zeros(1,length(time));
t_0_end = 0;

flag = 0;
%%
for t = 1:length(time)
    %% effective current    
    if flag == 0
        I_eff = sum(input_exc(:,t)*(E_exc-V_d)) + sum(input_inh(:,t)*(E_inh-V_d));
    elseif flag == 1
        I_eff = sum(input_exc(:,t)*(E_exc-V_d)) + sum(input_inh(:,t)*(E_inh-V_d)) + PIC*10^-3;
    end
    
    %%
    k_1_s = f_dV_s(V_s,V_d,m,h,n,q,param_s);
    y_1_s = V_s + k_1_s*step/2;
    k_2_s = f_dV_s(y_1_s,V_d,m,h,n,q,param_s);
    y_2_s = V_s + k_2_s*step/2;
    k_3_s = f_dV_s(y_2_s,V_d,m,h,n,q,param_s);
    y_3_s = V_s + k_3_s*step/2;
    k_4_s = f_dV_s(y_3_s,V_d,m,h,n,q,param_s);
    V_s = V_s + step/6*(k_1_s+2*k_2_s+2*k_3_s+k_4_s);
    
    %%
    k_1_d = f_dV_d(V_d,V_s,param_d,I_eff);
    y_1_d = V_d + k_1_d*step/2;
    k_2_d = f_dV_d(y_1_d,V_s,param_d,I_eff);
    y_2_d = V_d + k_2_d*step/2;
    k_3_d = f_dV_d(y_2_d,V_s,param_d,I_eff);
    y_3_d = V_d + k_3_d*step/2;
    k_4_d = f_dV_d(y_3_d,V_s,param_d,I_eff);
    V_d = V_d + step/6*(k_1_d+2*k_2_d+2*k_3_d+k_4_d);
    
    %% Spike detection and pulse generation
    V_th =  param_s.V_th;
    if V_s >= V_th
        index_th(t) = 1;
    end
    if t > 1
        if index_th(t) ==1  && index_th(t-1) ==0 && index_s(:,t) == 0
            binary(t) = 1;
            index_s(t:t+round(0.0006*Fs)) = 1; %pulse duration of 0.6 ms
            flag = 1;
        end
    end
    %% Update state variables for ionic conductances
    if t > 1
        if index_s(t) == 1 && index_s(t-1) == 0
            t_0_start = t/Fs;
            t_0_end = (t+round(0.0006*Fs))/Fs;
            m_0 = m;
            h_0 = h;
            n_0 = n;
            q_0 = q;
        elseif index_s(t) == 0 && index_s(t-1) == 1
            m_0 = m;
            h_0 = h;
            n_0 = n;
            q_0 = q;
        end
        if index_s(t) == 1
            T = t/Fs;
            m = 1 + (m_0-1)*exp(-alpha_m*(T-t_0_start));
            h = h_0*exp(-beta_h*(T-t_0_start));
            n = 1 + (n_0-1)*exp(-alpha_n*(T-t_0_start));
            q = 1 + (q_0-1)*exp(-alpha_q*(T-t_0_start));
        else
            T = t/Fs;
            m = m_0*exp(-beta_m*(T-t_0_end));
            h = 1 + (h_0-1)*exp(-alpha_h*(T-t_0_end));
            n = n_0*exp(-beta_n*(T-t_0_end));
            q = q_0*exp(-beta_q*(T-t_0_end));
        end
    end
    
    %% Store variables
    V_s_vec(t) = V_s;
    V_d_vec(t) = V_d;
    m_vec(t) =  m;
    h_vec(t) =  h;
    n_vec(t) =  n;
    q_vec(t) =  q;
       
    I_eff_vec(t) = I_eff;
    
end

%% Outout data
output.binary = binary;
output.V_s = V_s_vec;
output.I_eff = I_eff_vec;
%%
if pltOpt == 1
    figure()
    ax1 = subplot(3,1,1);
    plot(time,I_eff_vec,'LineWidth',1,'Color','k')
    xlabel('Time (s)')
    ylabel({'Effective current';'(nA)'})
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax2 = subplot(3,1,2);
    plot(time,V_s_vec,'LineWidth',1,'Color','k')
    hold on
    plot(time,binary*50,'LineWidth',1,'Color','b')
    xlabel('Time (s)')
    ylabel({'Soma potential';'(mV)'})
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax3 = subplot(3,1,3);
    plot(time,V_d_vec,'LineWidth',1,'Color','k')
    xlabel('Time (s)')
    ylabel({'Dendritic potential';'(mV)'})
    set(gca,'TickDir','out');
    set(gca,'box','off')
    linkaxes([ax1,ax2,ax3],'x')
end
    %% Eq 32
    function dx = f_dV_s(V_s,V_d,m,h,n,q,param_s)
        I_Na = param_s.g_Na*m^3*h*(V_s-param_s.V_Na); % sodium channel
        I_Kf = param_s.g_Kf*n^4*(V_s-param_s.V_K); %  fast potassium delayed rectifier 
        I_Ks = param_s.g_Ks*q^2*(V_s-param_s.V_K); %  slow potassium delayed rectifier
        
        I_l = param_s. g_l*(V_s-param_s.V_l); % leak channel
        I_c = param_s.g_c*(V_s-V_d); % coupling conductance
        dx = 1/param_s.C_m*(-I_Na - I_Kf - I_Ks...
            - I_l - I_c)*1000; %
        
    end
    
    %% Eq 31
    function dx = f_dV_d(V_d,V_s,param_d,I_eff)
        
        I_l = param_d.g_l*(V_d-param_d.V_l); % leak channel
        I_c = param_d.g_c*(V_d-V_s); % coupling conductance
      
        dx = 1/param_d.C_m*(- I_l - I_c + I_eff)*1000; 
              
    end
end