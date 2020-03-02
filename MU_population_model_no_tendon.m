%==========================================================================
% MU_population_model_no_tendon.m
% Author: Akira Nagamori
% Last update: 3/1/20
% Model descriptions:
%   input: Fs = sampling frequency
%             time = time vector
%             synaptic_drive = synaptic drive to a motor unit population 
%             modelParameter = a structure that contains all model
%             parameters
%             fitOpt = figure display option (0: no figure output, 1: display
%             figures)
%   output: a data structure ('output') that contains spike trains, motor
%   unit forces and tendon force (see l.242)
%==========================================================================

function [output] = MU_population_model_no_tendon(Fs,time,synaptic_drive,modelParameter,figOpt)
%% Muscle architectural parameters
Lce = 1;
Vce = 0;

%% Motor unit architecture
N_MU = modelParameter.N_MU; % number of motor units
i_MU = modelParameter.i_MU; % index for motor units
index_slow = modelParameter.index_slow; % index for slow-twitch units
PTi = modelParameter.PTi; % peak tetanic force [N]
U_th = modelParameter.U_th; % recruitment threshold [0-1]
FR_half = modelParameter.FR_half; % f_0.5 (frequency at which a motor unit reaches half the maximum force)
MDR = modelParameter.MDR; % minimum discharge rate
PDR = modelParameter.PDR; % peak discharge rate

% coefficients for Eq. 2-7
g_e = modelParameter.g_e;
index_saturation = modelParameter.index_saturation;
lamda = modelParameter.lamda;
k_e = modelParameter.k_e;
U_th_t = modelParameter.U_th_t;

% matrix of coefficients used for Eq. 11, 13-15, 20
parameter_Matrix = modelParameter.parameterMatrix;

% coefficients for Eq. 12
tau_1 = parameter_Matrix(:,9);
R_temp = exp(-time./tau_1);
gamma = parameter_Matrix(:,15);

% parameters for Eq. 8
cv_MU = modelParameter.CV_MU; % coefficient of variation for interspike intervals

%% Initilization
spike_time = zeros(N_MU,1);
spike_train = zeros(N_MU,length(time));
force = zeros(N_MU,length(time));
Force = zeros(1,length(time));

R = zeros(N_MU,length(time));
c = zeros(N_MU,1);
cf = zeros(N_MU,1);
A = zeros(N_MU,1);
c_mat = zeros(N_MU,length(time));
cf_mat = zeros(N_MU,length(time));
A_tilde_mat = zeros(N_MU,length(time));
A_mat = zeros(N_MU,length(time));

S_i = zeros(N_MU,1);
Y_i = zeros(N_MU,1);
S_mat = zeros(N_MU,length(time));
Y_mat = zeros(N_MU,length(time));

FL = zeros(N_MU,1);
FV = zeros(N_MU,1);

DR_temp = zeros(N_MU,1);
DR_mat = zeros(N_MU,length(time));
%% Simulation
rng('shuffle')
for t = 1:length(time)
    
    if t > 1
        %% Calculate discharge rate (Eq.2-7)
        U_eff = synaptic_drive(t); % effective synaptic input 
        DR_MU = g_e.*(U_eff-U_th)+MDR;  % compute discharge rates of fast-twitch units
        for m = 1:length(index_saturation) % go through all slow-twitch units
            index = index_saturation(m);
            if U_eff <= U_th_t(index)
                DR_temp(index) = MDR(index) + lamda(index).*k_e(index)*(U_eff-U_th(index));
            else
                DR_temp(index) = PDR(index)-k_e(index)*(1-U_eff);
            end
        end
        DR_MU(index_saturation) = DR_temp(index_saturation); 
        DR_MU(DR_MU<MDR) = 0; % zero out discharge rate < minimum discharge rate
        DR_MU(DR_MU>PDR) = PDR(DR_MU>PDR); % set discharge rate = peak discharge rate if discharge rate > peak discharge rate
        DR_mat(:,t) = DR_MU;
        %% Sag & Yield (Eq.16-18) (Song et al., 2008)
        f_eff = DR_MU./FR_half;
        S_i = sag_function(S_i,f_eff,a_s,Fs);
        S_i(1:index_slow) = 1;
        S_mat(:,t) = S_i;
        Y_i = yield_function(Y_i,Vce,Fs);
        Y_i(index_slow+1:end) = 1;
        Y_mat(:,t) = Y_i;
        
        %% Convert activation into spike trains
        index_1 = i_MU(DR_MU >= MDR & DR_mat(:,t-1) == 0);
        index_2 = i_MU(DR_MU >= MDR & spike_time ==t);
        index = [index_1;index_2];
        
        for j = 1:length(index) % loop through motor units whose firing rate is greater than minimum firing rate defined by the user
            n = index(j);
            spike_train_temp = zeros(1,length(t));
            if ~any(spike_train(n,:)) % when the motor unit fires at the first time
                spike_train(n,t) = 1; % add a spike to the vector
                spike_train_temp(t) = 1;
                mu = 1/DR_MU(n); % interspike interval
                Z = randn(1);
                Z(Z>3.9) = 3.9;
                Z(Z<-3.9) = -3.9;
                spike_time_temp = (mu + mu*cv_MU*Z)*Fs;  % compute next spike time
                if spike_time_temp <= 0.002*Fs
                    spike_time_temp = 0.002*Fs;
                end
                spike_time(n) = round(spike_time_temp) + t; % save spike time
                
                % convert spike train into R in Eq. 12
                temp = conv(spike_train_temp,R_temp(n,:)*(1+2*A(n)^gamma(n)));
                R(n,:) = R(n,:) + temp(1:length(time));
            else % when the motor unit have already fired at least once
                if spike_time(n) == t % when the motor unit fires
                    spike_train(n,t) = 1;
                    spike_train_temp(t) = 1;    
                    mu = 1/DR_MU(n); % interspike interval
                    Z = randn(1);
                    Z(Z>3.9) = 3.9;
                    Z(Z<-3.9) = -3.9;
                    spike_time_temp = (mu + mu*cv_MU*Z)*Fs; % compute next spike time
                    if spike_time_temp <= 0.002*Fs
                        spike_time_temp = 0.002*Fs;
                    end
                    spike_time(n) = round(spike_time_temp) + t; % save spike time
                    
                    % convert spike train into R in Eq. 12
                    temp = conv(spike_train_temp,R_temp(n,:)*(1+2*A(n)^gamma(n)));
                    R(n,:) = R(n,:) + temp(1:length(time));
                elseif t > spike_time(n) + round(1/DR_MU(n)*Fs)
                    spike_train(n,t) = 1;
                    spike_train_temp(t) = 1;
                    spike_time(n) = t;
                    mu = 1/DR_MU(n); % interspike interval
                    Z = randn(1);
                    Z(Z>3.9) = 3.9;
                    Z(Z<-3.9) = -3.9;
                    spike_time_temp = (mu + mu*cv_MU*Z)*Fs; % compute next spike time
                    if spike_time_temp <= 0.002*Fs
                        spike_time_temp = 0.002*Fs;
                    end
                    spike_time(n) = round(spike_time_temp) + t; % save spike time
                    
                    % convert spike train into R in Eq. 12
                    temp = conv(spike_train_temp,R_temp(n,:)*(1+2*A(n)^gamma(n)));
                    R(n,:) = R(n,:) + temp(1:length(time));
                end
            end
        end
        
        %% Convert spikes into activation (Eq. 10-19)
        [c,cf,A_tilde,A] = spike2activation(R(:,t),c,cf,A,parameter_Matrix,Lce,S_i,Y_i,Fs);
        
        c_mat(:,t) = c;
        cf_mat(:,t) = cf;
        A_tilde_mat(:,t) = A_tilde;
        A_mat(:,t) = A;
        
        %% Force-length and force-velocity
        FL(1:index_slow) = FL_slow_function(Lce);
        FL(index_slow+1:end) = FL_fast_function(Lce);
        
        if Vce > 0
            FV(1:index_slow) = FVecc_slow_function(Lce,Vce);
            FV(index_slow+1:end) = FVecc_fast_function(Lce,Vce);
        else
            FV(1:index_slow) = FVcon_slow_function(Lce,Vce);
            FV(index_slow+1:end) = FVcon_fast_function(Lce,Vce);
        end
        %% Muscle force output
        f_i = A.*PTi.*FL.*FV;
        force(:,t) = f_i; % motor unit force
        
        Force(t) = sum(f_i); % muscle force
    end
end

%% Plotting figures
if figOpt == 1
    figure(1)
    plot(time,Force)
    xlabel('Time (s)')
    ylabel('Force (N)')
    hold on
end

%% Save data in the data structure, output
output.Force = Force; % muscle force
output.force = force; % individual motor unit forces
output.spike_train = spike_train; % spike trains

%% Convert spike trian into activation
    function [c,cf,A_tilde,A] = spike2activation(R,c,cf,A,parameter_Matrix,Lce,S_i,Y_i,Fs)
        %% Stage 1: calcium kinetics (Eq.10-11)
        % parameters
        S = parameter_Matrix(:,1); 
        C = parameter_Matrix(:,2); 
        k_1 = parameter_Matrix(:,3);
        k_2 = parameter_Matrix(:,4); 
        k_3 = parameter_Matrix(:,5)*Lce + parameter_Matrix(:,6); 
        k_4 = parameter_Matrix(:,7)*Lce + parameter_Matrix(:,8); 
        tau_2 = parameter_Matrix(:,10);
        N = parameter_Matrix(:,11)*Lce + parameter_Matrix(:,12); 
        K = parameter_Matrix(:,13)*Lce + parameter_Matrix(:,14); 
        
        c_dot = k_1.*(C-c-cf).*R - k_2.*c.*(S-C+c+cf)-(k_3.*c-k_4.*cf).*(1-cf);
        cf_dot = (1-cf).*(k_3.*c-k_4.*cf);
        c = c_dot/Fs + c;
        cf = cf_dot/Fs + cf;
        
        %% Stage 2: cooperativity and saturation + sag and yield (Eq.15-18)
        if cf < 0
            cf_temp = 0;
        else
            cf_temp = cf.*S_i.*Y_i;
        end
        A_tilde = cf_temp.^N./(cf_temp.^N+K.^N);
        
        %% Stage 3: first-order dynamics to muscle activation, A (Eq. 19)
        A_dot = (A_tilde-A)./tau_2;
        A = A_dot./Fs + A;
        
    end

    %% Sag  (Eq.16-17)
    function [S] = sag_function(S,f_eff,a_s,Fs)
        
        a_s(f_eff<0.1) = 20;
        
        T_s = 0.015;
        S_dot = (a_s - S)./T_s;
        S = S_dot/Fs + S;
        
    end

    %% Yield (Eq.18)
    function [Y] = yield_function(Y,V,Fs)
        c_y = 0.35;
        V_y = 0.1;
        T_y = 0.2;
        Y_dot = (1-c_y.*(1-exp(-abs(V)./V_y))-Y)./T_y;
        Y = Y_dot/Fs + Y;
        
    end

    %% Force-length relationship for slow twitch (Eq.24)
    function FL = FL_slow_function(L)
        beta = 2.3;
        omega = 1.12;
        rho = 1.62;
        
        FL = exp(-abs((L^beta - 1)/omega)^rho);
    end

    %% Force-length relationship for fast twitch (Eq.24)
    function FL = FL_fast_function(L)
     
        beta = 1.55;
        omega = 0.75;
        rho = 2.12;
        
        FL = exp(-abs((L^beta - 1)/omega)^rho);
    end

    %% Concentric force-velocity relationship for slow twitch (Eq.25)
    function FVcon = FVcon_slow_function(L,V)
       
        Vmax = -7.88;
        cv0 = 5.88;
        cv1 = 0;
        
        FVcon = (Vmax - V)/(Vmax + (cv0 + cv1*L)*V);
    end

    %% Concentric force-velocity relationship for fast twitch (Eq.25)
    function FVcon = FVcon_fast_function(L,V)
       
        Vmax = -9.15;
        cv0 = -5.7;
        cv1 = 9.18;
        
        FVcon = (Vmax - V)/(Vmax + (cv0 + cv1*L)*V);
    end

    %% Eccentric force-velocity relationship for slow twitch (Eq.25)
    function FVecc = FVecc_slow_function(L,V)
        
        av0 = -4.7;
        av1 = 8.41;
        av2 = -5.34;
        bv = 0.35;
        FVecc = (bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end

    %% Eccentric force-velocity relationship for slow twitch  (Eq.25)
    function FVecc = FVecc_fast_function(L,V)
       
        av0 = -1.53;
        av1 = 0;
        av2 = 0;
        bv = 0.69;
        FVecc = (bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end
end