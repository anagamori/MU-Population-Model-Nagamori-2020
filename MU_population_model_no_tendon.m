%==========================================================================
% muscleModel_noTendon.m
% Author: Akira Nagamori
% Last update: 7/10/2020
% Descriptions:
%   Full model without tendon
%==========================================================================
function [output] = MU_population_model_noTendon(Fs,time,synaptic_input,modelParameter,figOpt)

%% Motor unit architecture
N_MU = modelParameter.N_MU; % number of motor units
i_MU = modelParameter.i_MU; % index for motor units
index_slow = modelParameter.index_slow;

%% Peak tetanic force
PTi = modelParameter.PTi;

%% Recruitment threshold
U_th = modelParameter.U_th;

%% Minimum and maximum firing rate
FR_half = modelParameter.FR_half;
MDR = modelParameter.MDR;
PDR = modelParameter.PDR;

g_e = modelParameter.g_e;
index_saturation = modelParameter.index_saturation;
lamda = modelParameter.lamda;
k_e = modelParameter.k_e;
U_th_t = modelParameter.U_th_t;

Z = randn(N_MU,length(time));
Z(Z>3.9) = 3.9;
Z(Z<-3.9) = -3.9;

%% Motor unit parameters
parameter_Matrix = modelParameter.parameterMatrix;

%% Activation dynamics (Song et al., 2008)
tau_1 = parameter_Matrix(:,7);
tau_2 = parameter_Matrix(:,8);
R_temp = 1-exp(-time./tau_1);
R_temp_2 = exp(-time./tau_2);

%% Initilization
spike_time = zeros(N_MU,1);
spike_time_mat = zeros(N_MU,length(time));
spike_train = zeros(N_MU,length(time));
force = zeros(N_MU,length(time));
Force = zeros(1,length(time));

% Module 2 parameteres
R = zeros(N_MU,length(time));
c = zeros(N_MU,1);
cf = zeros(N_MU,1);
A = zeros(N_MU,1);
c_mat = zeros(N_MU,length(time));
cf_mat = zeros(N_MU,length(time));
A_tilde_mat = zeros(N_MU,length(time));
A_mat = zeros(N_MU,length(time));

a_s = ones(N_MU,1)*0.96;
S_i = zeros(N_MU,1);
Y_i = zeros(N_MU,1);
S_mat = zeros(N_MU,length(time));
Y_mat = zeros(N_MU,length(time));

% Muscle length
Lce = 1;
Vce = 0;
FL = zeros(N_MU,1);
FV = zeros(N_MU,1);

DR_temp = zeros(N_MU,1);
DR_mat = zeros(N_MU,length(time));
%% Simulation
rng('shuffle')
for t = 1:length(time)    
    if t > 1
       %% Module 1       
        U_eff = synaptic_input(t);
        
        CV_ISI = 10+20*exp(-(U_eff*100-U_th*100)/2.5);
        CV_ISI = CV_ISI./100;      
        
        % for constant CoV of ISI
        % CV_ISI = ones(N_MU)*0.1;
              
       % compute discharge rate (DR_MU) 
        DR_MU = g_e.*(U_eff-U_th)+MDR;
        for m = 1:length(index_saturation)
            index = index_saturation(m);
            if U_eff <= U_th_t(index)
                DR_temp(index) = MDR(index) + lamda(index).*k_e(index)*(U_eff-U_th(index));
            else
                DR_temp(index) = PDR(index)-k_e(index)*(1-U_eff);
            end
        end
        DR_MU(index_saturation) = DR_temp(index_saturation);
        DR_MU(DR_MU<MDR) = 0;
        DR_MU(DR_MU>PDR) = PDR(DR_MU>PDR);
        
        DR_mat(:,t) = DR_MU;
        %% Sag & Yield (Song et al., 2008)
        f_eff = DR_MU./FR_half;
        S_i = sag_function(S_i,f_eff,a_s,Fs);
        S_i(1:index_slow) = 1;
        S_mat(:,t) = S_i;
        Y_i = yield_function(Y_i,Vce,Fs);
        Y_i(index_slow+1:end) = 1;
        Y_mat(:,t) = Y_i;
        
        % Convert activation into spike trains
        index_1 = i_MU(DR_MU >= MDR & DR_mat(:,t-1) == 0);
        index_2 = i_MU(DR_MU >= MDR & spike_time ==t);
        index = [index_1;index_2];
        
        % Generate spike trains
        index_1 = i_MU(DR_MU >= MDR & DR_mat(:,t-1) == 0); % find index of units that discharge for the first time
        index_2 = i_MU(DR_MU >= MDR & spike_time ==t); % find index of units whose spike time is at time = t
        index = [index_1;index_2];
        
        for j = 1:length(index) % loop through motor units whose firing rate is greater than minimum firing rate defined by the user
            n = index(j);
            spike_train_temp = zeros(1,length(t));
            if ~any(spike_train(n,:)) % when the motor unit fires for the first time
                spike_train(n,t) = 1; % add a spike to the vector
                spike_train_temp(t) = 1;
                
                % compute the spike time of the next spike
                mu = 1/DR_MU(n); % interspike interval                  
                spike_time_temp = (mu + mu*CV_ISI(n)*Z(n,t))*Fs; % add variabiltiy 
                if spike_time_temp <= 0.002*Fs
                    spike_time_temp = 0.002*Fs;
                end
                spike_time(n) = round(spike_time_temp) + t;
                
                % assign the value of R 
                temp = conv(spike_train_temp,R_temp_2(n,:).*R_temp(n,:));
                R(n,:) = R(n,:) + temp(1:length(time));
            else % when the motor unit have already fired at least once
                if spike_time(n) == t % when the motor unit fires
                    spike_train(n,t) = 1;
                    spike_train_temp(t) = 1;
                    
                    % compute the spike time of the next spike
                    mu = 1/DR_MU(n); % interspike interval                  
                    spike_time_temp = (mu + mu*CV_ISI(n)*Z(n,t))*Fs; % add variabiltiy 
                    if spike_time_temp <= 0.002*Fs
                        spike_time_temp = 0.002*Fs;
                    end
                    spike_time(n) = round(spike_time_temp) + t;
                    
                    % assign the value of R 
                    temp = conv(spike_train_temp,R_temp_2(n,:).*R_temp(n,:));
                    R(n,:) = R(n,:) + temp(1:length(time));
                elseif t > spike_time(n) + round(1/DR_MU(n)*Fs) % after the motor unit stops firing                   
                    spike_train(n,t) = 1;
                    spike_train_temp(t) = 1;
                    spike_time(n) = t;
                    
                    % compute the spike time of the next spike
                    mu = 1/DR_MU(n); % interspike interval
                    spike_time_temp = (mu + mu*CV_ISI(n)*Z(n,t))*Fs; % interspike interval
                    if spike_time_temp <= 0.002*Fs
                        spike_time_temp = 0.002*Fs;
                    end
                    spike_time(n) = round(spike_time_temp) + t;
                    
                    % assign the value of R 
                    temp = conv(spike_train_temp,R_temp_2(n,:).*R_temp(n,:));
                    R(n,:) = R(n,:) + temp(1:length(time));
                end
            end
        end
        
        %% Module 2: Convert spikes into activation
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
        %%
        f_i = A.*PTi.*FL.*FV;
        force(:,t) = f_i;
        
        Force(t) = sum(f_i);
    end
    spike_time_mat(:,t) = spike_time;
end

%%
if figOpt == 1
    figure(1)
    plot(time,Force)
    xlabel('Time (s)')
    ylabel('Force (N)')
    hold on
end

output.Force = Force;
output.spike_train = spike_train;

%% Convert spike trian into activation
    function [c,cf,A_tilde,A] = spike2activation(R,c,cf,A,parameter_Matrix,Lce,S_i,Y_i,Fs)
        
        S = parameter_Matrix(:,1);
        C = parameter_Matrix(:,2);
        k_1 = parameter_Matrix(:,3);
        k_2 = parameter_Matrix(:,4);
        k_3 = parameter_Matrix(:,5);
        k_4_i = parameter_Matrix(:,6);
        N = parameter_Matrix(:,9);
        K = parameter_Matrix(:,10);
        tau_3 = parameter_Matrix(:,11);
        gamma = parameter_Matrix(:,12);
        phi_1 = parameter_Matrix(:,13);
        phi_2 = parameter_Matrix(:,14);
        
        if Lce >= 1
            k_3 = (phi_1.*k_3)*(Lce-1) + k_3; 
            N = (-phi_1.*N)*(Lce-1) + N; 
            K = (-phi_1.*K)*(Lce-1) + K;
            gamma = (phi_1.*gamma)*(Lce-1) + gamma; 
        elseif Lce < 1
            k_3 = (phi_2.*k_3)*(Lce-1) + k_3; 
            N = (-phi_2.*N)*(Lce-1) + N;
            K = (-phi_2.*K)*(Lce-1) + K;
            gamma = (phi_2.*gamma)*(Lce-1) + gamma; 
        end
        
        %% Stage 1
        k_4 = k_4_i./(1+gamma.*A);
        c_dot = k_1.*(C-c-cf).*R - k_2.*c.*(S-C+c+cf)-(k_3.*c-k_4.*cf).*(1-cf);
        cf_dot = (1-cf).*(k_3.*c-k_4.*cf);
        c = c_dot/Fs + c;
        cf = cf_dot/Fs + cf;
        
        %% Stage 2
        if cf < 0
            cf_temp = 0;
        else
            cf_temp = cf.*S_i.*Y_i;
        end
        A_tilde = cf_temp.^N./(cf_temp.^N+K.^N);
        
        %% Stage 3
        % First-order dynamics to muscle activation, A
        A_dot = (A_tilde-A)./tau_3;
        A = A_dot./Fs + A;
        
    end

%% Sag
    function [S] = sag_function(S,f_eff,a_s,Fs)
        
        a_s(f_eff<0.1) = 20;
        
        T_s = 0.015;
        S_dot = (a_s - S)./T_s;
        S = S_dot/Fs + S;
        
    end

%% Yield
    function [Y] = yield_function(Y,V,Fs)
        c_y = 0.35;
        V_y = 0.1;
        T_y = 0.2;
        Y_dot = (1-c_y.*(1-exp(-abs(V)./V_y))-Y)./T_y;
        Y = Y_dot/Fs + Y;
        
    end

%% Force-length relationship for slow twitch
    function FL = FL_slow_function(L)
        
        beta = 2.3;
        omega = 1.12;
        rho = 1.62;
        
        FL = exp(-abs((L^beta - 1)/omega)^rho);
    end

%% Force-length relationship for fast twitch
    function FL = FL_fast_function(L)
       
        beta = 1.55;
        omega = 0.75;
        rho = 2.12;
        
        FL = exp(-abs((L^beta - 1)/omega)^rho);
    end

%% Concentric force-velocity relationship for slow twitch
    function FVcon = FVcon_slow_function(L,V)
       
        Vmax = -7.88;
        cv0 = 5.88;
        cv1 = 0;
        
        FVcon = (Vmax - V)/(Vmax + (cv0 + cv1*L)*V);
    end

%% Concentric force-velocity relationship for fast twitch
    function FVcon = FVcon_fast_function(L,V)
       
        Vmax = -9.15;
        cv0 = -5.7;
        cv1 = 9.18;
        
        FVcon = (Vmax - V)/(Vmax + (cv0 + cv1*L)*V);
    end

%% Eccentric force-velocity relationship for slow twitch
    function FVecc = FVecc_slow_function(L,V)
       
        av0 = -4.7;
        av1 = 8.41;
        av2 = -5.34;
        bv = 0.35;
        FVecc = (bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end

%% Eccentric force-velocity relationship for slow twitch
    function FVecc = FVecc_fast_function(L,V)
        
        av0 = -1.53;
        av1 = 0;
        av2 = 0;
        bv = 0.69;
        FVecc = (bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end


end