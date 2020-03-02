%==========================================================================
% MU_population_model.m
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

function [output] = MU_population_model(Fs,time,synaptic_drive,modelParameter,figOpt)
%% Muscle architectural parameters
L0 = modelParameter.optimalLength; % optimal muscle length [cm]
F0 = modelParameter.F0; % maximal force

L0T = modelParameter.L0T; % optima tendon length [cm]
alpha = modelParameter.pennationAngle; % pennation angle [radians]
Lmt =modelParameter.Lmt; % intial musculotendon length [cm]
L_ce = modelParameter.L_ce; % normalized muscle length [L0]
L_se = modelParameter.L_se; % normalized tendon length [L0T]
Lmax = modelParameter.Lmax; % maximal muscle length [cm]

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
Z = randn(N_MU,length(time));
Z(Z>3.9) = 3.9;
Z(Z<-3.9) = -3.9;

%% Parameter initilization
DR_temp = zeros(N_MU,1);
DR_mat = zeros(N_MU,length(time));

spike_time = zeros(N_MU,1);
spike_train = zeros(N_MU,length(time));
force = zeros(N_MU,length(time));
F_se = zeros(1,length(time));

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

FL = zeros(N_MU,1);
FV = zeros(N_MU,1);

MuscleVelocity = zeros(1,length(time));
MuscleLength = zeros(1,length(time));
MuscleLength(1) = L_ce*L0/100;

V_ce = 0;

%% Simulation
rng('shuffle')
h = 1/Fs; % time step
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
        Y_i = yield_function(Y_i,V_ce,Fs);
        Y_i(index_slow+1:end) = 1;
        Y_mat(:,t) = Y_i;
        
        %% Convert activation into spike trains
        index_1 = i_MU(DR_MU >= MDR & DR_mat(:,t-1) == 0);
        index_2 = i_MU(DR_MU >= MDR & spike_time == t);
        index = [index_1;index_2];
        
        for j = 1:length(index) % loop through motor units whose firing rate is greater than minimum discharge rate
            n = index(j);
            spike_train_temp = zeros(1,length(t));
            if ~any(spike_train(n,:)) % when the motor unit fires at the first time
                spike_train(n,t) = 1; % add a spike to the vector
                spike_train_temp(t) = 1;
                mu = 1/DR_MU(n); % interspike interval
                spike_time_temp = (mu + mu*cv_MU*Z(n,t))*Fs; % compute next spike time
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
                    % update mean firing rate of the motor unit given the current value of input
                    mu = 1/DR_MU(n); % interspike interval
                    spike_time_temp = (mu + mu*cv_MU*Z(n,t))*Fs; % compute next spike time
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
                    spike_time_temp = (mu + mu*cv_MU*Z(n,t))*Fs; % compute next spike time
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
    end
    %% Convert spikes into activation (Eq. 10-19)
    [c,cf,A_tilde,A] = spike2activation(R(:,t),c,cf,A,parameter_Matrix,L_ce,S_i,Y_i,Fs);
    
    c_mat(:,t) = c;
    cf_mat(:,t) = cf;
    A_tilde_mat(:,t) = A_tilde;
    A_mat(:,t) = A;
    
    %% Muscle contraction dynamics (Eq. 21-22)
    [force(:,t),F_se(t)] = contraction_dynamics_v2(A,L_se,L_ce,V_ce,FL,FV,index_slow,Lmax,PTi,F0);
    
    k_0_de = h*MuscleVelocity(t);
    l_0_de = h*contraction_dynamics(A,L_se,L_ce,V_ce,FL,FV,modelParameter,index_slow,Lmax,PTi,F0);
    k_1_de = h*(MuscleVelocity(t)+l_0_de/2);
    l_1_de = h*contraction_dynamics(A,(Lmt - (L_ce+k_0_de/L0)*L0*cos(alpha))/L0T,L_ce+k_0_de/L0,V_ce+l_0_de/L0,FL,FV,modelParameter,index_slow,Lmax,PTi,F0);
    k_2_de = h*(MuscleVelocity(t)+l_1_de/2);
    l_2_de = h*contraction_dynamics(A,(Lmt - (L_ce+k_1_de/L0)*L0*cos(alpha))/L0T,L_ce+k_1_de/L0,V_ce+l_1_de/L0,FL,FV,modelParameter,index_slow,Lmax,PTi,F0);
    k_3_de = h*(MuscleVelocity(t)+l_2_de);
    l_3_de = h*contraction_dynamics(A,(Lmt - (L_ce+k_2_de/L0)*L0*cos(alpha))/L0T,L_ce+k_2_de/L0,V_ce+l_2_de/L0,FL,FV,modelParameter,index_slow,Lmax,PTi,F0);
    MuscleLength(t+1) = MuscleLength(t) + 1/6*(k_0_de+2*k_1_de+2*k_2_de+k_3_de);
    MuscleVelocity(t+1) = MuscleVelocity(t) + 1/6*(l_0_de+2*l_1_de+2*l_2_de+l_3_de);
    
    % normalize each variable to optimal muscle length or tendon length
    V_ce = MuscleVelocity(t+1)/(L0/100);
    L_ce = MuscleLength(t+1)/(L0/100);
    L_se = (Lmt - L_ce*L0*cos(alpha))/L0T;
end

%% Plotting figures
if figOpt == 1
    figure(1)
    plot(time,F_se)
    xlabel('Time (s)')
    ylabel('Force (N)')
    hold on
end

%% Save data in the data structure, output
output.spike_train = spike_train; % spike trains
output.ForceTendon = F_se; % tendon force
output.force = force; % individual motor unit forces

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

    %% Sag (Eq.16-17)
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

    %% Eccentric force-velocity relationship for slow twitch (Eq.25)
    function FVecc = FVecc_fast_function(L,V)
        av0 = -1.53;
        av1 = 0;
        av2 = 0;
        bv = 0.69;
        FVecc = (bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end

    %% Force-length relationship for passive element 1 (Song et al. 2008)
    function Fpe1 = Fpe1_function(L,V)
     
        c1_pe1 = 23;
        k1_pe1 = 0.046;
        Lr1_pe1 = 1.17;
        eta = 0.01;
        
        Fpe1 = c1_pe1 * k1_pe1 * log(exp((L - Lr1_pe1)/k1_pe1)+1) + eta*V;
        
    end

    %% Force-length relationship for passive element 2  (Song et al. 2008)
    function Fpe2 = Fpe2_function(L)
      
        c2_pe2 = -0.02;
        k2_pe2 = -21;
        Lr2_pe2 = 0.70;
        
        Fpe2 = c2_pe2*exp((k2_pe2*(L-Lr2_pe2))-1);
        
    end

    %% Force-length relationship for series-elastic element  (Song et al. 2008)
    function Fse = Fse_function(LT)
      
        cT_se = 27.8; 
        kT_se = 0.0047;
        LrT_se = 0.964;
        
        Fse = cT_se * kT_se * log(exp((LT - LrT_se)/kT_se)+1);
        
    end
    %% Muscle contraction dynamics (Eq. 21)
    function ddx = contraction_dynamics(A,L_s,L_m,L_m_dot,FL_vec,FV_vec,modelParameter,index_slow,Lmax,PT,F0)
        %% Force-length and force-velocity
        FL_vec(1:index_slow) = FL_slow_function(L_m);
        FL_vec(index_slow+1:end) = FL_fast_function(L_m);
        
        if L_m_dot > 0
            FV_vec(1:index_slow) = FVecc_slow_function(L_m,L_m_dot);
            FV_vec(index_slow+1:end) = FVecc_fast_function(L_m,L_m_dot);
        else
            FV_vec(1:index_slow) = FVcon_slow_function(L_m,L_m_dot);
            FV_vec(index_slow+1:end) = FVcon_fast_function(L_m,L_m_dot);
        end
        
        %% Passive element 1
        F_pe1 = Fpe1_function(L_m/Lmax,L_m_dot);
        
        %% Passive element 2
        F_pe2 = Fpe2_function(L_m);
        if F_pe2 > 0
            F_pe2 = 0;
        end
        
        f_i = A.*PT.*(FL_vec.*FV_vec+F_pe2);
        
        F_m_temp = sum(f_i);
        F_m = F_m_temp + F_pe1*F0;
        
        F_t = Fse_function(L_s) * F0;
        
        M = modelParameter.mass;
        rho = modelParameter.pennationAngle;
        
        ddx = (F_t*cos(rho) - F_m*(cos(rho)).^2)/(M) ...
            + (L_m_dot).^2*tan(rho).^2/(L_m);
    end
    
    %% Muscle contraction dynamics (Eq. 21)
    function [f_i,F_t] = contraction_dynamics_v2(A,L_s,L_m,L_m_dot,FL_vec,FV_vec,index_slow,Lmax,PT,F0)
        %% Force-length and force-velocity
        FL_vec(1:index_slow) = FL_slow_function(L_m);
        FL_vec(index_slow+1:end) = FL_fast_function(L_m);
        
        if L_m_dot > 0
            FV_vec(1:index_slow) = FVecc_slow_function(L_m,L_m_dot);
            FV_vec(index_slow+1:end) = FVecc_fast_function(L_m,L_m_dot);
        else
            FV_vec(1:index_slow) = FVcon_slow_function(L_m,L_m_dot);
            FV_vec(index_slow+1:end) = FVcon_fast_function(L_m,L_m_dot);
        end
        
        %% Passive element 1
        F_pe1 = Fpe1_function(L_m/Lmax,L_m_dot);
        
        %% Passive element 2
        F_pe2 = Fpe2_function(L_m);
        if F_pe2 > 0
            F_pe2 = 0;
        end
        
        f_i = A.*PT.*(FL_vec.*FV_vec+F_pe2);
        
        F_t = Fse_function(L_s) * F0;
        
    end
end