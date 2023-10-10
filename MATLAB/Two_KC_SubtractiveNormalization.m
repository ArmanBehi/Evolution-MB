
clear all;
close all;

%% conditioned and unconditioned stimuli, over time
CS1 = [0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0];    % neutral odor 1
% CS1 = repmat(CS1,1,10);

US1 = [0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0];    % reward
US2 = [0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0];    % punishment
% US1 = repmat(US1,1,10);
% US2 = repmat(US2,1,10);

alpha = 5;        % gain factor when US and CS coinicide

Kmax  = 5;        % maximal K-cell activity
Khalf = 3;        % non-linearity of activation function

K_0   = 0;        % initial K-cell activity        

tau_K = 0.5;      % time-constant of K-cell activity

Mmax = 5;         % maximal M-cell activity
Mhalf = 3;        % non-linearity of activation function

M_0  = 0;         % initial M-cell activity

tau_M = 0.5;      % time-constant of M-cell activity

beta = 0.8;       % initial strength of excitatory projection from K-cells (plastic synapses)

gamma = -0.5;     % strength of inhibitory projection from M-cells (fixed synapses)

tau_Theta = 10;   % time-constant of plasticity threshold
tau_W     = 10;   % time-constant of plasticity strength

Theta_0 = 0.5;    % initial value of plasticity threshold

% time increments and time axis

dt = 0.1;
N_i = round(numel(CS1)/dt);
t_i = (1:N_i)*dt; 

% replicate stimuli along time axis

n = round(1/dt);  %number of replications

CS1_i = Replicate(CS1,n);
US1_i = Replicate(US1,n);
US2_i = Replicate(US2,n);
Half_way = round(N_i/2);
% US1_i(:,Half_way:end) =0;
% US2_i(:,Half_way:end) =0;
% quarter_way = round(N_i/5);
% US1_i(:,2*quarter_way:4*quarter_way) =0;
% US2_i(:,2*quarter_way:4*quarter_way) =0;

% allocate activity of K cells

K_i = nan(2,N_i);  % odor 1, rewarded
                   % odor 1, punished
                  

% allocate activity of M cells

M_i = nan(2,N_i);  % approach% avoidance

%% obtain time-course of combined inputs to K-cells

K11_input_i = max( [CS1_i; US1_i; alpha*CS1_i.*US1_i], [], 1 );    % categorical values
K12_input_i = max( [CS1_i; US2_i; alpha*CS1_i.*US2_i], [], 1 );

K_input_i   = [K11_input_i; K12_input_i];   % combine for matrix multiplication

% steady-state value (this step is not really needed), but by changing
% Khalf the influence of the US can be amplified or attenuated

K_ss_i = ActivationFunction( K_input_i, Kmax, Khalf );

% obtain time-course of K-cell activities

efactor = exp( - dt/tau_K );

K_ini = K_0*ones(2,1);

for i = 1 : N_i
    
    K_i(:,i) = K_ss_i(:,i) + (K_ini - K_ss_i(:,i)) * efactor;
    
    K_ini = K_i(:,i);
end

% 
% figure;
% subplot(4,1,1);
% hold on;
% plot(t_i, CS1_i, 'LineWidth', 2);
% hold off;
% title('CS');
% 
% subplot(4,1,2);
% hold on;
% plot(t_i, US1_i, 'LineWidth', 2);
% hold off;
% title('US1');
% 
% subplot(4,1,3);
% hold on;
% plot(t_i, US2_i, 'LineWidth', 2);
% hold off;
% title('US2');
% 
% subplot(4,1,4);
% hold on;
% plot(t_i, K_i, 'LineWidth', 2 );
% hold off;
% 
% title('Activity of K-cells');


%% obtain time-course of combined inputs to M-cells


           %inhibition, C1 attracts, C2 attracts
M_matrix = [0 gamma beta 0;gamma 0 beta 0];
        
efactor = exp( - dt/tau_M );

K_ini = K_0*ones(2,1);
M_ini = M_0*ones(2,1);

for i = 1 : N_i
    
    M_input_i = M_matrix * [M_ini', K_ini']';
    
    M_ss_i = ActivationFunction( M_input_i, Mmax, Mhalf );
    
    M_i(:,i) = M_ss_i + (M_ini - M_ss_i) * efactor;
    
    M_ini = M_i(:,i);
    
    K_ini = K_i(:,i);
end

% figure;
% subplot(2,1,1);
% hold on;
% plot(t_i, K_i, 'LineWidth', 2);
% hold off;
% 
% title('Activity of K-cells');
% 
% subplot(2,1,2);
% hold on;
% plot(t_i, M_i, 'LineWidth', 2 );
% hold off;
% 
% title('Activity of M-cells');
%   

%% obtain time course of synaptic strength, but without influencing M-cell activity

dW_i    = nan(2,2,N_i);            % allocate synpatic weight change and synaptic threshold
Theta_i = nan(2,N_i);            % K11 to M1


%odor1/rewarded, odor1/punished, odor2/rewarded, odor2/punished
W_ini = [beta, 0; 0, beta];
   
W_connectivity = [W_ini~=0];
   
efactor_Theta = exp( -dt/tau_Theta );

Theta_ini = Theta_0 * ones(2,1);

alpha = 1; %alpha term for Oja rule

for i = 1 : N_i
    
%     Theta_ss_i = M_i(:,i).^2;   % instantaneous steady state
%     
%     Theta_i(:,i) = Theta_ss_i + (Theta_ini - Theta_ss_i) * efactor_Theta;
%     
%     Theta_ini = Theta_i(:,i); %theta(threshhold) at time t
    
    U_i = repmat( K_i(:,i), 1, 1 ); % presynaptic activity
    
    V_i = repmat(M_i(:,i),1,1);      % postsynaptic activity
    
%     T_i = repmat(Theta_i(:,i),1,1);  % threshold
    
%     dW_i(:,:,i) =  V_i .* U_i -(W_connectivity .* alpha.*(V_i.^2))/ tau_W;   % Oja Rule
     dW_i(:,:,i) =  2*(V_i.^2).*(1-(W_connectivity.^2).*(alpha))/ tau_W;   % Oja Rule2
end

dW_i = reshape( dW_i, 4, N_i );

W_i = W_ini(:) + cumsum(dW_i,2) * dt;

% figure;
% subplot(2,1,1);
% hold on;
% plot(t_i, dW_i, 'LineWidth', 2);
% hold off;
% 
% title('Weight change');
% 
% subplot(2,1,2);
% hold on;
% plot(t_i, W_i, 'LineWidth', 2);
% hold off;
% 
% title('Weights');

%% coupled dynamics: activity changes synaptic strength, and synaptic strenght influences activity

efactor_M = exp( - dt/tau_M );

efactor_Theta = exp( -dt/tau_Theta );

K_ini = K_0*ones(2,1);
M_ini = M_0*ones(2,1);

Theta_ini = Theta_0*ones(2,1);

            % mutual inhibition
M_matrix = [0, gamma; ...
             gamma, 0];
 
            % plastic feeforward input
W_matrix = [beta, 0, ;0, beta];
   
W_connectivity = [W_matrix~=0];

dW_i   = nan(2,2,N_i);

W_matrix_i   = nan(2,2,N_i);


for i = 1 : N_i
    
    M_input_i = W_matrix * K_i(:,i);           % use stored values of K-activity
    
    M_input_i = M_input_i + M_matrix * M_ini;    % use computed values of M-activity
    
    M_ss_i = ActivationFunction( M_input_i, Mmax, Mhalf );
    
    M_i(:,i) = M_ss_i + (M_ini - M_ss_i) * efactor_M;
    
    M_ini = M_i(:,i);
    
    
    U_i = repmat( K_i(:,i)', 1, 1 ); % presynaptic activity
    
    V_i = repmat(M_i(:,i),1,1);      % postsynaptic activity
    
    
    
%     dW_i(:,:,i) =  V_i .* U_i -(W_connectivity .* alpha.*(V_i.^2))/ tau_W;   % Oja Rule
    dW_i(:,:,i) =  2*(V_i.^2).*(1-(W_connectivity.^2).*(alpha))/ tau_W;   % Oja Rule2
    
    W_matrix = W_matrix + dW_i(:,:,i);
    
    W_matrix_i(:,:,i) = W_matrix;
    
end


figure;
subplot(4,1,1);
hold on;
plot(t_i, K_i, 'LineWidth', 2 );
hold off;
title('Activity of K-cells');

dW_i = reshape( dW_i, 4, N_i );
subplot(4,1,2);
hold on;
plot(t_i, dW_i, 'LineWidth', 2);
hold off;
title('Weight change');

subplot(4,1,3);
hold on;
plot(t_i, M_i, 'LineWidth', 2);
hold off;
title('Activity of M-cells');

W_matrix_i = reshape( W_matrix_i, 4, N_i );

subplot(4,1,4);
hold on;
plot(t_i, W_matrix_i, 'LineWidth', 2);
hold off;
title('Weights with variable activity');
