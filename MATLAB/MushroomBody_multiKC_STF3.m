
clear all;
close all;

%% conditioned and unconditioned stimuli, over time
CS_1 = [0,1,0,0,0,1,0,0,0,1,0,0,0,1,0];    % neutral odor 1
CS_2 = [0,0,0,1,0,0,0,1,0,0,0,1,0,0,0];    % neutral odor 2

US_1 = [0,1,0,0,0,1,0,0,0,1,0,0,0,1,0];    % reward
US_2 = [0,0,0,1,0,0,0,1,0,0,0,1,0,0,0];    % punishment

US11=[US_1,US_1,zeros(1,length(US_1)),US_1]; %rewarding US to K11
US12=[zeros(1,length(US_1)),zeros(1,length(US_1)),zeros(1,length(US_1)),zeros(1,length(US_1))]; %rewarding US to K21
US21=[zeros(1,length(US_1)),zeros(1,length(US_1)),zeros(1,length(US_1)),zeros(1,length(US_1))]; %punishing US to K12
US22=[zeros(1,length(US_1)),zeros(1,length(US_1)),zeros(1,length(US_1)),zeros(1,length(US_1))]; %punishing US to K22

CS1=[US_1,US_1,zeros(1,length(US_1)),US_1]; 
CS2=[zeros(1,length(US_1)),zeros(1,length(US_1)),US_2,zeros(1,length(US_1))];

replicate_number = 2;

CS1 = repmat(CS1,1,replicate_number);
CS2 = repmat(CS2,1,replicate_number);

US11 = repmat(US11,1,replicate_number);
US12 = repmat(US12,1,replicate_number);
US21 = repmat(US21,1,replicate_number);
US22 = repmat(US22,1,replicate_number);

alpha = 2;        % gain factor when US and CS coinicide

Kmax  = 5;        % maximal K-cell activity
Khalf = 3;        % non-linearity of activation function

K_0   = 0;        % initial K-cell activity        

tau_K = 0.5;      % time-constant of K-cell activity

Mmax = 5;         % maximal M-cell activity
Mhalf = 3;        % non-linearity of activation function

M_0  = 0;         % initial M-cell activity

tau_M = 0.5;      % time-constant of M-cell activity

beta = 0.5;       % initial strength of excitatory projection from K-cells (plastic synapses)

gamma = -0.8;     % strength of inhibitory projection from M-cells (fixed synapses)

tau_Theta = 10.0;   % time-constant of plasticity threshold
tau_W     = 10;   % time-constant of plasticity strength

Theta_0 = 0.5;    % initial value of plasticity threshold

% time increments and time axis

dt = 1;
N_i = round(numel(CS1)/dt);
t_i = (1:N_i)*dt; 

% replicate stimuli along time axis

n = round(1/dt);  %number of replications

CS1_i = Replicate(CS1,n);
CS2_i = Replicate(CS2,n);
US1_1 = Replicate(US11,n); %US to K11
US1_2 = Replicate(US12,n); %US to K12
US2_1 = Replicate(US21,n); %US to K21
US2_2 = Replicate(US22,n); %Us to K22

CS1_i(1,90:end) = CS1_i(1,1:31);
US1_1(1,90:end) = 0;
% allocate activity of K cellss

K_i = nan(4,N_i);  % odor 1, rewarded
                   % odor 1, punished
                   % odor 2, rewarded
                   % oder 2, punished

% allocate activity of M cells

M_i = nan(2,N_i);  % approach% avoidance

%% obtain time-course of combined inputs to K-cells

K11_input_i = max( [CS1_i; US1_1; alpha*CS1_i.*US1_1], [], 1 );    % categorical values
K12_input_i = max( [CS1_i; US1_2; alpha*CS1_i.*US1_2], [], 1 );
K21_input_i = max( [CS2_i; US2_1; alpha*CS2_i.*US2_1], [], 1 );
K22_input_i = max( [CS2_i; US2_2; alpha*CS2_i.*US2_2], [], 1 );

K_input_i   = [K11_input_i; K12_input_i; K21_input_i; K22_input_i];   % combine for matrix multiplication

% steady-state value (this step is not really needed), but by changing
% Khalf the influence of the US can be amplified or attenuated

K_ss_i = ActivationFunction( K_input_i, Kmax, Khalf );

% obtain time-course of K-cell activities

efactor = exp( - dt/tau_K );

K_ini = K_0*ones(4,1);

for i = 1 : N_i
    
    K_i(:,i) = K_ss_i(:,i) + (K_ini - K_ss_i(:,i)) * efactor;
    
    K_ini = K_i(:,i);
    
end

% figure;
% subplot(2,1,1);
% hold on;
% plot(t_i, K_input_i, 'LineWidth', 2);
% hold off;
% 
% title('Inputs to K-cells');
% 
% subplot(2,1,2);
% hold on;
% plot(t_i, K_i, 'LineWidth', 2 );
% hold off;
% 
% title('Activity of K-cells');


%% obtain time-course of combined inputs to M-cells

% 
%            %inhibition, C1 attracts, C2 attracts
% M_matrix = [0, gamma, beta, 0, beta, 0; ...
%            % inhibition, C1 avoids, C2 avoids
%             gamma, 0, 0, beta, 0, beta];
%         
% efactor = exp( - dt/tau_M );
% 
% K_ini = K_0*ones(4,1);
% M_ini = M_0*ones(2,1);
% 
% for i = 1 : N_i
%     
%     M_input_i = M_matrix * [M_ini', K_ini']';
%     
%     M_ss_i = ActivationFunction( M_input_i, Mmax, Mhalf );
%     
%     M_i(:,i) = M_ss_i + (M_ini - M_ss_i) * efactor;
%     
%     M_ini = M_i(:,i);
%     
%     K_ini = K_i(:,i);
% end

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
% 
% dW_i    = nan(2,4,N_i);            % allocate synpatic weight change and synaptic threshold
% Theta_i = nan(2,N_i);            % K11 to M1
% 
% 
% %odor1/rewarded, odor1/punished, odor2/rewarded, odor2/punished
% W_ini = [beta, 0, beta, 0; ...
%          0, beta, 0, beta];
%    
% W_connectivity = [W_ini~=0];
%    
% Pmax = 5;         % maximal M-cell activity
% Phalf = 3;        % non-linearity of activation function
% P_eff_1 = 0;
% P_eff_2 = 0;
% P_eff_3 = 0;
% P_eff_4 = 0;
% tau_P = 0.9; 
% P_efactor = exp(-1/tau_P );
% f_D = 0.5;
% 
% for i = 1:N_i
%     %
%     if K_i(1,i)>1
%         P_eff_1 = P_eff_1 + f_D*(1-P_eff_1);
%         beta11_1(:,i) = P_eff_1;
%     else
%        P_0 = ActivationFunction( P_eff_1, Pmax, Phalf );
%         P_eff_1 = P_0 + (P_eff_1 - P_0) * P_efactor;
%         beta11_1(:,i) = P_eff_2;
%     end
%     %
%      if K_i(2,i)>1
%         P_eff_2 = P_eff_2 + f_D*(1-P_eff_2);
%         beta12_2(:,i) = P_eff_2;
%     else
%         P_0 = ActivationFunction( P_eff_2, Pmax, Phalf );
%         P_eff_2 = P_0 + (P_eff_2 - P_0) * P_efactor;
%         beta12_2(:,i) = P_eff_2;
%      end
%     %
%       if K_i(3,i)>1
%         P_eff_3 = P_eff_3 + f_D*(1-P_eff_3);
%         beta21_1(:,i) = P_eff_3;
%     else
%         P_0 = ActivationFunction( P_eff_3, Pmax, Phalf );
%         P_eff_3 = P_0 + (P_eff_3 - P_0) * P_efactor;
%         beta21_1(:,i) = P_eff_3;
%       end
%     %
%       if K_i(4,i)>1
%         P_eff_4 = P_eff_1 + f_D*(1-P_eff_4);
%         beta22_2(:,i) = P_eff_4;
%     else
%         P_0 = ActivationFunction( P_eff_4, Pmax, Phalf );
%         P_eff_4 = P_0 + (P_eff_4 - P_0) * P_efactor;
%         beta22_2(:,i) = P_eff_4;
%       end
% end

% figure;
% subplot(4,1,1);
% hold on;
% plot(t_i, beta11_1, 'LineWidth', 2);
% hold off;
% 
% title('K11_M1 (Rewarding)');
% 
% subplot(4,1,2);
% hold on;
% plot(t_i, beta12_2, 'LineWidth', 2);
% hold off;
% 
% title('K12-M2 (Punishing)');
% 
% subplot(4,1,3);
% hold on;
% plot(t_i, beta21_1, 'LineWidth', 2);
% hold off;
% 
% title('K21-M2 (Rewarding)');
% 
% subplot(4,1,4);
% hold on;
% plot(t_i, beta22_2, 'LineWidth', 2);
% hold off;
% 
% title('K22-M2 (Punishing)');

%% coupled dynamics: activity changes synaptic strength, and synaptic strenght influences activity

efactor_M = exp( - 1/tau_M );

efactor_Theta = exp( -dt/tau_Theta );

K_ini = K_0*ones(4,1);
M_ini = M_0*ones(2,1);

Theta_ini = Theta_0*ones(2,1);

            % mutual inhibition
M_matrix = [0, gamma; ...
             gamma, 0];
 
            % plastic feeforward input
W_matrix = [beta, 0, beta, 0; ...
            0, beta, 0, beta];
   
W_connectivity = [W_matrix~=0];

dW_i   = nan(2,4,N_i);

W_matrix_i   = nan(2,4,N_i);

% Values for P_eff
Pmax = 5;         % maximal M-cell activity
Phalf = 3;        % non-linearity of activation function
P_eff_1 = 0;
P_eff_2 = 0;
P_eff_3 = 0;
P_eff_4 = 0;
tau_P = 10; 
P_efactor = exp(-dt/tau_P );
f_D = 0.5;

%initial Values for P_eff
P_eff_1 = 0.3;
P_eff_2 = 0.3;
P_eff_3 = 0.3;
P_eff_4 = 0.3;

P_0 = 0.3;

for i = 1 : length(CS1)
    
%     M_input_i = W_matrix * K_i(:,i);           % use stored values of K-activity
    M_input_i = W_matrix * K_input_i(:,i);            %In this case efactor_M = exp( - 1/tau_M )
    
    M_input_i = M_input_i + M_matrix * M_ini;    % use computed values of M-activity
    
    M_ss_i = ActivationFunction( M_input_i, Mmax, Mhalf );
    
    M_i(:,i) = M_ss_i + (M_ini - M_ss_i) * efactor_M;
    
    M_ini = M_i(:,i);
    
    if K_i(1,i)>1
        
        P_eff_1 = P_eff_1 + f_D*(1-P_eff_1);
        beta11_1(:,i) = P_eff_1;
    else
%         P_0 = ActivationFunction( P_eff_1, Pmax, Phalf );
        P_eff_1 = P_0 + (P_eff_1 - P_0) * P_efactor;
        beta11_1(:,i) = P_eff_1;
    end
    %
    if K_i(2,i)>1
        
        P_eff_2 = P_eff_2 + f_D*(1-P_eff_2);
        beta12_2(:,i) = P_eff_2;
    else
%         P_0 = ActivationFunction( P_eff_2, Pmax, Phalf );
        P_eff_2 = P_0 + (P_eff_2 - P_0) * P_efactor;
        beta12_2(:,i) = P_eff_2;
    end
    %
    if K_i(3,i)>1
        
        P_eff_3 = P_eff_3 + f_D*(1-P_eff_3);
        beta21_1(:,i) = P_eff_3;
    else
%         P_0 = ActivationFunction( P_eff_3, Pmax, Phalf );
        P_eff_3 = P_0 + (P_eff_3 - P_0) * P_efactor;
        beta21_1(:,i) = P_eff_3;
    end
    %
    if K_i(4,i)>1
        
        P_eff_4 = P_eff_1 + f_D*(1-P_eff_4);
        beta22_2(:,i) = P_eff_4;
    else
%         P_0 = ActivationFunction( P_eff_4, Pmax, Phalf );
        P_eff_4 = P_0 + (P_eff_4 - P_0) * P_efactor;
        beta22_2(:,i) = P_eff_4;
    end
    
    W_matrix = [beta11_1(:,i), 0, beta21_1(:,i), 0; 0, beta12_2(:,i), 0, beta22_2(:,i)];
    
end




