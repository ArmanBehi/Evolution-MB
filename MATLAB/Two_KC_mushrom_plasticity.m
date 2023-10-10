clear all;
close all;

%% Introducing inputs to system

% conditioned and unconditioned stimuli, over time
CS1 = [0,1,0,0,0,1,0,0,0];    % neutral odor 1

US1 = [0,1,0,0,0,1,0,0,0];    % reward
US2 = [0,0,0,1,0,0,0,1,0];    % punishment

% time increments and time axis

dt = 0.1;
N_i = round(numel(CS1)/dt);
t_i = (1:N_i)*dt; 

% replicate stimuli along time axis

n = round(1/dt);  %number of replications

CS1_i = Replicate(CS1,n);

US1_i = Replicate(US1,n);
US2_i = Replicate(US2,n);

%combination of inputs to Kenyon cells
S_i   = [CS1_i; US1_i; US2_i];   % combine stimulus activities for matrix multiplication

S_0   = zeros(3,1);  % initial values

%%  Kenyon cell activities

% allocate activity of K cells

K11_i = nan(1,N_i);  % odor 1, rewarded
K12_i = nan(1,N_i);  % odor 1, punished

K_0 = zeros(2,1);

K_i   = [K11_i; K12_i];   % combine K-cell activities for matrix multiplication


% obtain time-course of combined inputs to K-cells

alpha = 5;        % gain factor when US and CS coinicide

Kmax  = 5;
Khalf = 3;

K11_input_i = max( [CS1_i; US1_i; alpha*CS1_i.*US1_i], [], 1 );    % categorical values
K12_input_i = CS1_i;


K_input_i   = [K11_input_i; K12_input_i];   % combine for matrix multiplication

% steady-state value (this step is not really needed), but by changing
% K_half the influence of the US can be amplified or attenuated

K_ss_i = ActivationFunction( K_input_i, Kmax, Khalf );
%now we have the acticity of Kenyon cells based on the inputs we gave them

% obtain time-course of K-cell activities

tau_K = 0.5; %time-constant for kenyon cells

efactor_M = exp( - dt/tau_K ); %e^(-dt/tau)

for i = 1 : N_i
    
    K_i(:,i) = K_ss_i(:,i) + (K_0 - K_ss_i(:,i)) * efactor_M; %x(i+1) = F + (-F +x(i)).e(-dt/tau)
    
    K_0 = K_i(:,i); %redefining x(i)
end



%% Connection between KC's and MBON's
% obtain time-course of combined inputs to M-cells

beta1_1 = 1; % strength of excitatory projection from K-cells(plasticity shoul apply here)
beta1_2 = 1; %defining initial weight

%% MBON's activity
% allocate activity of M cells

M1_i = nan(1,N_i);  % approach
M2_i = nan(1,N_i);  % avoidance

M_0 = zeros(2,1);

M_i   = [M1_i; M2_i];   % combine M-cell activities for matrix multiplication
gamma = 1;  % strength of inhibitory projection from M-cells

           %inhibition, C1 attracts, C2 attracts
M_matrix = [beta1_1 , 0 , 0, gamma; 0 , beta1_2 , 0 , gamma];
           %[kc1_1 , kc1_2  , M1 , M2]
           
tau_M = 0.5; %time-constant for MBON's
tau_W = 1;
tau_theta = 1;
efactor_M = exp( - dt/tau_M );  %e^(-dt/tau)
efactor_theta = exp( - dt/tau_theta);  %e^(-dt/tau)

Mmax = 5;
Mhalf = 3;

K_0 = zeros(2,1);
M_0 = zeros(2,1);
theta = [1;1];
W=[1 ; 1]; %initial weight

for i = 1 : N_i
    
    M_input_i = M_matrix * [K_0',M_0']';
    z(:,i)=M_input_i;
    
    M_ss_i = ActivationFunction( M_input_i, Mmax, Mhalf );
    
    M_i(:,i) = M_ss_i + (M_0 - M_ss_i) * efactor_M;
    
    M_0 = M_i(:,i);
    
    K_0 = K_i(:,i);
    
    %imply plasticity
    
    %     delta_theta = ((M_0.^2) - theta)./tau_theta;
    %     theta = theta + delta_theta;
    %     TT(i,:)=theta;
    %     delta_W = ((M_0*K_0')*(M_0 - theta))./tau_W;
    %     DD(i,:) = delta_W;
    %     W = W + delta_W;
    %     QZ(:,i)=W;
    %     W(1,1) = beta1_1;
    %     W(2,1) = beta1_2;
    
    theta(:,i+1) = (M_0).^2 - ((M_0).^2 + theta(:,i))*(efactor_theta);
    theta(:,i+1) = ActivationFunction( theta(:,i+1), Mmax, Mhalf );
    delta_W(:,i)= ((M_0*K_0')*(M_0 - theta(:,i+1)))./tau_W;
    beta1_1 = sum(delta_W(1,1:i+1));
    beta1_2 = sum(delta_W(2,1:i+1));
    
    M_matrix = [beta1_1 , 0 , 0, gamma; 0 , beta1_2 , 0 , gamma];
%     xx(:,:,i) = M_matrix;
end

figure;

subplot(3,1,1);
hold on;
plot(t_i, K_input_i, 'LineWidth', 2);
hold off;
legend('( KC1_1) odour1 + reward','(KC1_2) odour1 + punishment');
title('Inputs to K-cells');

subplot(3,1,2);
hold on;
plot(t_i, K_i, 'LineWidth', 2 );
hold off;
legend('K1_1','K1_2');
title('Activity of K-cells');

subplot(3,1,3);
hold on;
plot(t_i, M_i, 'LineWidth', 2 );
hold off;
legend('Attractor MBON','Aversion MBON');
title('Activity of M-cells');
        
return;