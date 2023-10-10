clear all;
close all;

dt = 0.1;
N_i = round(numel(10)/dt);
t_i = (1:N_i)*dt; 

Mmax = 1.5;         % maximal M-cell activity
Mhalf = 0.75; 




%% M cell Activity

%assigning parameter
alpha =  1;
beta  =  2;
w11 = 1;
w12 = 1;
w21 = 1;
w22 = 1;
tau_M = 0.1;

M_matrix = [alpha, -beta; -beta, alpha]; %M cells level connection

M_input = [w11 0 w21 0 ; 0 w12 0 w22]; %matrix of K cells input to M cells
        
efactor = exp( - dt/tau_M );

k11 = linspace(0,1,N_i);
k12 = linspace(0,1,N_i);
k21 = linspace(0,1,N_i);
k22 = linspace(0,1,N_i);

K_ini = [k11;k12;k21;k22];
y= K_ini(1,:) + K_ini(3,:);
x= K_ini(2,:) + K_ini(4,:);

[X,Y] = meshgrid(x,y);

% M_ini = zeros(2,1);
M_ini = [0.8;0.1];

M_out = zeros(2,10,N_i);

for i = 1:N_i
    
    for j=1:N_i
        
        %Initial value for M's
        
        
        M_mutual = M_matrix * [M_ini];
        K_to_M = M_input * [K_ini(:,j)];
        
        M_input_i  = M_mutual + K_to_M;
        
        M_ss_i = ActivationFunction( M_input_i, Mmax, Mhalf );
        
        M_non(:,j) = M_ss_i + (M_ini - M_ss_i) * efactor;
        
        M_ini = M_non(:,j);
    end
    M_out(:,:,i) = M_non;
end

M1(:,:) = M_out(1,:,:);

M2(:,:) = M_out(2,:,:);
% Z = M_non(1,:);
figure;
subplot(2,1,1)
contourf(X,Y,M1,'LineWidth',2)
title('Activitiy of M1')
xlabel('K12.W12 + k22.W22');
ylabel('K11.W11 + k21.W21');
axis equal;


subplot(2,1,2)
contour(X,Y,M2,'LineWidth',2,'ShowText','on')
title('Activitiy of M2')
xlabel('K12.W12 + k22.W22');
ylabel('K11.W11 + k21.W21');
axis equal;
annotation('textbox', [0.1, 0.85, 0.1, 0.1], 'String', "alpha(self-excitation) = " + alpha +", beta(mutual-inhibition) ="+beta)
