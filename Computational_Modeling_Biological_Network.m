clc;
clear;
close all;

%"Computational Modeling of a Biological Network Implementing the
%Rescorla-Wagner Learning Rule: A Simulation Study"
% This Script is written by Arman Behrad under supervision of
% Prof.Dr.Juchen Braun
% the sections consistenet with the project manuscript

%%  Introducing Constant CS*US
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

%plotting
figure()


subplot(3,1,1)
plot(CS1_i , 'b','LineWidth',4)
title("CS - Odour")
ylabel("Intensity")

% Get the current axes handle
ax = gca;
% Set the font size for all text elements
fontSizeValue = 14; % Adjust the font size as needed
ax.FontSize = fontSizeValue;

subplot(3,1,2)
plot(US2_i , 'r','LineWidth',4)
title("US - Punishing")
ylabel("Intensity")

% Get the current axes handle
ax = gca;
% Set the font size for all text elements
fontSizeValue = 14; % Adjust the font size as needed
ax.FontSize = fontSizeValue;

subplot(3,1,3)
plot(US1_i , 'g','LineWidth',4)
title("US - Rewarding")
ylabel("Intensity")
xlabel('Time (s)')

% Get the current axes handle
ax = gca;
% Set the font size for all text elements
fontSizeValue = 14; % Adjust the font size as needed
ax.FontSize = fontSizeValue;
%% 




