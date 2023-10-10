clc;
clear;

alpha =  1;
beta  =  3.4;
w11 = 1;
w12 = 1;
w21 = 1;
w22 = 1;


A = [alpha,-beta; -beta,alpha]; %system matrix
B=[w11 0 w21 0; 0 w12 0 w22]; %constant input


fs = 13;    % fontsize for legends

%% solve for steady-state

% steady-state condition is 0 = dv/dt = A dot Xss + B or, equivalently, Xss = A^(-1) dot (-B)
Xss = A \ (-B);

%% obtain eigenvectors and eigenvalues
[V, D] = eig( A );

%% alternative variable names
I_1 = w11 + w21;
I_2 = w12 + w22;

% b1 = B(1,1);  %alternative variable names for constant input
% b2 = B(2,1);
% b3 = B(3,1);
% b4 = B(4,1);


m1ss = Xss(1,1);
m2ss = Xss(1,2);

lambda1  = D(1,1);  % eigenvalues
lambda2  = D(1,2);
lambda3  = D(2,1);
lambda4  = D(2,2);

E1  = V(1,1);       % eigenvectors
E2  = V(1,2);
E2  = V(2,1);
E2  = V(2,2);

statesize = 10;

X0 = [0 0]';

tend = 10;


% xmin = floor(xss - 0.5*statesize);
% xmax = ceil(xss + 0.5*statesize);
xmin=-10;
xmax=10;

%% nullclines dot x = 0 and dot y = 0
xi = linspace(xmin, xmax, 100);

m1_ss = -beta*xi + alpha*xi + I_1; %considering the input(k) constant
m2_ss = -( -beta*xi + alpha*xi + I_2); 


% ymin = floor(min([ydx0i ydy0i]));
% ymax = ceil(max([ydx0i ydy0i]));
ymin=-10;
ymax=10;

%% gradient vectors [dot x, dot y]'
nv = 10;

xv = linspace(xmin, xmax, nv );
yv = linspace(ymin, ymax, nv );

[XV, YV] = meshgrid( xv, yv );

Dm1V = ( alpha*XV - beta*YV + I_1);
Dm2V = (-beta*XV + alpha*YV+ I_2);

%% plot state space
figure;

hold on;

% plot( xi, k1dx0i, 'r', 'LineWidth', 2);   % nullclines
% plot( xi, k2dy0i, 'b', 'LineWidth', 2);
plot( xi, m1_ss, 'g', 'LineWidth', 2);   % nullclines
plot( xi, m2_ss, 'r', 'LineWidth', 2);

% set(h,'FontSize', fs );
%  
% quiver( XV, YV, Dk1V, Dk2V,'k', 'LineWidth', 1 );   % gradients
quiver( XV, YV, Dm1V, Dm2V,'b', 'LineWidth', 1 ); 

% plot( m1ss, m2ss, 'k+', 'MarkerSize', 10, 'LineWidth', 2 );   % steady-state
plot( 0, 0, 'ko', 'MarkerSize', 5, 'LineWidth', 1 );   % coordinate origin

plot( [xmin xmax], [0 0], 'k:', 'LineWidth', 2);    % coordinate axes
plot( [0 0], [ymin ymax], 'k:', 'LineWidth', 2);

hold off;
axis 'square';
axis([xmin xmax ymin ymax]);
xlabel( 'I_1', 'FontSize', fs );
ylabel( 'I_2', 'FontSize', fs );

h=legend( 'MBON1', 'MBON2', 'MBON1-MBON2 System','coordinate origin');
annotation('textbox', [0.3, 0.90, 0.1, 0.1], 'String', "alpha(self-excitation) = " + alpha +", beta(mutual-inhibition) ="+beta)
