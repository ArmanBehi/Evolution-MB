function   m2 = runga4m2(m1,m2,i,k2)

%a function to calculate the value of MBON1 cell during an specefic time
%period (one time step)
%

tspan = [i i+1];
y0 = m1;
m1 = ode45(@(t) - (3 * m2 * t) - ( 2*m1*t) + k2 , tspan, y0);