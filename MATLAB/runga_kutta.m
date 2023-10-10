function M1 = estimate_MBON(initial,K) 
% Runge-Kutta(Order 4) Algorithm for MBON1


 f1 = @(t,m1,m2) (-3*m1 - 2*m2);
 a =    %'Enter left end ponit
 b = %'Enter right end point
 n = %'Enter no. of subintervals          % n=(b-a)/h
 alpha =0 %'Enter the initial condition, alpha
 
 h = (b-a)/n;
 t = a;
 w = alpha;



   k1 = h*(f(t,m1,m2) + k1(t));
   k2 = h*(f(t+h/2.0, w+k1/2.0) + k1(t+h));
   k3 = h*(f(t+h/2.0, w+k2/2.0) + k1(t+h/2.0));
   k4 = h*(f(t+h,w+k3) + k1(t+h));
   w = w+(k1+2.0*(k2+k3)+k4)/6.0;
   t = a+i*h;
   fprintf('%5.3f %11.7f\n', t, w);
   plot(t,w,'b*'); grid on; 
   xlabel('t values'); ylabel('w values');
   hold on;
 end