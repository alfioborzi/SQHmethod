function [dydt] = stateDyn(DYN,y,u)

% Dynamics parameters
sigma=DYN.sigmaKC;
b=DYN.b;
alpha0=DYN.alpha0;
K=DYN.K;

% y = x,y,v,theta
% u = alpha, eta 

Cx = 1 - cos(2*alpha0)*cos(2*u(1));
Cy = K * sin(2*alpha0)*sin(2*u(1));

dydt(1) = y(3)*cos(y(4));
dydt(2) = y(3)*sin(y(4));
dydt(3) = -sigma * y(3)^2 * Cx * (1 + u(2)*b ) - sin(y(4));
dydt(4) =  sigma * y(3)   * Cy * (1 + u(2)*b ) - cos(y(4))/y(3);

dydt=dydt' ;

end