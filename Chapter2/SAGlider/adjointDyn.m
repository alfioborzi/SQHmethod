function [dpdt] = adjointDyn(DYN,p,y,u)


% Dynamics parameters
sigma=DYN.sigmaKC;
b=DYN.b;
alpha0=DYN.alpha0;
K=DYN.K;


% y = x,y,v,theta
% u = alpha, eta 

Cx = 1 - cos(2*alpha0)*cos(2*u(1));
Cy = K * sin(2*alpha0)*sin(2*u(1));

dpdt(1) = 0 ;
dpdt(2) = 0 ;
dpdt(3) = p(1)*cos(y(4)) + p(2)*sin(y(4)) ... 
  - 2*sigma*p(3)*y(3)*Cx*(1 + u(2)*b ) + sigma* p(4)*Cy*(1 + u(2)*b )... 
  + p(4)*cos(y(4))/y(3)^2 ;
dpdt(4) = -p(1)*y(3)*sin(y(4)) + p(2)*y(3)*cos(y(4))- p(3)*cos(y(4)) ... 
    +p(4)*sin(y(4))/y(3);

 dpdt = dpdt' ;

end