function [H] = HPfunction(A,B,y,p,u,OCP)

% H=-(OCP.nu/2)*u^2 - OCP.beta*abs(u) + (p')*(A+u*B)*y; 

u1 = u(1);


H=-(OCP.nu/2)*(u1^2) - OCP.beta*(abs(u1)) ... 
    + (p')*stateDyn(A,B,y,u); 

end