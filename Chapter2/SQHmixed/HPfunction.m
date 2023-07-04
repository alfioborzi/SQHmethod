function [H] = HPfunction(A,B1,B2,gv,y,p,u,ct,OCP)

% H=-(OCP.nu/2)*u^2 - OCP.beta*abs(u) + (p')*(A+u*B)*y; 

u1 = u(1);
u2 = u(2);

alpha=OCP.alpha;

H = -(OCP.nu/2)*(u1^2 + u2^2) - OCP.beta*(abs(u1)+abs(u2)) ... 
    - 0.5 * alpha * (max(0, gc(y,u) + ct/alpha ) ).^2 ...
    + (p')*stateDyn(A,B1,B2,gv,y,u); 
end