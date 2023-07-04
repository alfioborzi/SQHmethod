function [dydt] = stateDyn(A,B,y,u)

u1 = u(1);

dydt = A*y + B*u1 ; 

end