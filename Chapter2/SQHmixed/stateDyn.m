function [dydt] = stateDyn(A,B1,B2,gv,y,u)

u1 = u(1);
u2 = u(2);

dydt=( A + u1*B1 + u2*B2 )*y + gv; 

end