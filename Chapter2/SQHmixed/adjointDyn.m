function [dpdt] = adjointDyn(A,B,p,y,u)

dpdt=transpose(A+u*B)*p; 

end