function [dpdt] = adjointDyn(A,B,p,y,yd,u)

dpdt=transpose(A)*p - (y-yd); 

end