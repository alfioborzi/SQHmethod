function [H] = HPfunction(y,p,u,OCP,DYN)

H = (p')*stateDyn(DYN,y,u);

end