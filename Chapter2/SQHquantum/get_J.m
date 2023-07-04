function [J] = get_J(y,u,yd,OCP)
nu=OCP.nu;
dt=OCP.dt;
s=OCP.s;
beta=OCP.beta;
J=0.5*(norm(y-yd)^2) + (nu/2)*u*(u')*dt + beta*sum((abs(u)>s).*abs(u))*dt;
end

