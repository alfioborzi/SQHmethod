function [J] = get_J(y,u,yd,OCP)

T =OCP.T;
Nt=OCP.Nt;
dt=OCP.T/OCP.Nt;

nu=OCP.nu;
beta=OCP.beta;

J=0.5*((y-yd)*(y-yd)')*dt+(nu/2)*u(1,:)*u(1,:)'*dt ;
end

