function [J] = get_J(y,u,yd,OCP)

T =OCP.T;
Nt=OCP.Nt;
dt=OCP.T/OCP.Nt;

nu=OCP.nu;
beta=OCP.beta;

yT = y(:,end);

J=0.5*((yT-yd)'*(yT-yd))+(nu/2)*u(1,:)*u(1,:)'*dt... 
    +beta*sum(abs(u(1,:)))*dt ...
    +(nu/2)*u(2,:)*u(2,:)'*dt+beta*sum(abs(u(2,:)))*dt ;
end

