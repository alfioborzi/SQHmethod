function [J] = get_J1(y,u,yd,OCP)
N = round(OCP.T /OCP.dt);
nu=OCP.nu1;
dt=OCP.dt;
gamma = OCP.gamma1;
alpha = OCP.alpha1;
beta = OCP.beta1;


J=0;
for i=1:N-1
J = J + 0.5*dt*u(:,i)'*nu*u(:,i)  + dt*beta*sum(abs(u(:,i))) ...
    + 0.5*dt*(y(:,i)-yd(:,i))'*alpha*(y(:,i)-yd(:,i)) ;
end
J = J+ 0.5*((y(:,end)-yd(:,end))'*gamma*(y(:,end)-yd(:,end)));
end


