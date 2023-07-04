function [J] = get_J2(y,u,yd,OCP)
N = round(OCP.T/OCP.dt);
nu=OCP.nu2;
dt=OCP.dt;
gamma = OCP.gamma2;
alpha = OCP.alpha2;
beta = OCP.beta2;

J=0;
for i=1:N-1
J = J + 0.5*dt*u(:,i)'*nu*u(:,i) + dt*beta*sum(abs(u(:,i))) ... 
    + 0.5*dt*(y(:,i)-yd(:,i))'*alpha*(y(:,i)-yd(:,i)) ;
end
J = J+ 0.5*((y(:,end)-yd(:,end))'*gamma*(y(:,end)-yd(:,end)));
end

