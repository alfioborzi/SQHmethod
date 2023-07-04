function [J] = get_J1(y,u,yd,OCP)
N = round(OCP.T /OCP.dt);
nu=OCP.nu1;
dt=OCP.dt;
gamma = OCP.gamma1;
alpha = OCP.alpha1;
beta = OCP.beta1;


J=0;
for i=1:N-1
J = J + 0.5*dt*u(:,i)'*nu*u(:,i) + 0.5*dt*y(:,i)'*alpha*y(:,i) ... 
      + dt*beta*sum(abs(u(:,i))) ; 
end
J = J+ 0.5*(y(:,end)'*gamma*y(:,end));
end

