function [J] = get_J1(y,u,yd,OCP)
N = round(OCP.T /OCP.dt);
nu=OCP.nu1;
dt=OCP.dt;
gamma = OCP.gamma1;
alpha = OCP.alpha1;

%J=0.5*((y-yd)'*(y-yd)) + (nu/2)*u*(u')*dt + beta*sum(abs(u))*dt; 
% J = 0.5* gamma*sum((y(:,end)).^2) + 0.5*dt*nu*sum(sum(u.^2)) + 0.5*alpha*dt*sum(sum(y.^2)); 
J=0;
for i=1:N-1
J = J + 0.5*dt*u(:,i)'*nu*u(:,i) + 0.5*dt*y(:,i)'*alpha*y(:,i); 
end
J = J+ 0.5*(y(:,end)'*gamma*y(:,end));
end

