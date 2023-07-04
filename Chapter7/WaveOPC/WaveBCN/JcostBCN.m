function J = JcostBCN(y,yd,u,alpha,beta,nu,gamma,Nx,Nt)
global L T v

dx = L/Nx; 
dt = T/Nt; 

J1 = 0;
for j=1:Nt-1
J1 = J1 + (y(j,:)-yd(j,:))*(y(j,:)-yd(j,:))';
end
J1 = 0.5*alpha*dx*dt*J1;

J2 = 0;
for j=1:Nt
J2 = J2 + (u(j,:))*(u(j,:))';
end
J2 = 0.5*nu*dt*J2;

J3 = 0;
for j=1:Nt
J3 = J3 + sum(abs(u(j,:)));
end
J3 = gamma*dt*J3;


J4 = 0.5*beta*dx*(y(end,:)-yd(end,:))*(y(end,:)-yd(end,:))';

J = J1 + J2 + J3 + J4 ;

end