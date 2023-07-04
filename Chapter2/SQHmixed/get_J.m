function [J] = get_J(y,u,ct,yd,OCP)

T =OCP.T;
Nt=OCP.Nt;
dt=OCP.T/OCP.Nt;

nu=OCP.nu;
beta=OCP.beta;
alpha=OCP.alpha;

yT = y(:,end);

sumj = 0;
for j=1:Nt
    sumj = sumj + (max(0, gc(y(:,j),u(:,j)) + ct(1,j)/alpha )).^2;
end

J=0.5*((yT-yd)'*(yT-yd)) ... 
    +(nu/2)*u(1,:)*u(1,:)'*dt+beta*sum(abs(u(1,:)))*dt ... 
    +(nu/2)*u(2,:)*u(2,:)'*dt+beta*sum(abs(u(2,:)))*dt ... 
    + 0.5 * alpha*sumj*dt;
end

