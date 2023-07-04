function [p] = backward2(A,B1,B2,y,yd,u1,u2,OCP)
dt=OCP.dt;
Nt=round(OCP.T/dt);
gamma= OCP.gamma2; 
alpha = OCP.alpha2;

N=2;
p = zeros(N,Nt+1);
I=eye(N,N);

p(:,Nt+1)= - gamma*(y(:,end)-yd(:,end));

for n=Nt:-1:1
   
p(:,n)=(I-0.5*dt*A')\((I+0.5*dt*A')* p(:,n+1) ...
  -  0.5*dt*alpha*I * ( y(:,n) +y(:,n+1) - yd(:,n) - yd(:,n+1)) );
     
end

end

