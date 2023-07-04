function [p] = backward2(A,B1,B2,y,yd,u1,u2,OCP)
dt=OCP.dt;
Nt=round(OCP.T/dt);
gamma= OCP.gamma2; 
alpha = OCP.alpha2;

N=max(size(yd));
p = zeros(N,Nt+1);

p(1:N,Nt+1)= - gamma*(y(1:N,end)-yd);
I=eye(N,N);


for n=Nt:-1:1
   
    p(1:N,n)=(I-0.5*dt*A')\((I+0.5*dt*A')* p(1:N,n+1)  -  0.5*dt*alpha *I* (y(1:N,n) +y(1:N,n+1) ));
   
    
end


end

