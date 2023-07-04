function [y] = forward(A,B1, B2 ,y0,u1, u2,OCP)
dt=OCP.dt;
Nt=round(OCP.T/dt);
N=max(size(y0));

y=zeros(N,Nt+1);

y(1:N,1)=y0;
I=eye(N,N);

for n=1:Nt
    C=0.5* ( B1*( u1(1:N,n)+u1(1:N,n+1) )  + B2*( u2(1:N,n)+u2(1:N,n+1) )  );
    y(1:N,n+1)= (I-0.5*dt*A) \ ( (I + 0.5*dt*A)*y(1:N,n) + dt*C );
 
end
