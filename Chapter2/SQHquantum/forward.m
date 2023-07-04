function [y] = forward(A,B,y0,u,OCP)
dt=OCP.dt;
Nt=round(OCP.T/dt);
N=max(size(y0));
y=zeros(N,Nt+1);
y(1:N,1)=y0;
I=eye(N,N);
for n=1:Nt
    C=2*A+(u(n)+u(n+1))*B;
    y(1:N,n+1)=(I-0.25*dt*C)\(y(1:N,n)+0.25*dt*(C)*y(1:N,n));
end
