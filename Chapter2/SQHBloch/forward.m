function [y] = forward(A,B1,B2,Dv,y0,u,OCP)

T =OCP.T;
Nt=OCP.Nt;
dt=OCP.T/OCP.Nt;

N=max(size(y0));
y=zeros(N,Nt+1);
y(1:N,1)=y0;
I=eye(N,N);

% Midpoint in the bilinear case 
for n=1:Nt
    C=2*A + (u(1,n)+u(1,n+1))*B1 + (u(2,n)+u(2,n+1))*B2;
    y(1:N,n+1)=(I-0.25*dt*C)\(y(1:N,n)+0.25*dt*(C)*y(1:N,n));
    y(1:N,n+1) = y(1:N,n+1)+ dt*Dv(1:N) ; % splitting
end

end
